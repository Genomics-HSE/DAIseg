import numpy as np
from numba import jit, prange, config, set_num_threads
import json
import hmm
import gc
import sys, os
import pandas as pd

# restruc numba 
if 'NUMBA_NUM_THREADS' not in os.environ:
    os.environ['NUMBA_NUM_THREADS'] = '32'
    print(f"NUMBA_NUM_THREADS={os.environ['NUMBA_NUM_THREADS']}")

# NORMALIZED FORWARD-BACKWARD 
@jit(nopython=True)
def forward_backward_normalized(emit, trans, start):
    """
    Forward-Backward using scaling variables.
    """
    N, n_states = emit.shape

    alpha = np.zeros((N, n_states))
    scales = np.zeros(N) 

    for s in range(n_states):
        alpha[0, s] = start[s] * emit[0, s]

    scales[0] = 1.0 / (np.sum(alpha[0]) + 1e-300)
    alpha[0] *= scales[0]

    for t in range(1, N):
        for s in range(n_states):
            acc = 0.0
            for p in range(n_states):
                acc += alpha[t-1, p] * trans[p, s]
            alpha[t, s] = acc * emit[t, s]

        scales[t] = 1.0 / (np.sum(alpha[t]) + 1e-300)
        alpha[t] *= scales[t]

    log_lik = -np.sum(np.log(scales + 1e-300))

    beta = np.zeros((N, n_states))
    beta[N-1, :] = scales[N-1]

    for t in range(N-2, -1, -1):
        for s in range(n_states):
            acc = 0.0
            for next_s in range(n_states):
                acc += trans[s, next_s] * emit[t+1, next_s] * beta[t+1, next_s]
            beta[t, s] = acc * scales[t]

    gamma = alpha * beta
    for t in range(N):
        gamma[t] /= (np.sum(gamma[t]) + 1e-300)

    return gamma, log_lik

# E-STEP 
@jit(nopython=True, parallel=True)
def e_step_normalized(emit, trans, start, O1, O2, L1, L2):
    M, N, n_states = emit.shape
    numerators = np.zeros((M, 3)) 
    denominators = np.zeros((M, 3))
    total_log_lik = 0.0

    for m in prange(M):
        gamma, log_lik = forward_backward_normalized(emit[m], trans, start)
        total_log_lik += log_lik

        num_n = np.sum(gamma[:, 0] * O2[m]) + np.sum(gamma[:, 1] * O1[m])
        den_n = np.sum(gamma[:, 0] * L2[m]) + np.sum(gamma[:, 1] * L1[m])

        num_af = np.sum(gamma[:, 0] * O1[m])
        den_af = np.sum(gamma[:, 0] * L1[m])

        num_i = np.sum(gamma[:, 1] * O2[m])
        den_i = np.sum(gamma[:, 1] * L2[m])

        numerators[m, 0] = num_n
        numerators[m, 1] = num_af
        numerators[m, 2] = num_i

        denominators[m, 0] = den_n
        denominators[m, 1] = den_af
        denominators[m, 2] = den_i

    return numerators, denominators, total_log_lik

# MAIN EM 
def train_em_normalized(O1, O2, L1, L2, init_lmbd, rr, Ti, init_a, max_iter=20, tol=1e-8):
    M, N = O1.shape
    curr_lmbd = np.array(init_lmbd, dtype=float)
    prev_log_lik = -np.inf

    print(f"Starting EM... Max Iter: {max_iter}")

    for it in range(max_iter):

        log_emissions = hmm.compute_emissions_custom(O1, O2, L1, L2, curr_lmbd)
        log_A = hmm.get_log_A(1000, rr, Ti, init_a)
        log_start = np.log([1.0 - init_a, init_a])

        emit_linear = np.exp(log_emissions)
        trans_linear = np.exp(log_A)
        start_linear = np.exp(log_start)

        nums, dens, log_lik = e_step_normalized(
            emit_linear, trans_linear, start_linear, O1, O2, L1, L2
        )

        total_nums = np.sum(nums, axis=0)
        total_dens = np.sum(dens, axis=0)

        # Update Rates
        new_lmbd = total_nums / (total_dens + 1e-10)

        diff = log_lik - prev_log_lik
        print(f"  Iter {it+1}: LL={log_lik:.2f}, Delta={diff:.4f} | Rates: N={new_lmbd[0]:.5f}, AF={new_lmbd[1]:.5f}, I={new_lmbd[2]:.5f}")

        if abs(diff) < tol and it > 0:
            print("  Converged.")
            break

        prev_log_lik = log_lik
        curr_lmbd = new_lmbd

    return curr_lmbd

def run_batch_em_pipeline(json_files_list, output_combined_file=None, max_iter=15, tol=1e-6):
    """
    Loads ALL chromosomes into memory .
    Runs EM optimization.
    Runs Viterbi on EACH chromosome
    """

    # Storage
    batch_data = []

    print(f"Loading {len(json_files_list)} files for Global EM & Inference...")

    init_params = None
    rr_val = 0

    for j_file in json_files_list:
        with open(j_file, 'r') as f:
            data = json.load(f)

        if init_params is None:
            prms = data["parameters_initial"]
            gen_time, mu, rr, l = prms['generation_time'], prms['mutation'], prms['rr'], prms['window_length']
            d = mu * l / gen_time
            Ti_generations = prms['t_introgression'] / gen_time
            init_params = [d*prms['t_archaic_c'], d*prms['t_split_c'], d*prms['t_introgression_c'], Ti_generations, prms['admixture_proportion']]
            rr_val = rr
        # Load Observations
        try:
            obs_seq, _ = hmm.create_observations(
                os.path.join(data["prefix"], data["data"]), 
                os.path.join(data["prefix"], data["window_callability"]["Thousand_genomes"])
            )
        except Exception as e:
            print(f"Critical Error loading observations for {j_file}: {e}")
            sys.exit(1)

        # Load 1kG Mask
        path_1kg = os.path.join(data['prefix'], data["window_callability"]["Thousand_genomes"])
        if not os.path.exists(path_1kg):
            print(f"Error: 1kG mask missing: {path_1kg}")
            sys.exit(1)

        cal_1kG = np.loadtxt(path_1kg, usecols=-1)

        # Load ND Mask 
        path_nd = os.path.join(data['prefix'], data.get("window_callability", {}).get("Nd_1k_genomes", ""))

        # Determine sequence length N from obs or 1kg mask
        O1_temp, _, _ = hmm.prepare_matrices_from_dict(obs_seq)
        M_check, N_check = O1_temp.shape

        if path_nd and os.path.exists(path_nd):
            try:
                cal_nd_1kG = np.loadtxt(path_nd, usecols=-1)
                # Resize if length mismatch (e.g. empty file or wrong chrm)
                if len(cal_nd_1kG) != N_check:
                    print(f"Warning: ND mask length mismatch in {j_file}. Resizing/Zeroing.")
                    cal_nd_1kG = np.zeros(N_check)
            except:
                cal_nd_1kG = np.zeros(N_check)
        else:

            cal_nd_1kG = np.zeros(N_check)

        # Prepare Matrices
        O1, O2, names = hmm.prepare_matrices_from_dict(obs_seq)
        M, N = O1.shape

        L1_mat = np.tile(cal_1kG, (M, 1))
        L2_mat = np.tile(cal_nd_1kG, (M, 1))

        # Если покрытие 0, наблюдения ОБЯЗАНЫ быть 0.
        # Для Неандертальцев (L2)
        zero_mask_nd = (L2_mat <= 1e-9)
        if np.any(zero_mask_nd):
            O2[zero_mask_nd] = 0

        # Для Людей (L1) - на всякий случай
        zero_mask_1kg = (L1_mat <= 1e-9)
        if np.any(zero_mask_1kg):
            O1[zero_mask_1kg] = 0

        batch_data.append({
            "em": (O1, O2, L1_mat, L2_mat),
            "viterbi": (cal_1kG, cal_nd_1kG, names),
            "json": data
        })

        gc.collect()

    print("Data loaded. Starting Global EM optimization...")

    # GLOBAL TRAINING

    curr_lmbd = np.array(init_params[:3], dtype=float)
    Ti = init_params[3]
    a_param = init_params[4]

    prev_log_lik = -np.inf

    for it in range(max_iter):
        total_nums = np.zeros(3)
        total_dens = np.zeros(3)
        iter_log_lik = 0.0

        log_A = hmm.get_log_A(1000, rr_val, Ti, a_param)
        log_start = np.log([1.0 - a_param, a_param])

        trans_linear = np.exp(log_A)
        start_linear = np.exp(log_start)

        for item in batch_data:
            O1, O2, L1, L2 = item["em"]

            log_emissions = hmm.compute_emissions_custom(O1, O2, L1, L2, curr_lmbd)
            emit_linear = np.exp(log_emissions)

            nums, dens, ll = e_step_normalized(emit_linear, trans_linear, start_linear, O1, O2, L1, L2)

            total_nums += np.sum(nums, axis=0)
            total_dens += np.sum(dens, axis=0)
            iter_log_lik += ll


        new_lmbd = curr_lmbd.copy()
        for k in range(3):
            # Обновляем только если есть данные (покрытие > 0)
            if total_dens[k] > 1e-6:
                new_lmbd[k] = total_nums[k] / total_dens[k]
            else:
                pass 

        diff = iter_log_lik - prev_log_lik
        print(f"Iter {it+1}: LL={iter_log_lik:.2f} | Rates: N={new_lmbd[0]:.5f}, AF={new_lmbd[1]:.5f}, I={new_lmbd[2]:.5f}")

        if abs(diff) < tol and it > 0:
            print("Converged.")
            break

        prev_log_lik = iter_log_lik
        curr_lmbd = new_lmbd

    final_params = np.concatenate([curr_lmbd, [Ti, a_param]])
    print(f"Final Global Parameters: {final_params}")

    # INFERENCE

    print("Starting Inference on all files...")

    all_results = []

    for item in batch_data:
        O1, O2, _, _ = item["em"]
        cal_1kG, cal_nd_1kG, names = item["viterbi"]
        jsn_data = item["json"]

        result = hmm.run_hmm(O1, O2, cal_1kG, cal_nd_1kG, lmbd=final_params, rr=rr_val)

        dictionary = {k: v for k, v in zip(names, result)}
        out_dict = {name: hmm.get_tracts(dictionary[name]) for name in names}

        # clean gaps
        out_dict_new = hmm.clean_gaps(out_dict, jsn_data.get("gaps"), jsn_data["CHROM"])

        # Saving results
        output_tsv = f"{jsn_data['prefix']}/{jsn_data['output']}.em.tsv"
        with open(output_tsv, "w", encoding="utf-8") as f:
            f.write("Sample\tCHROM\tStart\tEnd\tLength\n")
            for sample_name, tracks in out_dict_new.items():
                for start, end in tracks.get('Archaic', []):
                    f.write(f"{sample_name}\t{jsn_data['CHROM']}\t{start}\t{end}\t{end - start + 1}\n")
                    if output_combined_file:
                        all_results.append({
                            "CHR": jsn_data["CHROM"],
                            "Sample": sample_name,
                            "Start": start,
                            "End": end,
                            "Length": end - start + 1
                        })

    print(f"Done! Processed {len(batch_data)} files.")

    # SAVING MERGED FILE 
    if output_combined_file:
        print(f"Saving merged results to {output_combined_file}...")
        df = pd.DataFrame(all_results, columns=[ "Sample", "CHR",  "Start", "End", "Length"])

        if not df.empty:
            try:
                df['CHR'] = df['CHR'].astype(int)
                df = df.sort_values(by=['CHR', 'Sample', 'Start'])
            except:
                pass
        df.to_csv(output_combined_file, sep='\t', index=False)





@jit(nopython=True)
def forward_backward_xi_normalized(emit, trans, start):
    """
    Выполняет Forward-Backward и возвращает:
    Gamma, Xi_sum (для обновления переходов T и a), Log-Likelihood
    """
    N, n_states = emit.shape

    #  Forward
    alpha = np.zeros((N, n_states))
    scales = np.zeros(N)

    # t=0
    for s in range(n_states):
        alpha[0, s] = start[s] * emit[0, s]

    scales[0] = 1.0 / (np.sum(alpha[0]) + 1e-300)
    alpha[0] *= scales[0]

    # t > 0
    for t in range(1, N):
        for s in range(n_states):
            acc = 0.0
            for p in range(n_states):
                acc += alpha[t-1, p] * trans[p, s]
            alpha[t, s] = acc * emit[t, s]

        scales[t] = 1.0 / (np.sum(alpha[t]) + 1e-300)
        alpha[t] *= scales[t]

    log_lik = -np.sum(np.log(scales + 1e-300))

    # Backward Pass
    beta = np.zeros((N, n_states))
    beta[N-1, :] = scales[N-1]

    for t in range(N-2, -1, -1):
        for s in range(n_states):
            acc = 0.0
            for next_s in range(n_states):
                acc += trans[s, next_s] * emit[t+1, next_s] * beta[t+1, next_s]
            beta[t, s] = acc * scales[t]

    #  Gamma & Xi 
    gamma = np.zeros((N, n_states))
    xi_sum = np.zeros((2, 2)) 

    # Gamma
    for t in range(N):
        denom = 0.0
        for s in range(n_states):
            val = alpha[t, s] * beta[t, s]
            gamma[t, s] = val
            denom += val
        gamma[t] /= (denom + 1e-300)

    # Xi 
    for t in range(N - 1):
        current_xi_sum = 0.0
        temp_xi = np.zeros((2, 2))

        for i in range(n_states):
            for j in range(n_states):
                val = alpha[t, i] * trans[i, j] * emit[t+1, j] * beta[t+1, j]
                temp_xi[i, j] = val
                current_xi_sum += val

        # Нормализуем, чтобы сумма по всем переходам в момент t была равна 1
        inv_sum = 1.0 / (current_xi_sum + 1e-300)
        for i in range(n_states):
            for j in range(n_states):
                xi_sum[i, j] += temp_xi[i, j] * inv_sum

    return gamma, xi_sum, log_lik

@jit(nopython=True, parallel=True)
def e_step_v2_full(emit, trans, start, O1, O2, L1, L2):
    """
    Возвращает для Lambda - (nums, dens) И статистику для матрицы A (xi).
    """
    M, N, n_states = emit.shape

    # Для Lambda
    numerators = np.zeros((M, 3))
    denominators = np.zeros((M, 3))
    xi_accum = np.zeros((M, 2, 2))
    total_log_lik = 0.0

    for m in prange(M):
        gamma, xi, log_lik = forward_backward_xi_normalized(emit[m], trans, start)
        total_log_lik += log_lik

        # Сохраняем Xi
        xi_accum[m] = xi

        numerators[m, 0] = np.sum(gamma[:, 0] * O2[m]) + np.sum(gamma[:, 1] * O1[m])
        denominators[m, 0] = np.sum(gamma[:, 0] * L2[m]) + np.sum(gamma[:, 1] * L1[m])

        # African (Rate 1)
        numerators[m, 1] = np.sum(gamma[:, 0] * O1[m])
        denominators[m, 1] = np.sum(gamma[:, 0] * L1[m])

        # Introgression (Rate 2)
        numerators[m, 2] = np.sum(gamma[:, 1] * O2[m])
        denominators[m, 2] = np.sum(gamma[:, 1] * L2[m])

    return numerators, denominators, xi_accum, total_log_lik

def update_structural_params_logic(xi_sum, r, L, current_T, current_a):
    """
    Вычисляет новые T и a 
    """
    # Эмпирические вероятности переходов
    n00, n01 = xi_sum[0, 0], xi_sum[0, 1]
    n10, n11 = xi_sum[1, 0], xi_sum[1, 1]

    denom0 = n00 + n01
    p01 = n01 / denom0 if denom0 > 0 else 0.0

    denom1 = n10 + n11
    p10 = n10 / denom1 if denom1 > 0 else 0.0

    # Оценка rho (общая вероятность рекомбинации)
    rho_est = p01 + p10

    # Защита
    if rho_est <= 1e-15 or rho_est >= 1.0:
        return current_T, current_a

    # 3. Восстановление T
    # rho = 1 - exp(-T * r * L)
    try:
        new_T = -np.log(1.0 - rho_est) / (r * L)
    except:
        new_T = current_T

    new_a = p01 / rho_est

    # Constraints
    new_T = max(new_T, 100.0)
    new_a = np.clip(new_a, 0.0001, 0.999)

    return new_T, new_a






# PIPELINE V2: Updates Rates AND Transition Matrix
def run_batch_em_pipeline_v2(json_files_list, output_combined_file=None, max_iter=15, tol=1e-6, threads=4):
    """
    V2 Pipeline :
    Updates Emission Rates (Lambdas) 
    Updates Transition Matrix (T, a)
    Handles missing data (e.g. no Neanderthal) gracefully.
    """

    # 1. Load Data
    batch_data = []
    print(f"[V2] Loading {len(json_files_list)} files...")

    init_params_full = None # [rates... , T, a]
    rr_val = 0
    L_val = 0

    for j_file in json_files_list:
        with open(j_file, 'r') as f:
            data = json.load(f)

        if init_params_full is None:
            prms = data["parameters_initial"]
            gen_time = prms['generation_time']
            mu = prms['mutation']
            rr_val = prms['rr']
            L_val = prms['window_length']

            # Initial Rates
            d = mu * L_val / gen_time
            lmbd_init = [d*prms['t_archaic_c'], d*prms['t_split_c'], d*prms['t_introgression_c']]

            T_init = prms['t_introgression'] / gen_time # generations
            a_init = prms['admixture_proportion']

            init_params_full = {
                'rates': np.array(lmbd_init, dtype=float),
                'T': T_init,
                'a': a_init
            }

        # Load Observations
        try:
            obs_seq, _ = hmm.create_observations(
                os.path.join(data["prefix"], data["data"]), 
                os.path.join(data["prefix"], data["window_callability"]["Thousand_genomes"])
            )
        except Exception as e:
            print(f"Skipping {j_file} due to load error: {e}")
            continue

        if not obs_seq: continue


        # 1kG Mask
        path_1kg = os.path.join(data['prefix'], data["window_callability"]["Thousand_genomes"])
        if not os.path.exists(path_1kg):
            print(f"Error: 1kG mask missing for {j_file}")
            continue
        cal_1kG = np.loadtxt(path_1kg, usecols=-1)

        # Determine sequence length
        O1_temp, _, _ = hmm.prepare_matrices_from_dict(obs_seq)
        M_check, N_check = O1_temp.shape

        # ND Mask
        path_nd = os.path.join(data['prefix'], data.get("window_callability", {}).get("Nd_1k_genomes", ""))

        if path_nd and os.path.exists(path_nd):
            try:
                cal_nd_1kG = np.loadtxt(path_nd, usecols=-1)
                if len(cal_nd_1kG) != N_check:
                    print(f"Warning: ND mask length mismatch in {j_file}. Using Zeros.")
                    cal_nd_1kG = np.zeros(N_check)
            except:
                cal_nd_1kG = np.zeros(N_check)
        else:
            # Missing file or key -> Zero coverage
            cal_nd_1kG = np.zeros(N_check)

        O1, O2, names = hmm.prepare_matrices_from_dict(obs_seq)
        M, N = O1.shape

        L1_mat = np.tile(cal_1kG, (M, 1))
        L2_mat = np.tile(cal_nd_1kG, (M, 1))

        # If L=0, force O=0.
        zero_mask_nd = (L2_mat <= 1e-9)
        if np.any(zero_mask_nd):
            O2[zero_mask_nd] = 0

        zero_mask_1kg = (L1_mat <= 1e-9)
        if np.any(zero_mask_1kg):
            O1[zero_mask_1kg] = 0

        batch_data.append({
            "em": (O1, O2, L1_mat, L2_mat),
            "viterbi": (cal_1kG, cal_nd_1kG, names),
            "json": data
        })
        gc.collect()

    if not batch_data:
        print("No valid data loaded. Exiting.")
        return

    # EM Loop
    print(f"[V2] Starting EM Optimization (Rates + T + a)...")

    curr_rates = init_params_full['rates']
    curr_T = init_params_full['T']
    curr_a = init_params_full['a']

    prev_log_lik = -np.inf

    for it in range(max_iter):
        total_nums = np.zeros(3)
        total_dens = np.zeros(3)
        total_xi = np.zeros((2, 2))
        iter_log_lik = 0.0

        # Update Transition Matrix A
        full_lambda_vec = np.concatenate([curr_rates, [curr_T, curr_a]])

        log_A = hmm.get_log_A(L_val, rr_val, curr_T, curr_a)
        log_start = np.log(np.array([1.0 - curr_a, curr_a]) + 1e-300)

        trans_linear = np.exp(log_A)
        start_linear = np.exp(log_start)

        # E-Step over all batches
        for item in batch_data:
            O1, O2, L1, L2 = item["em"]

            # Recalculate Emissions B
            log_emissions = hmm.compute_emissions_custom(O1, O2, L1, L2, full_lambda_vec)
            emit_linear = np.exp(log_emissions)

            # Run V2 E-step
            nums, dens, xi_accum, ll = e_step_v2_full(
                emit_linear, trans_linear, start_linear, O1, O2, L1, L2
            )

            total_nums += np.sum(nums, axis=0)
            total_dens += np.sum(dens, axis=0)
            total_xi += np.sum(xi_accum, axis=0)
            iter_log_lik += ll

        # M-Step 

        # Update Rates 
        new_rates = curr_rates.copy()
        for k in range(3):
            if total_dens[k] > 1e-6:
                new_rates[k] = total_nums[k] / total_dens[k]
            else:
                pass # Keep old value

        #  Update T and a 
        new_T, new_a = update_structural_params_logic(total_xi, rr_val, L_val, curr_T, curr_a)

        # Check convergence
        diff_ll = iter_log_lik - prev_log_lik
        print(f"Iter {it+1}: LL={iter_log_lik:.2f}, dLL={diff_ll:.2f}")

        if abs(diff_ll) < tol and it > 0:
            print("Converged.")
            break

        prev_log_lik = iter_log_lik
        curr_rates = new_rates
        curr_T = new_T
        curr_a = new_a

    # RESULT OUTPUT
    print("\n" + "="*50)
    print("       EM OPTIMIZATION COMPLETED")
    print("="*50)
    print(f"Emission Rates (Poisson lambda per window):")
    print(f"  > Neutral (Neanderthal): {curr_rates[0]:.6f}")
    print(f"  > African (Background):  {curr_rates[1]:.6f}")
    print(f"  > Introgressed:          {curr_rates[2]:.6f}")
    print("-" * 50)
    print(f"Structural Parameters:")
    print(f"  > Time since Admixture: {curr_T:.1f} generations")
    print(f"  > Admixture Proportion: {curr_a:.4%} ({curr_a:.6f})")
    print("="*50 + "\n")

    # Inference & Save
    final_full_params = np.concatenate([curr_rates, [curr_T, curr_a]])

    all_results = []

    for item in batch_data:
        O1, O2, _, _ = item["em"]
        cal_1kG, cal_nd_1kG, names = item["viterbi"]
        jsn = item["json"]

        result = hmm.run_hmm(O1, O2, cal_1kG, cal_nd_1kG, lmbd=final_full_params, rr=rr_val)

        dictionary = {k: v for k, v in zip(names, result)}
        out_dict = {name: hmm.get_tracts(dictionary[name]) for name in names}
        out_dict_new = hmm.clean_gaps(out_dict, jsn.get("gaps"), jsn["CHROM"])

        # Save individual
        out_path = f"{jsn['prefix']}/{jsn['output']}.em_v2.tsv"
        with open(out_path, "w") as f:
            f.write("Sample\tCHROM\tStart\tEnd\tLength\n")
            for sample_name, tracks in out_dict_new.items():
                for start, end in tracks.get('Archaic', []):
                    f.write(f"{sample_name}\t{jsn['CHROM']}\t{start}\t{end}\t{end - start + 1}\n")
                    if output_combined_file:
                        all_results.append({
                            "CHR": jsn["CHROM"],
                            "Sample": sample_name,
                            "Start": start,
                            "End": end,
                            "Length": end - start + 1
                        })

    if output_combined_file:
        print(f"Saving merged results to {output_combined_file}")
        df = pd.DataFrame(all_results, columns=["Sample", "CHR", "Start", "End", "Length"])
        if not df.empty:
            try:
                df['CHR'] = pd.to_numeric(df['CHR'], errors='coerce')
                df = df.sort_values(by=['CHR', 'Sample', 'Start'])
            except: pass
        df.to_csv(output_combined_file, sep='\t', index=False)

    # Save optimized parameters 
    if output_combined_file:
        meta_file = output_combined_file + ".meta.json"
        meta_data = {
            "rate_neutral": curr_rates[0],
            "rate_african": curr_rates[1],
            "rate_introgressed": curr_rates[2],
            "T_introgression": curr_T,
            "admixture_prop": curr_a
        }
        with open(meta_file, 'w') as f:
            json.dump(meta_data, f, indent=4)
