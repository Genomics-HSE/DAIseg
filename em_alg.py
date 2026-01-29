import numpy as np
from numba import jit, prange, config, set_num_threads
import json
import hmm  
import gc
import sys,os

import pandas as pd

# ==========================================
# ОГРАНИЧЕНИЕ ПОТОКОВ NUMBA 
# ==========================================
if 'NUMBA_NUM_THREADS' not in os.environ:
    os.environ['NUMBA_NUM_THREADS'] = '32'
    print(f"Установлено NUMBA_NUM_THREADS={os.environ['NUMBA_NUM_THREADS']}")
else:
    print(f"NUMBA_NUM_THREADS уже установлено: {os.environ['NUMBA_NUM_THREADS']}")

# ==========================================
# 1. NORMALIZED FORWARD-BACKWARD 
# ==========================================

@jit(nopython=True)
def forward_backward_normalized(emit, trans, start):
    """
    Performs Forward-Backward using the Scaling Method.
    Inputs are in LINEAR space (probabilities), NOT logarithms.
    """
    N, n_states = emit.shape
    
    # --- FORWARD PASS ---
    alpha = np.zeros((N, n_states))
    scales = np.zeros(N) 
    
    # Init (t=0)
    for s in range(n_states):
        alpha[0, s] = start[s] * emit[0, s]
        
    scales[0] = 1.0 / (np.sum(alpha[0]) + 1e-300)
    alpha[0] *= scales[0]
        
    # Induction
    for t in range(1, N):
        for s in range(n_states):
            acc = 0.0
            for p in range(n_states):
                acc += alpha[t-1, p] * trans[p, s]
            alpha[t, s] = acc * emit[t, s]
            
        scales[t] = 1.0 / (np.sum(alpha[t]) + 1e-300)
        alpha[t] *= scales[t]

    log_lik = -np.sum(np.log(scales + 1e-300))
        
    # --- BACKWARD PASS ---
    beta = np.zeros((N, n_states))
    beta[N-1, :] = scales[N-1]

    for t in range(N-2, -1, -1):
        for s in range(n_states):
            acc = 0.0
            for next_s in range(n_states):
                acc += trans[s, next_s] * emit[t+1, next_s] * beta[t+1, next_s]
            beta[t, s] = acc * scales[t] 

    # --- GAMMA ---
    gamma = alpha * beta
    for t in range(N):
        gamma[t] /= (np.sum(gamma[t]) + 1e-300)
            
    return gamma, log_lik

# ==========================================
# 2. E-STEP (Parallel Statistics Collection)
# ==========================================

@jit(nopython=True, parallel=True)
def e_step_normalized(emit, trans, start, O1, O2, L1, L2):
    M, N, n_states = emit.shape
    numerators = np.zeros((M, 3)) 
    denominators = np.zeros((M, 3))
    total_log_lik = 0.0
    
    for m in prange(M):
        gamma, log_lik = forward_backward_normalized(emit[m], trans, start)
        total_log_lik += log_lik
        
        # 1. Lambda Neutral (lmbd[0])
        num_n = np.sum(gamma[:, 0] * O2[m]) + np.sum(gamma[:, 1] * O1[m])
        den_n = np.sum(gamma[:, 0] * L2[m]) + np.sum(gamma[:, 1] * L1[m])
        
        # 2. Lambda AF (lmbd[1])
        num_af = np.sum(gamma[:, 0] * O1[m])
        den_af = np.sum(gamma[:, 0] * L1[m])
        
        # 3. Lambda Intro (lmbd[2])
        num_i = np.sum(gamma[:, 1] * O2[m])
        den_i = np.sum(gamma[:, 1] * L2[m])
        
        # --- ИСПРАВЛЕНИЕ ЗДЕСЬ ---
        # Numba не любит присваивание списков срезам (numerators[m,:] = list).
        # Делаем присваивание поэлементно:
        
        numerators[m, 0] = num_n
        numerators[m, 1] = num_af
        numerators[m, 2] = num_i
        
        denominators[m, 0] = den_n
        denominators[m, 1] = den_af
        denominators[m, 2] = den_i
        
    return numerators, denominators, total_log_lik

# ==========================================
# 3. MAIN EM TRAINING LOOP
# ==========================================

def train_em_normalized(O1, O2, L1, L2, init_lmbd, rr, Ti, init_a, max_iter=20, tol=1e-4):
    M, N = O1.shape
    curr_lmbd = np.array(init_lmbd, dtype=float)
    prev_log_lik = -np.inf
    
    print(f"Starting EM... Max Iter: {max_iter}")
    
    for it in range(max_iter):
        # Use helper functions from hmm.py
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

# ==========================================
# 4. PIPELINE WRAPPER
# ==========================================

def run_daiseg_em(json_file):
    """
    Loads data using hmm.py utils, runs EM training, then runs Viterbi.
    """
    with open(json_file, 'r') as f:
        data = json.load(f)
        
    print(f"Loading data from {json_file}...")
    
    # 1. Prepare Data using hmm.py
    obs_seq, nst = hmm.create_observations(data["data"], data["window_callability"]["Thousand_genomes"])
    prms = data["parameters_initial"]
    gen_time, mu, rr, l = prms['generation_time'], prms['mutation'], prms['rr'], prms['window_length']
    
    # Initial Params
    d = mu * l / gen_time
    lambda_0 = [d*prms['t_archaic_c'],  d*prms['t_split_c'], d*prms['t_introgression_c'], d*prms['t_introgression'], prms['admixture_proportion']]
    
    cal_1kG = np.loadtxt(data["window_callability"]["Thousand_genomes"], usecols=-1)
    cal_nd_1kG = np.loadtxt(data["window_callability"]["Nd_1k_genomes"], usecols=-1)  
    
    O1, O2, names = hmm.prepare_matrices_from_dict(obs_seq)   
    M, N = O1.shape
    
    # 2. Expand Masks for EM (needs MxN matrices)
    L1_mat = np.tile(cal_1kG, (M, 1))
    L2_mat = np.tile(cal_nd_1kG, (M, 1))
    
    # 3. Run EM Optimization
    init_rates = lambda_0[:3] # [N, AF, I]
    Ti_param = lambda_0[3]
    a_param = lambda_0[4]
    
    optimized_rates = train_em_normalized(
        O1, O2, L1_mat, L2_mat, 
        init_lmbd=init_rates, 
        rr=rr, 
        Ti=Ti_param, 
        init_a=a_param, 
        max_iter=20
    )
    
    final_lambda = np.concatenate([optimized_rates, [Ti_param, a_param]])
    
    # 4. Run Standard Inference (Viterbi) with NEW parameters via hmm.py
    print("Running Viterbi with optimized parameters...")
    result = hmm.run_hmm(O1, O2, cal_1kG, cal_nd_1kG, lmbd=final_lambda, rr=rr)
    
    # 5. Save results
    dictionary = {k: v for k, v in zip(names, result)}
    out_dict = {}
    for name in names:
        out_dict[name] = hmm.get_tracts(dictionary[name])
        
    out_dict_new = hmm.clean_gaps(out_dict, data["gaps"], data["CHROM"])

    output_tsv = f"{data['prefix']}.em.tsv"
    print(f"Saving EM-optimized TSV results to: {output_tsv}")
    
    with open(output_tsv, "w", encoding="utf-8") as f:
        f.write("Sample\tStart\tEnd\n")
        for sample_name, tracks in out_dict_new.items():
            archaic_intervals = tracks.get('Archaic', [])
            for start, end in archaic_intervals:
                f.write(f"{sample_name}\t{start}\t{end}\n")
              
    return result

def run_batch_em_pipeline(json_files_list, output_combined_file=None, max_iter=15, tol=1e-6):
    """
    1. Loads ALL chromosomes into memory.
    2. Runs GLOBAL EM optimization.
    3. Runs Viterbi on EACH chromosome using the globally optimized parameters.
    4. Saves individual TSV files.
    5. Saves COMBINED TSV file (if requested).
    """

    
    # Storage for data needed for Training AND Inference
    batch_data = []
    
    print(f"Loading {len(json_files_list)} files for Global EM & Inference...")

    # --- PHASE 1: LOAD EVERYTHING ---
    
    init_params = None
    rr_val = 0
    
    for j_file in json_files_list:
        with open(j_file, 'r') as f:
            data = json.load(f)
            
        if init_params is None:
            prms = data["parameters_initial"]
            gen_time, mu, rr, l = prms['generation_time'], prms['mutation'], prms['rr'], prms['window_length']
            d = mu * l / gen_time
            init_params = [d*prms['t_archaic_c'], d*prms['t_split_c'], d*prms['t_introgression_c'], d*prms['t_introgression'], prms['admixture_proportion']]
            rr_val = rr

        obs_seq, _ = hmm.create_observations(data["prefix"]+'/'+data["data"], data["prefix"]+'/'+data["window_callability"]["Thousand_genomes"])
        

        
        
        cal_1kG = np.loadtxt(data['prefix'] + '/' + data["window_callability"]["Thousand_genomes"], usecols=-1)
        cal_nd_1kG = np.loadtxt(data['prefix'] + '/' + data["window_callability"]["Nd_1k_genomes"], usecols=-1)
        O1, O2, names = hmm.prepare_matrices_from_dict(obs_seq)
        M, N = O1.shape
        
        L1_mat = np.tile(cal_1kG, (M, 1))
        L2_mat = np.tile(cal_nd_1kG, (M, 1))
        
        batch_data.append({
            "em": (O1, O2, L1_mat, L2_mat),
            "viterbi": (cal_1kG, cal_nd_1kG, names),
            "json": data
        })
        
        gc.collect()

    print("Data loaded. Starting Global EM optimization...")

    # --- PHASE 2: GLOBAL TRAINING ---
    
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
            
        new_lmbd = total_nums / (total_dens + 1e-10)
        
        diff = iter_log_lik - prev_log_lik
        print(f"Iter {it+1}: LL={iter_log_lik:.2f} | Rates: N={new_lmbd[0]:.5f}, AF={new_lmbd[1]:.5f}, I={new_lmbd[2]:.5f}")
        
        if abs(diff) < tol and it > 0:
            print("Converged.")
            break
            
        prev_log_lik = iter_log_lik
        curr_lmbd = new_lmbd

    final_params = np.concatenate([curr_lmbd, [Ti, a_param]])
    print(f"Final Global Parameters: {final_params}")
    
    # --- PHASE 3: INFERENCE & SAVING ---
    
    print("Starting Inference on all files...")
    
    all_results = [] 
    
    for item in batch_data:
        O1, O2, _, _ = item["em"]
        cal_1kG, cal_nd_1kG, names = item["viterbi"]
        jsn_data = item["json"]
        
        result = hmm.run_hmm(O1, O2, cal_1kG, cal_nd_1kG, lmbd=final_params, rr=rr_val)
        
        dictionary = {k: v for k, v in zip(names, result)}
        out_dict = {}
        for name in names:
            out_dict[name] = hmm.get_tracts(dictionary[name])
            
        out_dict_new = hmm.clean_gaps(out_dict, jsn_data["gaps"], jsn_data["CHROM"])
        
        output_tsv = f"{jsn_data['prefix']}/{jsn_data['output']}.em.tsv"
        
        with open(output_tsv, "w", encoding="utf-8") as f:
            f.write("Sample\tCHROM\tStart\tEnd\tLength\n")
            for sample_name, tracks in out_dict_new.items():
                archaic_intervals = tracks.get('Archaic', [])
                for start, end in archaic_intervals:
                    f.write(f"{sample_name}\t{jsn_data["CHROM"]}\t{start}\t{end}\t{end - start +1}\n")
                    
                    if output_combined_file:
                        all_results.append({
                            "CHR": jsn_data["CHROM"],
                            "Sample": sample_name,
                            "Start": start,
                            "End": end,
                            "Length": end - start + 1
                        })
                    
    print(f"Done! Processed {len(batch_data)} files.")
    
    # --- SAVING MERGED FILE ---
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
