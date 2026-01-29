import json
import multiprocessing
import sys
from pprint import pprint

import numpy as np
import pandas as pd
from numba import jit, prange
from scipy.stats import poisson

import obs  # module for creating observation sequence


def create_observations(tsv, bed):
    """Generates observation sequences from input data and determines state dimensions."""
    try:
        result = obs.process_data(tsv, bed)
        print(' Observation sequences for HMM created successfully!')
    except Exception as e:
        print(f"!!! Critical error: {e}")
        sys.exit(1)
    
    # Get number of states
    max_val, max_info = obs.get_number_states(result)
    print(f'Max differences in 1000bp window: {max_val + 1}')
    
    return result, max_val + 1


def prepare_matrices_from_dict(data_dict):
    """Converts haplotype dictionary into NumPy matrices (O1, O2) preserving insertion order."""
    
    # 1. Keep original insertion order
    hap_names = list(data_dict.keys())
    
    M = len(hap_names)  # Number of haplotypes
    if M == 0:
        raise ValueError("Dictionary contains no haplotypes.")
    
    N = len(data_dict[hap_names[0]])  # Number of windows
    
    # 2. Initialize matrices
    O1 = np.zeros((M, N), dtype=np.int32)
    O2 = np.zeros((M, N), dtype=np.int32)
    
    # 3. Fill matrices
    for i, hap in enumerate(hap_names):
        hap_data = np.array(data_dict[hap])  # shape (N, 2)
        O1[i, :] = hap_data[:, 0]  # Column 0 → O1
        O2[i, :] = hap_data[:, 1]  # Column 1 → O2
        
    return O1, O2, hap_names


# Ti: Introgression of Nd
# Taf: Time out of Africa
# Tn: Time of Split between Nd and Sapiens

# Transition probabilities
def initA(L, rr, Ti, a) -> np.array:
    """Calculates standard transition probability matrix based on recombination and admixture."""
    A = np.zeros((2, 2))
    
    A[0][1] = Ti * rr * L * a
    A[0][0] = 1 - A[0][1]
    
    A[1][0] = Ti * rr * L * (1 - a)
    A[1][1] = 1 - A[1][0]
    
    return A


# Log-transition probabilities
def get_log_A(L, rr, Ti, a):
    """Calculates transition probabilities in log-space to prevent numerical underflow."""
    A = np.zeros((2, 2))
    
    # Probability of recombination event
    prob = Ti * rr * L
    if prob > 0.5: 
        prob = 0.5  # Safety cap
    
    # Transitions: State 0 → 1
    A[0, 1] = prob * a
    A[0, 0] = 1.0 - A[0, 1]
    
    # Transitions: State 1 → 0
    A[1, 0] = prob * (1.0 - a)
    A[1, 1] = 1.0 - A[1, 0]
    
    return np.log(A + 1e-300)


def initB_arch_cover(lmbd, n_st, cover_1k, cover_nd):
    """Computes full Poisson emission probability matrix (reference implementation)."""
    
    # 1. Define lambdas (means)
    mean_n = lmbd[1] * cover_nd
    mean_n2 = lmbd[1] * cover_1k
    mean_af = lmbd[2] * cover_1k
    mean_i2 = lmbd[0] * cover_nd
    
    # 2. Helper function to generate probability vectors
    def get_prob_vec(mu):
        k = np.arange(1, n_st)  # Indices from 1 to n_st-1
        probs = poisson.pmf(k, mu)  # Vectorized Poisson calculation
        p0 = 1.0 - np.sum(probs)  # Residual mass for index 0
        return np.concatenate(([p0], probs))
    
    # 3. Generate vectors
    Paf = get_prob_vec(mean_af)
    Pn = get_prob_vec(mean_n)
    Pn2 = get_prob_vec(mean_n2)
    Pi2 = get_prob_vec(mean_i2)
    
    # 4. Construct emission matrix
    B = np.empty((2, n_st, n_st))
    B[0] = np.outer(Paf, Pn)
    B[1] = np.outer(Pn2, Pi2)
    
    return B


def compute_emissions_custom(O1, O2, L1, L2, lmbd):
    """Computes vectorized log-emission scores for all haplotypes using simplified Poisson model."""
    M, N = O1.shape
    n_states = 2
    
    # Output matrix (M x N x 2 states)
    log_emit = np.zeros((M, N, n_states))
    
    # Extract rate parameters
    rate_i = lmbd[2]
    rate_n = lmbd[0]
    rate_af = lmbd[1]
    rate_i = lmbd[2]
    
    # Epsilon to avoid log(0)
    eps = 1e-300
    ln_af = np.log(rate_af + eps)
    ln_n = np.log(rate_n + eps)
    ln_i = np.log(rate_i + eps)
    
    # Poisson Score = O * ln(rate) - rate * L (ln(O!) cancels)
    log_emit[:, :, 0] = (O1 * ln_af - rate_af * L1) + (O2 * ln_n - rate_n * L2)
    log_emit[:, :, 1] = (O1 * ln_n - rate_n * L1) + (O2 * ln_i - rate_i * L2)
    
    return log_emit


# VITERBI ALGORITHM
@jit(nopython=True, parallel=True)
def viterbi_fast(log_emit, log_trans, log_start):
    """Executes Viterbi algorithm in parallel using Numba to find most likely state path."""
    M, N, n_states = log_emit.shape
    paths = np.zeros((M, N), dtype=np.int32)
    
    # Parallel loop over all haplotypes
    for m in prange(M):
        # Local buffers for dynamic programming
        viterbi = np.zeros((N, n_states))
        backpointer = np.zeros((N, n_states), dtype=np.int32)
        
        # 1. Initialization
        for s in range(n_states):
            viterbi[0, s] = log_start[s] + log_emit[m, 0, s]
        
        # 2. Forward Pass
        for i in range(1, N):
            for s in range(n_states):
                # Find max(prev_prob + transition)
                max_val = -1e200
                best_prev = 0
                for p in range(n_states):
                    val = viterbi[i - 1, p] + log_trans[p, s]
                    if val > max_val:
                        max_val = val
                        best_prev = p
                
                # Add emission score for current window
                viterbi[i, s] = max_val + log_emit[m, i, s]
                backpointer[i, s] = best_prev
        
        # 3. Backtrace
        # Decide final state
        if viterbi[N - 1, 0] > viterbi[N - 1, 1]:
            paths[m, N - 1] = 0
        else:
            paths[m, N - 1] = 1
        
        # Reconstruct path backwards
        for i in range(N - 2, -1, -1):
            paths[m, i] = backpointer[i + 1, paths[m, i + 1]]
    
    return paths


def run_hmm(O1, O2, L1, L2, lmbd, rr):
    """Orchestrates HMM pipeline: emissions, transitions, and Viterbi."""
    
    print("Calculating emission scores...")
    log_emissions = compute_emissions_custom(O1, O2, L1, L2, lmbd)
    
    # Transitions
    log_A = get_log_A(1000, rr, lmbd[3], lmbd[4])
    
    # Initial probabilities
    log_start = np.log([1.0 - lmbd[4], lmbd[4]])
    
    print("Running Viterbi...")
    paths = viterbi_fast(log_emissions, log_A, log_start)
    
    return paths


def get_tracts(vector, step=1000):
    """Converts binary state vector into dictionary of genomic intervals (start, end)."""
    result = {0: [], 1: []}
    
    current_state = vector[0]
    start_index = 0
    
    for i, state in enumerate(vector):
        if state != current_state:
            result[current_state].append((start_index * step, i * step - 1))
            current_state = state
            start_index = i
    
    result[current_state].append((start_index * step, len(vector) * step))
    
    return {"Modern": result[0], "Archaic": result[1]}


def clean_gaps(dct, gap_file, target_chrom):
    """Filters genomic regions defined in gap file from inferred tracts."""
    
    print('Processing gaps...')
    raw_gaps = []
    
    try:
        with open(gap_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4 and parts[1] == 'chr' + target_chrom:
                    raw_gaps.append((int(parts[2]), int(parts[3]) - 1))
    except FileNotFoundError:
        print(" !!! Gap file not found.")
        return dct
    

    
    # Merge overlapping gaps
    merged_gaps = []
    if raw_gaps:
        raw_gaps.sort()
        merged_gaps = [raw_gaps[0]]
        for curr in raw_gaps[1:]:
            prev = merged_gaps[-1]
            if curr[0] <= prev[1] + 1:
                merged_gaps[-1] = (prev[0], max(prev[1], curr[1]))
            else:
                merged_gaps.append(curr)
    
    # Helper to subtract gaps from interval
    def subtract(interval, gaps):
        start, end = interval
        res = []
        curr = start
        for g_s, g_e in gaps:
            if g_e < curr:
                continue
            if g_s > end:
                break
            if curr < g_s:
                res.append((curr, g_s - 1))
            curr = max(curr, g_e + 1)
        if curr <= end:
            res.append((curr, end))
        return res
    
    # Process dictionary
    new_dct = {}
    for sample, categories in dct.items():
        new_dct[sample] = {}
        for cat, intervals in categories.items():
            cleaned_list = []
            for interval in intervals:
                if not merged_gaps:
                    cleaned_list.append(interval)
                else:
                    cleaned_list.extend(subtract(interval, merged_gaps))
            new_dct[sample][cat] = cleaned_list
    
    return new_dct


def run_daiseg(json_file):
    """Main pipeline: runs HMM and saves ONLY TSV output."""
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Create observation sequences
    tsv_path = data["prefix"] + '/' + data["data"]
    bed_path = data["prefix"] + '/' + data["window_callability"]["Thousand_genomes"]
    
    obs_seq, nst = create_observations(tsv_path, bed_path)
    
    if len(obs_seq) > 0:
        first_key = list(obs_seq.keys())[0]
#        print(f"Sample length: {len(obs_seq[first_key])}")

    
    # Extract parameters
    prms = data["parameters_initial"]
    gen_time = prms['generation_time']
    mu = prms['mutation']
    rr = prms['rr']
    l = prms['window_length']
    
    d = mu * l / gen_time
    lambda_0 = [
        d * prms['t_archaic_c'],
        d * prms['t_split_c'],
        d * prms['t_introgression_c'],
        d * prms['t_introgression'],
        prms['admixture_proportion']
    ]
    
    # Load callability files
    cal_1kG = np.loadtxt(data['prefix'] + '/' + data["window_callability"]["Thousand_genomes"], usecols=-1)
    cal_nd_1kG = np.loadtxt(data['prefix'] + '/' + data["window_callability"]["Nd_1k_genomes"], usecols=-1)
    
    # Prepare matrices
    O1, O2, names = prepare_matrices_from_dict(obs_seq)
    
    # Run HMM
    result = run_hmm(O1, O2, cal_1kG, cal_nd_1kG, lmbd=lambda_0, rr=rr)
    dictionary = {k: v for k, v in zip(names, result)}
    
    # Extract tracts
    out_dict = {}
    for name in names:
        out_dict[name] = get_tracts(dictionary[name])
    


    # Remove gaps
    out_dict_new = clean_gaps(out_dict, data["gaps"], data["CHROM"])
    
    # Save TSV results
    output_tsv = f"{data['prefix']}/{data['output']}.tsv"
    print(f" Saving TSV results to: {output_tsv}")
    
    rows = []
    with open(output_tsv, "w", encoding="utf-8") as f:
        f.write("Sample\tCHROM\tStart\tEnd\tLength\n")
        for sample_name, tracks in out_dict_new.items():
            archaic_intervals = tracks.get('Archaic', [])
            for start, end in archaic_intervals:
                f.write(f"{sample_name}\t{data['CHROM']}\t{start}\t{end}\t{end-start+1}\n")
                rows.append({
                    "Sample": sample_name,
                    "CHR": data['CHROM'],
                    "Start": start,
                    "End": end
                })
    
    df_result = pd.DataFrame(rows)
    return df_result, out_dict_new


# Store original logic for reference
_original_logic = run_daiseg


def _worker_proxy(filepath):
    """Helper function for pickle compatibility (must be defined globally)."""
    return _original_logic(filepath)


def run_daiseg(json_input):
    """
    Wrapper: handles single file (string) or list of files (parallel).
    """

    # Single file → run normally БЕЗ multiprocessing
    if not isinstance(json_input, list):
        print(f" Processing single file (no parallelization needed)...")
        return _original_logic(json_input)

    # Single file in list → тоже без multiprocessing
    if len(json_input) == 1:
        print(f" Processing single file (sequential)...")
        return [_original_logic(json_input[0])]

    # Уже в подпроцессе → выполнять последовательно
    if multiprocessing.current_process().daemon:
        return [_original_logic(f) for f in json_input]

    # Множество файлов → параллельно
    MAX_WORKERS = 64
    cpu_count = multiprocessing.cpu_count()
    pool_size = min(cpu_count - 1, MAX_WORKERS, len(json_input))
    pool_size = max(pool_size, 1)
    
    print(f" Parallelizing {len(json_input)} files on {pool_size} cores...")

    with multiprocessing.Pool(processes=pool_size) as pool:
        return pool.map(_worker_proxy, json_input)
