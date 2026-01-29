# file to make observation sequence

import csv
import sys


def parse_set_fast(val_str):
    """
    Parses TSV string like '{A,G}' or '{T}' into a Python set.
    Returns empty set for '.' or empty strings.
    """
    if not val_str or val_str == '{}' or val_str == '.':
        return set()
    # Remove braces and split by comma
    return set(val_str.strip('{}').split(','))


def process_data(tsv_path, bed_path):
    """
    Reads BED (windows) and TSV (genotypes).
    """
    print(" [obs.py] Loading BED windows...")
    
    windows_by_chrom = {}
    all_windows_flat = []

    # 1. Load BED file
    try:
        with open(bed_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                # Skip headers or comments
                if parts[0].startswith('#') or parts[0].lower() == 'chrom':
                    continue
                
                # Normalize chromosome name
                chrom = parts[0]
                if not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                
                start = int(parts[1])
                end = int(parts[2])

                window = {
                    's': start,
                    'e': end,
                    'stats': None
                }
                all_windows_flat.append(window)

                if chrom not in windows_by_chrom:
                    windows_by_chrom[chrom] = []
                windows_by_chrom[chrom].append(window)
    except FileNotFoundError:
        raise FileNotFoundError(f"BED file not found: {bed_path}")

    print(f"[obs.py] Loaded {len(all_windows_flat)} windows. Processing TSV...")

    # 2. Process TSV file
    try:
        with open(tsv_path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            
            # Read Header
            header = next(reader, None)
            if not header:
                raise ValueError("TSV file is empty")

            # Map Columns
            header_map = {name: i for i, name in enumerate(header)}

            try:
                # Identify chromosome column
                if '#CHROM' in header_map:
                    idx_c = header_map['#CHROM']
                elif 'chr' in header_map:
                    idx_c = header_map['chr']
                else:
                    idx_c = header_map['CHROM']  # fallback

                idx_p = header_map['POS']
                idx_anc = header_map['Ancestral']
                idx_out = header_map['Outgroup']
                idx_nean = header_map['Neand']
            except KeyError as e:
                raise ValueError(f"Missing mandatory column in TSV: {e}")

            # Identify Haplotype Columns
            exclude_names = {'#CHROM', 'chr', 'CHROM', 'POS', 'REF', 'ALT', 'Ancestral', 'Outgroup', 'Neand'}
            
            hap_indices = []
            for i, h in enumerate(header):
                if h not in exclude_names:
                    hap_indices.append((h, i))
            
            hap_names = [x[0] for x in hap_indices]

            # Initialize stats
            for w in all_windows_flat:
                w['stats'] = {name: [0, 0] for name in hap_names}

            # Linear Scan Variables
            current_chrom = None
            current_windows_list = []
            win_idx = 0
            max_win_idx = 0

            for row_num, row in enumerate(reader):
                if row_num % 100000 == 0 and row_num > 0:
                    print(f"   📊 Processed lines: {row_num}...", end='\r')

                try:
                    # Sync coordinates: BED is 0-based, TSV is 1-based
                    raw_chrom = row[idx_c]
                    
                    # Normalize to 'chr21' format
                    if not raw_chrom.startswith('chr'):
                        chrom = 'chr' + raw_chrom
                    else:
                        chrom = raw_chrom

                    pos = int(row[idx_p]) - 1
                except (ValueError, IndexError):
                    continue

                # Handle Chromosome Switch
                if chrom != current_chrom:
                    current_chrom = chrom
                    current_windows_list = windows_by_chrom.get(chrom, [])
                    win_idx = 0
                    max_win_idx = len(current_windows_list)

                if win_idx >= max_win_idx:
                    continue

                # Advance window index
                curr_win = current_windows_list[win_idx]
                while pos >= curr_win['e']:
                    win_idx += 1
                    if win_idx >= max_win_idx:
                        break
                    curr_win = current_windows_list[win_idx]
                
                if win_idx >= max_win_idx:
                    continue

                # Skip if SNP is before start of current window
                if pos < curr_win['s']:
                    continue

                # === INSIDE WINDOW: CALCULATE STATS ===
                anc = row[idx_anc]
                
                if not anc.isupper():
                    continue

                out_set = parse_set_fast(row[idx_out])
                nean_set = parse_set_fast(row[idx_nean])
                
                for hap_name, hap_idx in hap_indices:
                    val = row[hap_idx]
                    
                    if val == '.' or val == anc:
                        continue
                    if val.upper() == anc:
                        continue

                    is_diff_out = val not in out_set
                    is_diff_nean = val not in nean_set

                    if is_diff_out or is_diff_nean:
                        stats = curr_win['stats'][hap_name]
                        if is_diff_out:
                            stats[0] += 1
                        if is_diff_nean:
                            stats[1] += 1

    except FileNotFoundError:
        raise FileNotFoundError(f" Data file not found: {tsv_path}")
    
    print("\n [obs.py] Processing done. Aggregating results...")
    
    final_result = {name: [] for name in hap_names}
    for w in all_windows_flat:
        for name in hap_names:
            final_result[name].append(w['stats'][name])
            
    return final_result


def get_number_states(result_dict):
    """
    Finds maximum number of differences across all windows.
    """
    max_val = -1
    max_info = {}

    for hap, windows in result_dict.items():
        for i, stats in enumerate(windows):
            if stats[0] > max_val:
                max_val = stats[0]
                max_info = {
                    'hap': hap,
                    'win_idx': i,
                    'type': 'Outgroup',
                    'pair': stats
                }
            if stats[1] > max_val:
                max_val = stats[1]
                max_info = {
                    'hap': hap,
                    'win_idx': i,
                    'type': 'Neanderthals',
                    'pair': stats
                }
    return max_val, max_info

