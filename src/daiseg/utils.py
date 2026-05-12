import json
import os
import sys

def expand_path(path):
    
    if not path:
        return None
    return os.path.expanduser(path)

def load_config(json_path):
    """Load .json configuration."""
    if not os.path.exists(json_path):
        sys.stderr.write(f"Error: File '{json_path}' not found.\n")
        sys.exit(1)
    with open(json_path, 'r') as f:
        return json.load(f)

def map_columns(header_list, samples_cfg):
    """
    Map VCF header columns to sample groups.
    CRITICAL: Enforces strict order from JSON for Ingroup/Outgroup.
    """
    try:
        idx_chrom = header_list.index("CHROM")
        idx_pos = header_list.index("POS")
        idx_ref = header_list.index("REF")
        idx_alt = header_list.index("ALT")
    except ValueError as e:
        raise ValueError(f"Missing standard columns in header: {e}")

    # Create a lookup map: {Name: Index}
    header_map = {name.split(':')[0]: i for i, name in enumerate(header_list)}

    s1_names_list = samples_cfg.get("outgroup", [])
    s3_names_list = samples_cfg.get("ingroup", [])
    
    s1_set = set(s1_names_list)
    s3_set = set(s3_names_list)

    col_map_yri = [] 
    col_map_ibs = [] 
    col_map_nd = []  

    # 1. Outgroup (S1) - Follow JSON order
    for name in s1_names_list:
        if name in header_map:
            col_map_yri.append(header_map[name])
        else:
            sys.stderr.write(f"[WARN] Sample {name} (Outgroup) not found in VCF.\n")

    # 2. Ingroup (S3) - Follow JSON order
    for name in s3_names_list:
        if name in header_map:
            col_map_ibs.append(header_map[name])
        else:
            sys.stderr.write(f"[WARN] Sample {name} (Ingroup) not found in VCF.\n")

    # 3. Neanderthals (S2) - Catch remaining samples
    start_idx = idx_alt + 1
    for i in range(start_idx, len(header_list)):
        raw_name = header_list[i]
        name = raw_name.split(':')[0]
        
        if name not in s1_set and name not in s3_set:
            col_map_nd.append(i)
            
    return idx_chrom, idx_pos, idx_ref, idx_alt, col_map_yri, col_map_ibs, col_map_nd
