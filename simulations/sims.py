import msprime, tskit
import pandas as pd
import numpy as np
import os, json, gc, sys
import subprocess

# ==============================================================================
# 1. DEMOGRAPHY
# ==============================================================================

def history_archaic(prms, ne, t, p_admix, seed):
    """Simulate demography with introgression."""
    # Define pops
    demography = msprime.Demography()
    demography.add_population(name="AF", initial_size=ne['af'])
    demography.add_population(name="EU", initial_size=ne['eu'])
    demography.add_population(name="AMH", initial_size=ne['amh'])
    demography.add_population(name="ND", initial_size=ne['nd'])
    demography.add_population(name="ANCES", initial_size=ne['anc'])  
    demography.add_population(name="OOA", initial_size=ne['ooa'])
    
    # Events
    demography.add_population_parameters_change(time=0, initial_size=ne['eu'], population="EU", growth_rate=0.00202)
    demography.add_population_parameters_change(time=t['t_eu_growth'], initial_size=ne['eu_growth'], population="EU", growth_rate=0)
    
    # Admixture
    demography.add_admixture(time=t['t_nd_migration'], derived="EU", ancestral=["OOA", "ND"], proportions=[1-p_admix, p_admix])
    
    # Splits
    demography.add_population_split(time=t['t_ooa'], derived=["AF", "OOA"], ancestral="AMH")
    demography.add_population_split(time=t['t_amh'], derived=["AMH", "ND"], ancestral="ANCES")
    
    # Run sim
    ts = msprime.sim_ancestry(
        samples=[       
            msprime.SampleSet(prms['n_eu'], ploidy=prms['ploidy'], population='EU'), 
            msprime.SampleSet(prms['n_af'], ploidy=prms['ploidy'], population='AF'),
            msprime.SampleSet(prms['n_nd'], ploidy=prms['ploidy'], population='ND', time=t['t_nd_samples'])       
        ],    
        ploidy=prms['ploidy'],    
        sequence_length=prms['chrom_length'],
        recombination_rate=prms['recomb_rate'], 
        demography=demography,
        random_seed=seed, 
        record_migrations=True   
    )
    
    # Mutations
    ts = msprime.sim_mutations(ts, rate=prms['mut_rate'], random_seed=seed)    
    return ts


# ==============================================================================
# 2. GROUND TRUTH TRACTS
# ==============================================================================

def get_migrating_tracts_ind(ts, pop_name, ind_node, T_anc):
    """Get migration intervals for a specific node."""
    # Find pop ID
    pop_id = -1
    for p in ts.populations():
        if p.metadata.get('name') == pop_name: pop_id = p.id; break
    if pop_id == -1: raise ValueError(f"Pop '{pop_name}' not found")

    # Filter migrations
    tables = ts.tables
    mask = (tables.migrations.time == T_anc) & (tables.migrations.dest == pop_id)
    relevant = np.where(mask)[0]
    
    # Map nodes
    mig_lookup = {}
    for i in relevant:
        n = tables.migrations.node[i]
        interval = (tables.migrations.left[i], tables.migrations.right[i])
        if n not in mig_lookup: mig_lookup[n] = []
        mig_lookup[n].append(interval)

    # Traverse trees
    tracts = [] 
    for tree in ts.trees():
        anc = ind_node
        if ts.nodes_time[anc] > T_anc: continue
        
        parent = tree.parent(anc)
        while parent != tskit.NULL and ts.nodes_time[parent] <= T_anc:
            anc = parent
            parent = tree.parent(anc)
        
        if anc in mig_lookup:
            t_l, t_r = tree.interval
            for (m_l, m_r) in mig_lookup[anc]:
                s = max(t_l, m_l)
                e = min(t_r, m_r)
                if s < e:
                    if tracts and tracts[-1][1] == s: tracts[-1][1] = e
                    else: tracts.append([s, e])
    return tracts


def get_population_tracts_dataframe(ts, target_pop, source_pop, migration_time):
    """Return DataFrame of true tracts."""
    # Find IDs
    t_id = -1
    for p in ts.populations():
        if p.metadata.get('name') == target_pop: t_id = p.id; break
            
    nodes = ts.samples(population=t_id)
    if len(nodes) == 0: return pd.DataFrame(columns=["Sample", "Start", "End", "Length"])

    ind_ids = np.unique(ts.nodes_individual[nodes])
    ind_ids = ind_ids[ind_ids != -1]
    
    # Collect tracts
    data = []
    for ind in ind_ids:
        individual = ts.individual(ind)
        for i, node in enumerate(individual.nodes):
            name = f"{target_pop}_{ind}_{i + 1}"
            tracts = get_migrating_tracts_ind(ts, source_pop, node, migration_time)
            for s, e in tracts:
                data.append({"Sample": name, "Start": int(s), "End": int(e), "Length": int(e - s)})
                
    if not data: return pd.DataFrame(columns=["Sample", "Start", "End", "Length"])
    return pd.DataFrame(data)


# ==============================================================================
# 3. VARIANT TABLE
# ==============================================================================

def generate_haplotype_table(ts, chrom_name, assign_nucleotides=True):
    """Generate VCF-like TSV."""
    G = ts.genotype_matrix()
    pos = ts.tables.sites.position.astype(int)
    
    # Filter variants
    pop_map = {p.metadata['name']: p.id for p in ts.populations()}
    i_af = ts.samples(population=pop_map.get("AF"))
    i_nd = ts.samples(population=pop_map.get("ND"))
    i_eu = ts.samples(population=pop_map.get("EU"))

    eu_has_0 = np.any(G[:, i_eu] == 0, axis=1)
    eu_has_1 = np.any(G[:, i_eu] == 1, axis=1)
    af_has_0 = np.any(G[:, i_af] == 0, axis=1)
    af_has_1 = np.any(G[:, i_af] == 1, axis=1)
    nd_has_0 = np.any(G[:, i_nd] == 0, axis=1)
    nd_has_1 = np.any(G[:, i_nd] == 1, axis=1)
    
    mask = ((eu_has_0 & ~af_has_0) | (eu_has_1 & ~af_has_1)) | \
           ((eu_has_0 & ~nd_has_0) | (eu_has_1 & ~nd_has_1))
    
    G = G[mask]
    pos = pos[mask]
    n_vars = G.shape[0]
    
    if n_vars == 0:
        return pd.DataFrame(columns=["#CHROM", "POS", "REF", "ALT", "Ancestral", "Outgroup", "Neand"])

    # Assign bases
    bases = np.array(['A', 'C', 'G', 'T'])
    ref_idx = np.random.randint(0, 4, size=n_vars)
    refs = bases[ref_idx]
    alts = bases[(ref_idx + np.random.randint(1, 4, size=n_vars)) % 4]
    
    # Helper for pop columns
    def fmt_col(indices):
        if len(indices) == 0: return ["{}"] * n_vars
        sub = G[:, indices]
        h0 = np.any(sub == 0, axis=1)
        h1 = np.any(sub == 1, axis=1)
        res = []
        for i in range(n_vars):
            s = []
            if h0[i]: s.append(refs[i])
            if h1[i]: s.append(alts[i])
            s.sort()
            res.append("{" + ",".join(s) + "}")
        return res

    data = {
        "#CHROM": [chrom_name] * n_vars, "POS": pos, "REF": refs, "ALT": alts,
        "Ancestral": refs, "Outgroup": fmt_col(i_af), "Neand": fmt_col(i_nd)
    }
    
    # Add EU samples
    all_s = ts.samples()
    n_ind = ts.nodes_individual
    eu_map = {}
    for n in i_eu:
        ind = n_ind[n]
        if ind != -1: eu_map.setdefault(ind, []).append(n)
        
    for ind in sorted(eu_map.keys()):
        nodes = eu_map[ind]
        c1 = np.where(all_s == nodes[0])[0][0]
        data[f"EU_{ind}_1"] = np.where(G[:, c1] == 0, refs, alts)
        if len(nodes) > 1:
            c2 = np.where(all_s == nodes[1])[0][0]
            data[f"EU_{ind}_2"] = np.where(G[:, c2] == 0, refs, alts)
        else:
            data[f"EU_{ind}_2"] = ["."] * n_vars

    return pd.DataFrame(data)


# ==============================================================================
# 4. HELPERS
# ==============================================================================

def create_dummy_mask(d, name, length, w=1000):
    """Create perfect callability mask."""
    path = os.path.join(d, f"mask_{name}.bed")
    with open(path, 'w') as f:
        for i in range(int(length / w)):
            f.write(f"{name}\t{i*w}\t{(i+1)*w}\t1.0\n")
    return f"mask_{name}.bed"

# ==============================================================================
# 5. SIMULATION CONTROLLER
# ==============================================================================

def process_one_chromosome(seed, prms, ne, t, p_admix, out_dir):
    """Simulate data, save files, return ground truth."""
    np.random.seed(seed) 
    
    # Sim
    ts = history_archaic(prms, ne, t, p_admix, seed)    
    ts.dump(os.path.join(out_dir, f"sim_seed_{seed}.trees"))
    
    # Save TSV
    c_str = str(seed)
    df_v = generate_haplotype_table(ts, c_str)
    tsv_name = f"variants_seed_{seed}.tsv"
    df_v.to_csv(os.path.join(out_dir, tsv_name), sep='\t', index=False)
    del df_v; gc.collect() 
    
    # Aux files
    mask = create_dummy_mask(out_dir, c_str, prms['chrom_length'])
    gap_name = f"gaps_{seed}.txt"
    open(os.path.join(out_dir, gap_name), 'w').close()

    # IDs
    def get_ids(nm):
        pid = [p.id for p in ts.populations() if p.metadata['name']==nm]
        if not pid: return []
        return sorted([int(i) for i in set(ts.nodes_individual[ts.samples(population=pid[0])]) if i!=-1])

    # Config JSON
    out_base = f"inferred_seed_{seed}"
    cfg = {        
        "data": tsv_name,        
        "description": "Sim",
        "CHROM": c_str,
        "prefix": out_dir, 
        "output": out_base,
        "gaps": os.path.join(out_dir, gap_name),        
        "window_callability": {"Thousand_genomes": mask, "Nd_1k_genomes": mask},        
        "samples": {"outgroup": get_ids("AF"), "ingroup": get_ids("EU"), "neand": get_ids("ND")},        
        "parameters_initial": {
            "admixture_proportion": p_admix, "introgression_time": int(t['t_nd_migration'] * prms['gen_time']),
            "rr": prms['recomb_rate'], "mutation": prms['mut_rate'], "window_length": 1000,
            "generation_time": int(prms['gen_time']), "t_archaic_c": int(t['t_amh'] * prms['gen_time']),       
            "t_split_c": int(t['t_ooa'] * prms['gen_time']), "t_introgression_c": int(t['t_nd_migration'] * prms['gen_time']), 
            "t_introgression": int(t['t_nd_migration'] * prms['gen_time'])
        }
    }

    with open(os.path.join(out_dir, f"config_seed_{seed}.json"), 'w') as f:
        json.dump(cfg, f, indent=4)
    
    # Ground Truth
    df_t = get_population_tracts_dataframe(ts, "EU", "ND", t['t_nd_migration'])
    if not df_t.empty:
        df_t.insert(0, 'CHR', seed)        
        df_t[['Start','End','Length']] = df_t[['Start','End','Length']].astype(int)
    else:
        df_t = pd.DataFrame(columns=['CHR', 'Sample', 'Start', 'End', 'Length'])
    
    return df_t


# ==============================================================================
# 6. RUN DAISEG
# ==============================================================================

def run_daiseg_task(seed, out_dir):
    """Run daiseg.py via subprocess and return formatted results."""
    cfg = os.path.join(out_dir, f"config_seed_{seed}.json")
    res_file = os.path.join(out_dir, f"inferred_seed_{seed}.tsv")
    
    # FIX: Define path to daiseg.py in parent folder
    # sims.py is in ./simulations, daiseg.py is in ./
    sims_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.dirname(sims_dir)
    daiseg_path = os.path.join(root_dir, "daiseg.py")
    
    # Subprocess
    cmd = [sys.executable, daiseg_path, "run", "-json", cfg]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"Err seed {seed}: {e.stderr.decode()}")
        return None

    # Load Results
    df = None
    if os.path.exists(res_file):
        try:
            tmp = pd.read_csv(res_file, sep='\t')
            if not tmp.empty:
                if 'CHR' not in tmp.columns: tmp.insert(0, 'CHR', seed)
                tmp['Length'] = tmp['End'] - tmp['Start']
                df = tmp
        except: pass 
    return df


def calculate_accuracy(true_pth, inf_pth):
    """Calc Precision/Recall/F1."""
    # Load
    df_t = pd.read_csv(true_pth, sep='\t') if isinstance(true_pth, str) else true_pth
    try:
        df_i = pd.read_csv(inf_pth, sep='\t') if isinstance(inf_pth, str) else inf_pth
    except:
        df_i = pd.DataFrame(columns=df_t.columns)

    if df_i is None or df_i.empty:
        t_bp = df_t['Length'].sum() if not df_t.empty else 0
        return {"True_BP": int(t_bp), "Inferred_BP": 0, "Overlap_BP": 0, "Precision": 0, "Recall": 0, "F1": 0}

    # Stats
    t_bp = df_t['Length'].sum()
    i_bp = df_i['Length'].sum()
    if t_bp == 0 and i_bp == 0: return {"Precision": 0.0, "Recall": 0.0, "F1": 0.0}

    # Overlap
    overlap = 0
    i_grp = df_i.groupby(['CHR', 'Sample'])
    t_grp = df_t.groupby(['CHR', 'Sample'])
    i_dct = {n: g[['Start', 'End']].values for n, g in i_grp}

    for name, g_t in t_grp:
        if name in i_dct:
            invs = i_dct[name]
            for _, r in g_t.iterrows():
                ts, te = r['Start'], r['End']
                for iss, ie in invs:
                    s, e = max(ts, iss), min(te, ie)
                    if s < e: overlap += (e - s)

    # Metrics
    p = overlap / i_bp if i_bp > 0 else 0
    r = overlap / t_bp if t_bp > 0 else 0
    f1 = 2 * (p * r) / (p + r) if (p + r) > 0 else 0

    return {
        "True_BP": int(t_bp), "Inferred_BP": int(i_bp), "Overlap_BP": int(overlap),
        "Precision": round(p, 4), "Recall": round(r, 4), "F1": round(f1, 4)
    }


def calculate_class_metrics(true_pth, inf_pth, L):
    """Calc metrics per class (Archaic/Modern)."""
    # Load
    df_t = pd.read_csv(true_pth, sep='\t') if isinstance(true_pth, str) else true_pth
    try:
        df_i = pd.read_csv(inf_pth, sep='\t') if isinstance(inf_pth, str) else inf_pth
    except:
        df_i = pd.DataFrame(columns=df_t.columns)

    if df_i is None: df_i = pd.DataFrame(columns=df_t.columns)

    # Totals
    pairs = set(zip(df_t['CHR'], df_t['Sample'])).union(set(zip(df_i['CHR'], df_i['Sample'])))
    tot_bp = len(pairs) * L
    
    # Overlap (TP Archaic)
    tp_arch = 0
    i_grp = df_i.groupby(['CHR', 'Sample'])
    i_dct = {n: g[['Start', 'End']].values for n, g in i_grp}
    t_grp = df_t.groupby(['CHR', 'Sample'])

    for name, g_t in t_grp:
        if name in i_dct:
            invs = i_dct[name]
            for _, r in g_t.iterrows():
                ts, te = r['Start'], r['End']
                for iss, ie in invs:
                    s, e = max(ts, iss), min(te, ie)
                    if s < e: tp_arch += (e - s)

    # Confusion Matrix
    tot_true = df_t['Length'].sum()
    tot_pred = df_i['Length'].sum()
    
    fp_arch = tot_pred - tp_arch
    fn_arch = tot_true - tp_arch
    tn_mod = tot_bp - (tp_arch + fp_arch + fn_arch)
    
    def get_stats(tp, fp, fn):
        p = tp/(tp+fp) if (tp+fp)>0 else 0
        r = tp/(tp+fn) if (tp+fn)>0 else 0
        f = 2*p*r/(p+r) if (p+r)>0 else 0
        return round(p, 4), round(r, 4), round(f, 4)

    pa, ra, fa = get_stats(tp_arch, fp_arch, fn_arch)
    pm, rm, fm = get_stats(tn_mod, fn_arch, fp_arch)

    return {
        "Total_BP": tot_bp,
        "Archaic": {"Precision": pa, "Recall": ra, "F1": fa, "TP": int(tp_arch), "FP": int(fp_arch), "FN": int(fn_arch)},
        "Modern":  {"Precision": pm, "Recall": rm, "F1": fm, "TP": int(tn_mod)}
    }
