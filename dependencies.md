# DAISeg Module Dependencies

## Core Modules

daiseg.py
├── imports: hmm, em_alg, argparse
└── calls: hmm.run_daiseg(), em_alg.run_batch_em_pipeline()

hmm.py
├── imports: numpy, numba, scipy.stats, pandas, obs
├── called by: daiseg.py (run mode), em_alg.py (inference phase)
└── key functions: run_hmm(), viterbi_fast(), compute_emissions(), clean_gaps()

em_alg.py
├── imports: numpy, numba, hmm, pandas, gc, sys, os
├── called by: daiseg.py (run.with.EM mode)
└── key functions: run_batch_em_pipeline(), train_em_normalized(), e_step_normalized()
text


## Helper Scripts

```bash
extract.samples.sh
├── called by: daiseg.py (restrict_1kG mode)
├── requires: bcftools, jq
└── input: JSON config, outputs: filtered VCF

callability.sh
├── called by: daiseg.py (callability mode)
├── requires: bedtools, jq, sort
└── input: JSON config, outputs: coverage BED files

main.prep.py
├── called by: daiseg.py (main.prep mode)
├── imports: pysam, subprocess, time, preprocessing
├── requires: bcftools (external)
└── input: JSON config, outputs: final TSV for HMM
```


## External Dependencies

| Tool/Library | Used By | Purpose |
|--------------|---------|---------|
| numpy | hmm.py, em_alg.py | Array operations |
| numba | hmm.py, em_alg.py | JIT compilation, parallel loops |
| scipy | hmm.py | Poisson distribution |
| pandas | hmm.py, em_alg.py | DataFrame output |
| pysam | main.prep.py | Read ancestral FASTA |
| bcftools | extract.samples.sh, main.prep.py | VCF manipulation |
| bedtools | callability.sh | BED operations |
| jq | extract.samples.sh, callability.sh | JSON parsing |

## Call Graph
```bash
daiseg.py
  ├── run: hmm.py
  │     └── (obs.py, numpy, numba, scipy)
  ├── run.with.EM: em_alg.py
  │     └── (hmm.py, numpy, numba)
  ├── restrict_1kG: extract.samples.sh
  │     └── (bcftools)
  ├── callability: callability.sh
  │     └── (bedtools)
  └── main.prep: main.prep.py
        └── (preprocessing.py, pysam, bcftools)
```
