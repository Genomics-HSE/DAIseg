# DAISeg Module Architecture

## Call Graph
A visual overview of how the main entry point interacts with internal modules and external tools.

```text
daiseg.py (Main Entry Point)
 ├── run ⮕ hmm.py (Viterbi & Emissions)
 │    └── Dependencies: obs.py, numpy, numba, scipy
 ├── run.with.EM ⮕ em_alg.py (Parameter Training)
 │    └── Dependencies: hmm.py, numpy, numba
 ├── main.prep ⮕ main.prep.py (Data Wrangling)
 │    └── Dependencies: preprocessing.py, pysam, bcftools
 ├── restrict_1kG ⮕ extract.samples.sh (VCF Filtering)
 │    └── Dependencies: bcftools, jq
 └── callability ⮕ callability.sh (BED Processing)
      └── Dependencies: bedtools, jq
```

## Component Logic

Functional roles of the core modules and scripts:

- **daiseg.py**: Primary CLI interface for toggling between processing modes.

- **hmm.py**: Core computational engine. Implements Viterbi algorithm and emission calculations (Numba-optimized).

- **em_alg.py**: Training module. Uses the EM algorithm for automated parameter estimation.

- **main.prep.py**: Data preparation. Converts genomic formats into TSV files compatible with the HMM.

- **Shell Helpers**: Low-level processing of heavy genomic files (VCF/BED) using specialized bioinformatics tools.


## Environment & Dependencies

| Category | Tools / Libraries | Purpose |
|----------|-------------------|---------|
| Python Core | numpy, numba, scipy, pandas | Numerical computing and JIT acceleration |
| Genomics | pysam, bcftools, bedtools | Sequence alignment and VCF/BED manipulation |
| Utilities | jq, sort | JSON parsing and data ordering |
