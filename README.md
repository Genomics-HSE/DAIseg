# Workflow

Below is the complete execution pipeline for chromosome 22.

### 1. Data Restriction
Filters the 1000 Genomes VCFs based on the configuration.
```bash
python daiseg.py restrict_1kG -json all.chr22.json -threads 8
```

### 2. Callability Mask
Calculates the genomic windows accessible for analysis (filters masks).
```bash
python daiseg.py callability -json all.chr22.json -threads 8
```

### 3. Main Preprocessing
Merges VCFs, filters SNPs, and creates the observation matrix (TSV).
```bash
python daiseg.py main.prep -json all.chr22.json -threads 8
```

### 4. Running HMM
Runs the Hidden Markov Model to infer introgression tracts.
```bash
python daiseg.py run -json all.chr22.json
```

### 5. Using EM for estimation
Runs the Hidden Markov Model to infer introgression tracts.
```bash
python daiseg.py run.with.EM -json all.chr22.json
```

---

###  Full Pipeline Script

```bash
#!/bin/bash

# Settings
CONF="all.chr22.json"
THR=8

echo "--- [1/4] Restricting 1000 Genomes ---"
python daiseg.py restrict_1kG -json $CONF -threads $THR

echo "--- [2/4] Calculating Callability ---"
python daiseg.py callability -json $CONF -threads $THR

echo "--- [3/4] Preprocessing Observations (VCF -> TSV) ---"
python daiseg.py main.prep -json $CONF -threads $THR

echo "--- [4/4] Running HMM  ---"
python daiseg.py run -json $CONF

echo " Pipeline finished successfully."
```
