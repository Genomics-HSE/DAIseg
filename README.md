# DAIseg (grch38 is under construction)
**A Hidden Markov Model (HMM) for detecting archaic introgression in modern genomes.**

## 📌 Overview
DAIseg is a highly accurate method to identify genomic segments in modern humans inherited from archaic admixture. 


The introgression scenario describes a historical admixture event where an ancestral modern human population interbred with an archaic hominin group, such as Neanderthals, after their initial divergence. This event, which occurred at a specific time in the past (e.g., tens of thousands of years ago), introduced a small proportion of archaic DNA (typically 1-2%) into the gene pool of the non-African population. The resulting genetic signature—archaic genomic segments embedded within modern human lineages—is what DAIseg is designed to detect.

## 🔧 Core Method
- **Model:** Hidden Markov Model (HMM)
* **Dual-Reference Strategy:** Leverages both archaic and unadmixed modern reference data to minimize false positives.
- **Output:** Genomic tracts of archaic ancestry with high precision/recall.

## 📖 Reference
Planche, L., Ilina, A.V., & Shchur, V.L. (2024). Highly Accurate Method for Detecting Archaic Segments in the Modern Genomes. *Lobachevskii J Math*, 45, 2910–2917.  [https://doi.org/10.1134/S1995080224602959](https://doi.org/10.1134/S1995080224602959)


# General workflow 

### 1. Create file with general information in significant genome positions
Processing Logic

    Window Creation: Divides chromosome into non-overlapping window_size bp windows

    Variant Counting: For each window, counts variants whose coordinates overlap the window

    Mask Coverage: Calculates number of base positions overlapping the callability mask

    Coverage Calculation: Computes coverage = pos_in_mask / window_length

Processes a chromosome to create a BED file with callability coverage statistics calculated in 1000 bp windows.
**Output BED Format**
Tab-separated values with columns:
```bash
chr  start_i  end_i  num_variants  pos_in_mask  window_length  coverage
```


```bash
chr  start_i  end_i  num_variants  pos_in_mask  window_length  coverage
1    0       999    12    980    1000    0.980
1    1000    1999    8     995    1000    0.995
1    2000    2999    15    876    1000    0.876
```

### 2.




### 3. Running HMM
Runs the Hidden Markov Model to infer introgression tracts.
```bash
python daiseg.py run -json examle.json
```

### 4. Using EM for estimation
Runs the Hidden Markov Model to infer introgression tracts.
```bash
python daiseg.py run.with.EM -json example.json
```



# Workflow for samples (ingroup and outgroup) in 1000GP

Below is the complete execution pipeline for chromosome 22.

### 1. Data Restriction
Filters the 1000 Genomes VCFs based on the configuration.
```bash
python daiseg.py restrict_1kG -json example.json -threads 8
```

### 2. Callability Mask
Calculates the genomic windows accessible for analysis (filters masks).
```bash
python daiseg.py callability -json example.json -threads 8
```

### 3. Main Preprocessing
Merges VCFs, filters SNPs, and creates the observation matrix (TSV).
```bash
python daiseg.py main.prep -json example.json -threads 8
```

### 4. Running HMM
Runs the Hidden Markov Model to infer introgression tracts.
```bash
python daiseg.py run -json examle.json
```

### 5. Using EM for estimation
Runs the Hidden Markov Model to infer introgression tracts.
```bash
python daiseg.py run.with.EM -json example.json
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
