# DAIseg 
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
To run DAIseg in its simplest way you need callability and prepared .tsv file.

### 1. Create file with general callability information of modern genomes

Processes a chromosome to create a .BED file with callability coverage statistics calculated in 1000 bp windows. Output BED Format  is 
Tab-separated values with columns 

```bash
chr  start_i  end_i  num_variants  pos_in_mask  window_length  coverage
```


**Processing Logic.**

    - Window Creation – Divides chromosome into non-overlapping window_size bp windows.
    
    - Variant Counting – For each window, counts variants whose coordinates overlap with window.
    
    - Mask Coverage – Calculates number of base positions overlapping the callability mask.
    
    - Coverage Calculation – Computes coverage = pos_in_mask / window_length.


**Example.**
```bash
1    0       999    12    980    1000    0.980
1    1000    1999    8     995    1000    0.995
1    2000    2999    15    876    1000    0.876
```

### 1.1 Create file with general information in significant genome positions of archaic genomes

Processes a chromosome to create a .BED file with neanderthal callability calculated in 1000 bp windows. Output BED Format  is 
Tab-separated values with columns 

```bash
chr  start_i  end_i  num_variants  pos_in_mask  window_length  coverage
```


**Processing Logic.**
    - Window Creation – Divides chromosome into non-overlapping window_size bp windows.
    - Variant Counting – For each window, counts variants whose coordinates overlap.
    - Mask Coverage – Calculates number of base positions of **JOINED** neanderthals overlapping the callability mask.
    - Coverage Calculation – Computes coverage = pos_in_mask / window_length.

### 2. Create file  like
```bash
CHROM    POS	REF	ALT	Ancestral	Outgroup	Neand	Sample1_hap1	Sample1_hap2    ...   SampleN_hap1	SampleN_hap2 
```
where each row corresponds to a single biallelic position where at least one difference exists in the target samples {Sample1.. SampleN} relative to Africans(Outgroup) or Neanderthals(Neand). REF, ALT and Ancestral are reference, alternative and ancestral alleles respectively. 

### 3. Create config .json file

```bash
{
  "description": "DAIseg.simple configuration to run ",
  "CHROM": "chr",  
  "output": "out.assembly.chr", 
  "prefix": "/path/to/output/directory/prefix",  
  "files": {
    "neand_files": {  # Neanderthal data files
      "Vindija33.19": {
        "bed": "/path/to/neand/vindija33.19/.bed",  
        "vcf": "/path/to/neand/vindija33.19/.vcf.gz" 
      },
      "Altai": {
        "bed": "/path/to/neand/altai/.bed", 
        "vcf": "/path/to/neand/altai/.vcf.gz" 
      }
    },
    "1000GP_files": {  # 1000 Genomes Project files
      "bed": "/path/to/1000gp/strict/mask/.bed",  
      "vcf": "preprocessed.vcf.gz",  
      "vcf_initial": "/path/to/1000gp/.vcf.gz"  
    },
    "ancestral": {
      "fasta": "/path/to/ancestral/reference/.fa" 
    },
    "reference": {
      "fasta": "/path/to/reference/genome/.fa"  
    },
    "chr_lengths": "/path/to/chromosome/lengths/assembly.sizes"  
  },
  "samples": {  
    "outgroup": [  
      "OUT1", "OUT2", "OUT3", "OUT4", 
      ...
      ],
    "ingroup": [  
      "IN1", "IN2", "IN3", "IN4",
      ... 
    ],
    "neand": [  
      "Vindija33.19",
      "AltaiNeandertal"
    ]
  },
  "parameters_initial": {  
    "admixture_proportion": 0.02,  
    "introgression_time": 55000,  
    "rr": 1e-08,  
    "mutation": 1.25e-08,  
    "window_length": 1000, 
    "generation_time": 29,  
    "t_archaic_c": 550000,  
    "t_split_c": 70000,  
    "t_introgression_c": 55000,  
    "t_introgression": 55000  
  },
  "window_callability": { 
    "Thousand_genomes": "coverage_file.bed", 
    "Nd_1k_genomes": "coverage_file.bed"  
  },
  "data": "preprocessed_data.tsv",  
  "gaps": "/path/to/genome/gaps/gap.txt" 
}
```


### 4. Running HMM
Runs the Hidden Markov Model to infer introgression tracts.
```bash
python daiseg.py run -json examle.json
```

### 5. Using EM for estimation
Runs the Hidden Markov Model to infer introgression tracts without transition estimates:
```bash
python daiseg.py run.with.EM -json example.json
```

Runs the Hidden Markov Model to infer introgression tracts with transitions estimates:
```bash
python daiseg.py run.EM.v2 -json example.json
```



# 1000GP Workflow 

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
