# hg19/GRCh37 Files Setup Guide

## Reference fasta

```bash
# Download GRCh37 reference genome (hs37d5)
wget ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
```
Description: Downloads the hs37d5 reference genome, which is the GRCh37 assembly with added decoy sequences used by the 1000 Genomes Project to improve mapping accuracy.

## Ancestral fasta
```bash

# Download ancestral genome for GRCh37 (Ensembl release 71)
wget ftp://ftp.ensembl.org/pub/release-71/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2
```
Description: Downloads the inferred ancestral allele sequences for human genome assembly GRCh37 from Ensembl release 71, used for evolutionary comparisons and identifying derived alleles.

## Chromosome Sizes
```bash

# Download chromosome sizes for hg19
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

# Remove 'chr' prefix from chromosome names
sed 's/^chr//' hg19.chrom.sizes > hg19.chrom.len
```

Description: Downloads chromosome lengths for hg19 assembly from UCSC and creates a version without 'chr' prefix in chromosome names. This format (numeric-only chromosome identifiers) is often required by bioinformatics tools.




## Genome Gaps
```bash

# Download hg19 gap annotations
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz

# Uncompress the file
zcat gap.txt.gz > gap.txt
```
Description: Downloads gap annotations for hg19 assembly from UCSC database. This file contains information about centromeres, telomeres, heterochromatin regions, and other unsequenced or problematic regions of the genome



## Genomes

### 1000 Genomes Project

```bash
# Download phase 3 VCF for chromosome 22
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

# Download corresponding tabix index
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
```

Description: Downloads the 1000 Genomes Project Phase 3 data for chromosome 22, containing phased genotype calls for all 2,504 samples across 85 populations, generated using the Shapeit2 and MVNcall pipelines (2013 release).

### Vindija 
```bash
# Download the compressed Altai Neanderthal VCF file 
 wget  http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr22_mq25_mapab100.vcf.gz

# Download the Tabix index file 
 wget http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr22_mq25_mapab100.vcf.gz.tbi
```
Description: Compressed VCF file containing variant calls for the Vindija 33.19 Neanderthal genome (chromosome 22) with quality filters: minimum mapping quality 25 (`mq25`) and 100% mappable regions (`mapab100`).


### Altai Neanderthal
```bash
# Download the compressed Altai Neanderthal VCF file 
wget http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr22_mq25_mapab100.vcf.gz

# Download the Tabix index file 
wget http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr22_mq25_mapab100.vcf.gz.tbi
```
Description: Compressed VCF file containing variant calls for the Altai Neanderthal genome (chromosome 22) with quality filters: minimum mapping quality 25 (`mq25`) and 100% mappable regions (`mapab100`).


## Genomics Masks for Callability

### 1. 1000 Genomes Project Strict Mask
```bash
# Download the strict mask from 1000 Genomes
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.strict_mask.autosomes.bed

# Split by chromosome
for chr in {1..22}; do
    grep -w "chr${chr}" 20140520.strict_mask.autosomes.bed > "chr${chr}.bed"
done
```
Description: Downloads the strict accessibility mask from 1000 Genomes Project (20140520 version) and splits it into individual chromosome files.

### 2. Vindija Neanderthal Mask
```bash

# Download Vindija 33.19 mask for chromosome 22 (example)
wget http://ftp.eva.mpg.de/neandertal/Vindija/FilterBed/Vindija33.19/chr22_mask.bed.gz
```

Description: Downloads the accessibility mask for Vindija 33.19 Neanderthal sample. Repeat for all chromosomes (1-22, X).


### 3. Altai Neanderthal Mask
```bash

# Download Altai Neanderthal mask for chromosome 22 (example)
wget http://ftp.eva.mpg.de/neandertal/Vindija/FilterBed/Altai/chr22_mask.bed.gz
```
Description: Downloads the accessibility mask for Altai Neanderthal sample. Repeat for all chromosomes (1-22, X).










# GRCh38 Reference Files Setup Guide
## Reference fasta
```bash
# Download complete GRCh38 reference 
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```
Description: Downloads the complete GRCh38 reference genome assembly used by 1000 Genomes Project
## Ancestral fasta
```bash
# Download ancestral genome for GRCh38 (latest Ensembl release)
wget https://ftp.ensembl.org/pub/release-114/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz
```
Description: Downloads inferred ancestral allele sequences for human genome assembly GRCh38 from Ensembl release 114, used for evolutionary analyses and identifying derived variants.
## Chromosome Sizes
```bash
# Download chromosome sizes for GRCh38
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
```
Description: Downloads chromosome lengths for GRCh38 assembly from UCSC, used for genome browser tracks and coordinate calculations.
## Genome Gaps
```bash
# Download gap annotations for GRCh38/hg38 assembly
curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
```
Description: Downloads gap annotations for the GRCh38/hg38 assembly from UCSC Genome Browser database. This file contains coordinates of genomic gaps including centromeres, telomeres, heterochromatin regions, and other unsequenced/problematic areas.

## Genomes
### 1000 Genome project
```bash
# Download phased VCF nad index file for each chromosome 
    wget  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz 
    wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz.tbi

```
Description: Downloads high-coverage (30x) phased genotype data from the 1000 Genomes Project (2020 release) for 2,504 samples, organized by chromosome directories. Data is filtered and phased using Shapeit2 and DuoHMM pipelines.
### Vindija and Altai Neanderthal
Download [VCF_Vindija](https://github.com/ekrez/daiseg.simple/blob/main/real.data.md#vindija) and [VCF_Altai](https://github.com/ekrez/daiseg.simple/blob/main/real.data.md#altai-neanderthal) in GRCh37 and make CrossMap to GRCh38 with chain file [GRCh37_to_GRCh38.chain.gz](https://42basepairs.com/browse/web/ensembl/assembly_mapping/homo_sapiens?file=GRCh37_to_GRCh38.chain.gz&preview=).
```bash
# Lifting Over
CrossMap vcf GRCh37_to_GRCh38.chain.gz VCF GRCh38 OUT --no-comp-allele --chromid l

# sorting
bcftools sort OUT  -Oz -o OUT.gz

# Indexing
tabix -p vcf OUT.gz 
```
 in GRCh37 and make CrossMap to GRCh38 with chain file [GRCh37_to_GRCh38.chain.gz](https://42basepairs.com/browse/web/ensembl/assembly_mapping/homo_sapiens?file=GRCh37_to_GRCh38.chain.gz&preview=).
```bash
# Lifting Over
CrossMap vcf GRCh37_to_GRCh38.chain.gz VCF GRCh38 OUT --no-comp-allele --chromid l

# sorting
bcftools sort OUT  -Oz -o OUT.gz

# Indexing
tabix -p vcf OUT.gz 
```



## Genomics Masks for Callability
### 1. 1000 Genomes Project Strict Mask
```bash
# Download strict accessibility mask for GRCh38
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed
```
Description: Downloads the strict accessibility mask for GRCh38 (20160622 version), defining reliably callable genomic regions across all chromosomes for 1000 Genomes Project data.
### 2. Vindija and Altai Neanderthal Mask
Download [BED_Vindija](https://github.com/ekrez/daiseg.simple/blob/main/real.data.md#2-vindija-neanderthal-mask) and [BED_Altai](https://github.com/ekrez/daiseg.simple/blob/main/real.data.md#3-altai-neanderthal-mask) in GRCh37 and make CrossMap to GRCh38 with chain file [GRCh37_to_GRCh38.chain.gz](https://42basepairs.com/browse/web/ensembl/assembly_mapping/homo_sapiens?file=GRCh37_to_GRCh38.chain.gz&preview=).
```bash
# Lifting over
CrossMap bed GRCh37_to_GRCh38.chain.gz BED OUT
    
# sorting with bedtools
bedtools sort -i OUT > OUT.sorted
    
# add chr prefix
sed 's/^/chr/' OUT.sorted > FINAL_BED
```



# .json input description

Main configuration file of DAIseg.simple:

```bash
{
  "description": "DAIseg.simple configuration to run ",
  "CHROM": "chr",  #chr name  - use as in assembly
  "output": "out.assembly.chr", 
  "prefix": "/path/to/output/directory/prefix",  # Path prefix for output files
  "files": {
    "neand_files": {  # Neanderthal data files
      "Vindija33.19": {
        "bed": "/path/to/neand/vindija33.19/.bed",  # BED file of accessible regions for Vindija 33.19 
        "vcf": "/path/to/neand/vindija33.19/.vcf.gz"  # VCF file with variants 
      },
      "Altai": {
        "bed": "/path/to/neand/altai/.bed",  # BED file of accessible regions for Altai Neanderthal 
        "vcf": "/path/to/neand/altai/.vcf.gz"  # VCF file with variants 
      }
    },
    "1000GP_files": {  # 1000 Genomes Project files
      "bed": "/path/to/1000gp/strict/mask/.bed",  # Strict accessibility mask 
      "vcf": "preprocessed.vcf.gz",  # Filtered VCF file for analysis
      "vcf_initial": "/path/to/1000gp/.vcf.gz"  # Original phased VCF file
    },
    "ancestral": {
      "fasta": "/path/to/ancestral/reference/.fa"  # FASTA file with ancestral sequence
    },
    "reference": {
      "fasta": "/path/to/reference/genome/.fa"  # Reference genome in FASTA format
    },
    "chr_lengths": "/path/to/chromosome/lengths/grch38.sizes"  # File with chromosome lengths (chr<tab>length format)
  },
  "samples": {  # Sample lists for analysis
    "outgroup": [  # African samples (outgroup)
      "OUT1", "OUT2", "OUT3", "OUT4", 
      ...
      ],
    "ingroup": [  # Non-African samples (ingroup)
      "IN1", "IN2", "IN3", "IN4",
      ... 
    ],
    "neand": [  # Neanderthal samples
      "Vindija33.19",
      "AltaiNeandertal"
    ]
  },
  "parameters_initial": {  # Initial parameters for the model
    "admixture_proportion": 0.02,  # Assumed admixture proportion (2%)
    "introgression_time": 55000,  # Time of introgression
    "rr": 1e-08,  # Recombination rate per base pair per generation
    "mutation": 1.25e-08,  # Mutation rate per base pair per generation
    "window_length": 1000,  # Window length for analysis (in bp)
    "generation_time": 29,  # Generation time (in years)
    "t_archaic_c": 550000,  # Divergence time with archaic lineage 
    "t_split_c": 70000,  # Coalescence  of non-African and African lineages 
    "t_introgression_c": 55000,  # Coalescent time 
    "t_introgression": 55000  # Introgression time 
  },
  "window_callability": {  # Files with window coverage/callability
    "Thousand_genomes": "coverage_file.bed",  # Coverage based on 1000 Genomes data
    "Nd_1k_genomes": "coverage_file.bed"  # Combined coverage (Neanderthal + 1000G)
  },
  "data": "preprocessed_data.tsv",  # Preprocessed data for analysis
  "gaps": "/path/to/genome/gaps/gap.txt"  # File with known gaps in assembly
}
```

To compute neanderthal callability with 1000GP sample we need to calculate $1000GP \intersect (nd1 \union nd2 \union ...)$ 

To speed up preparation remain in 1kG .vcf only samples from .json file

