# Working with real data

## Download files 

### Genomics masks for callability 

1. wget  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.strict_mask.autosomes.bed for 1kG project and
   split it by chr 
       for chr in {1..22}; do
           grep -w "chr${chr}" 20140520.strict_mask.autosomes.bed > "chr${chr}.bed"
       done
2. http://ftp.eva.mpg.de/neandertal/Vindija/FilterBed/Vindija33.19/chr22_mask.bed.gz for Vindija
3. http://ftp.eva.mpg.de/neandertal/Vindija/FilterBed/Altai/chr22_mask.bed.gz for Altai neanderthal

### Chromosome sizes to make 
Download chr-length https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes | sed 's/^chr//' hg19.chrom.sizes>hg19.chrom.len


To compute neanderthal callability with 1000GP sample we need to calculate $1000GP \intersect (nd1 \union nd2 \union ...)$ 

To speed up preparation remain in 1kG .vcf only samples from .json file



Download hg19 gaps
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz
zcat gap.txt.gz > gap.txt



## .json input description

Main configuration file of DAIseg.simple:

```json
{
  "description": "DAISeg.simple configuration to run ",
  "CHROM": "22", #chr name
  "prefix": "out.chr22",
  "files": {
    "neand_files": {
      "Vindija33.19": {
        "bed": "/path/full/to/bed",
        "vcf": "/path/full/to/vcf"
      },
      "Altai": {
        "bed": "/path/full/to/bed",
        "vcf": "/path/full/to/vcf"
      }
    },
    "1000GP_files": {
      "bed": "/path/full/to/bed",
      "vcf": "/path/full/to/vcf"
    },
    "ancestral": {
      "fasta": "/path/full/to/homo_sapiens_ancestor_22.fa"
    },
    "reference": {
      "fasta": "/path/full/to/hs37d5.fa.gz"
    },
    "chr_lengths": "/path/full/to/hg19.chrom.len"
  },
  "samples": {
    "outgroup": [
      "NA18488", "NA18508", "NA18510", "NA18489", "NA18504",
      "NA18511", "NA18934", "NA18864", "NA18871", "NA18876",
      ...
    ],
    "ingroup": [
      "HG01503", "HG01510", "HG01515", "HG01522", "HG01527",
      "HG01680", "HG01685", "HG01697", "HG01700", "HG01705"
    ],
    "neand": [
      "Vindija33.19",
      "AltaiNeandertal"
    ]
  },
   //simulation parameters
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
    "Thousand_genomes": "/path/to/file",
    "Nd_1k_genomes": "/path/to/file"
  },
  "data": "/path/to/file",
  "gaps": "/path/to/file"
}

}
```
