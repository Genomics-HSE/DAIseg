#!/usr/bin/env bash

# download data for an example hg19 run

mkdir -p data/{nd,1000g,reference,ancestral,gaps,sizes}
mkdir -p data/nd/{vindija,altai}

# download GRCh37 reference genome (hs37d5)
wget -P data/reference ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

# download ancestral genome (Ensembl release 71)
wget -P data/ancestral ftp://ftp.ensembl.org/pub/release-71/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2
ouch decompress --dir data/ancestral --yes data/ancestral/homo_sapiens_ancestor_GRCh37_e71.tar.bz2

# download chromosome sizes for hg19
wget -P data/sizes https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
# remove 'chr' prefix from chromosome names
sed -i -e 's/^chr//' data/sizes/hg19.chrom.sizes

# download hg19 gap annotations
wget -P data/gaps http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz
ouch decompress --dir data/gaps --yes data/gaps/gap.txt.gz

# Download 1000g phase 3 VCF for chromosome 22
wget -P data/1000g http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
wget -P data/1000g http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

# download the 1000g strict mask
wget -P data/1000g http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.strict_mask.autosomes.bed
# get only the chr22 data
grep -w "chr22" data/1000g/20140520.strict_mask.autosomes.bed >data/1000g/20140520.strict_mask.chr22.bed

# download Vindija Neanderthal VCF
wget -P data/nd/vindija http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr22_mq25_mapab100.vcf.gz
wget -P data/nd/vindija http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/chr22_mq25_mapab100.vcf.gz.tbi
# Download Vindija 33.19 mask for chr22
wget -P data/nd/vindija http://ftp.eva.mpg.de/neandertal/Vindija/FilterBed/Vindija33.19/chr22_mask.bed.gz

# download Altai Neanderthal VCF
wget -P data/nd/altai http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr22_mq25_mapab100.vcf.gz
wget -P data/nd/altai http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr22_mq25_mapab100.vcf.gz.tbi
# download Altai Neanderthal mask for chr22
wget -P data/nd/altai http://ftp.eva.mpg.de/neandertal/Vindija/FilterBed/Altai/chr22_mask.bed.gz
