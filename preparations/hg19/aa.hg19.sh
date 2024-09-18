#!/bin/bash


#we obtain ancestral alleles file from 1000GP hg19
dir=Ancestral.Alleles
mkdir Ancestral.Alleles
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do

VCF1000=/media/scglab/T7/Work/data/1000GP/${i}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
bcftools query -f '%POS %REF %ALT %INFO\n' ${VCF1000}> ./${dir}/POS.REF.ALT.INFO.chr${i}.txt

python3 Ancestral.Alleles.py ${i}
rm ./${dir}/POS.REF.ALT.INFO.chr${i}.txt
done





