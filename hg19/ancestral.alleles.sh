#!/bin/bash

CHR=$1
VCF1000=$2
dir=Ancestral.Alleles


mkdir Ancestral.Alleles
for i in ${CHR} 
do
bcftools query -f '%POS %REF %ALT %INFO\n' ${VCF1000}> ./${dir}/POS.REF.ALT.INFO.chr${i}.txt
done

python3 Ancestral.Alleles.py ${CHR}

rm ./${dir}/POS.REF.ALT.INFO.chr${i}.txt

#/media/scglab/T7/Work/data/1000GP/${i}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
