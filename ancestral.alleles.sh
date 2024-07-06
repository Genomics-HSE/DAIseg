#!/bin/bash


dir=Ancestral.Alleles


mkdir Ancestral.Alleles
for i in 22 
do
bcftools query -f '%POS %REF %ALT %INFO\n' /media/scglab/T7/Work/data/1000GP/${i}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz> ./${dir}/POS.REF.ALT.INFO.chr${i}.txt
done

python3 Ancestral.Alleles.py

rm ./${dir}/POS.REF.ALT.INFO.chr${i}.txt
