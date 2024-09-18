#!/bin/bash

mkdir Ancestral.Alleles
#we obtain Ancestral allele file from 1000GP hg19

CHAIN=/media/scglab/T7/Work/data/liftover/hg19ToHg38.over.chain.gz
dir=Ancestral.Alleles







for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do


VCF=/media/scglab/T7/Work/data/1000GP/${i}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
FASTA=/media/scglab/T7/Work/data/grch38.fasta/chr${i}.fa




bcftools view -s HG00096 ${VCF} -Oz -o temp.chr${i}.vcf.gz 
tabix -p vcf temp.chr${i}.vcf.gz 

CrossMap vcf ${CHAIN} temp.chr${i}.vcf.gz ${FASTA}  temp2.chr${i}.grch38.vcf --no-comp-allele 

bcftools query -f '%POS %REF %ALT %INFO\n' temp2.chr${i}.grch38.vcf> ./${dir}/POS.REF.ALT.INFO.chr${i}.txt



rm temp.chr${i}.vcf.gz
rm temp.chr${i}.vcf.gz.tbi
rm temp2.chr${i}.grch38.vcf
rm temp2.chr${i}.grch38.vcf.unmap

python3 Ancestral.Alleles.py ${i}

rm ./${dir}/POS.REF.ALT.INFO.chr${i}.txt

echo "Done"
done




