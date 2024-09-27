#!/bin/bash

#change the following names and directories
CHR=$1
bed=$4
NAME1000=$5
n1=$6
n2=$7
n3=$8
panel=$9

cat $3 $2 > samples.for.hmm.txt

echo "DAIseg: Extracting the positions of 1000GP..."
bcftools query -R ${bed} -f '%CHROM\t%POS\n' ${NAME1000} > 1000GP.pos.chr${CHR}.txt


echo "DAIseg: Extracts non-1000GP positions that has alternates in archaic..."
bcftools query -T ^1000GP.pos.chr${CHR}.txt -R ${bed} -i 'ALT!="."' -f '%CHROM\t%POS\n' ${n1} > 1.txt
bcftools query -T ^1000GP.pos.chr${CHR}.txt -R ${bed} -i 'ALT!="."' -f '%CHROM\t%POS\n' ${n2} > 2.txt
bcftools query -T ^1000GP.pos.chr${CHR}.txt -R ${bed} -i 'ALT!="."' -f '%CHROM\t%POS\n' ${n3} > 3.txt

#join extra positions
cat 1.txt 2.txt 3.txt| sort -u >extra.pos.chr${CHR}.txt

echo "DAIseg: join extra positions with 1000GP positions..."
cat  1000GP.pos.chr${CHR}.txt extra.pos.chr${CHR}.txt|sort -u> positions.chr${CHR}.txt


#restrict 
bcftools view -R ${bed} -T positions.chr${CHR}.txt ${n1} -Oz -o filtered.snps.chr${CHR}.1.vcf.gz
bcftools view -R ${bed} -T positions.chr${CHR}.txt ${n2} -Oz -o filtered.snps.chr${CHR}.2.vcf.gz
bcftools view -R ${bed} -T positions.chr${CHR}.txt ${n3} -Oz -o filtered.snps.chr${CHR}.3.vcf.gz
tabix -p vcf filtered.snps.chr${CHR}.1.vcf.gz
tabix -p vcf filtered.snps.chr${CHR}.2.vcf.gz
tabix -p vcf filtered.snps.chr${CHR}.3.vcf.gz



echo "DAIseg: we are gluing archaic genomes in extra positions..."
bcftools merge filtered.snps.chr${CHR}.1.vcf.gz filtered.snps.chr${CHR}.2.vcf.gz filtered.snps.chr${CHR}.3.vcf.gz  -Oz -o merged.chr${CHR}.archaic.vcf.gz
tabix -p vcf merged.chr${CHR}.archaic.vcf.gz





echo "DAIseg: Extracting positions from merged archaic..."
bcftools query -R ${bed} -f '%CHROM\t%POS\n'  merged.chr${CHR}.archaic.vcf.gz > archaic.merged.chr${CHR}.txt #we need to do it due to the possibility of abscence of some variants


echo "DAIseg: Restricting on the existing positions of merged archaic vcf"
bcftools view -R ${bed} -v snps -S samples.for.hmm.txt -T archaic.merged.chr${CHR}.txt ${NAME1000} -Oz -o temp.1000.chr${CHR}.vcf.gz
tabix -p vcf temp.1000.chr${CHR}.vcf.gz


echo "DAIseg: we are glueing archaic genomes and 1000GP in the existing positions of archaic using 0/0 genotypes on the missing 1000GP positiions"
bcftools merge -0 merged.chr${CHR}.archaic.vcf.gz temp.1000.chr${CHR}.vcf.gz -Ou|bcftools view -v snps|bcftools annotate -x INFO -Oz -o 1.vcf.gz
tabix -p vcf 1.vcf.gz






echo "DAIseg: The next goal is to obtain vcf with positions of 1000GP which are missed in archaic samples"
ggrep -Fxf 1000GP.pos.chr${CHR}.txt archaic.merged.chr${CHR}.txt > int.txt
bcftools view -R ${bed} -v snps -T ^int.txt -S samples.for.hmm.txt ${NAME1000} -Oz -o 2.vcf.gz
tabix -p vcf 2.vcf.gz

echo "DAIseg: Extaract empty vcf with header only to merge it with vcf obtained on the previous step"
bcftools view -h merged.chr${CHR}.archaic.vcf.gz -Oz -o header.archaic.vcf.gz
tabix -p vcf header.archaic.vcf.gz


echo "DAIseg: Merging"
bcftools merge header.archaic.vcf.gz 2.vcf.gz| bcftools annotate -x INFO  -Oz -o 3.vcf.gz
tabix -p vcf 3.vcf.gz




echo "DAIseg: Concatenating"
bcftools concat 3.vcf.gz 1.vcf.gz|bcftools sort |bcftools view -v snps| bcftools annotate -x INFO  -Oz -o ${9}
tabix -p vcf ${9} 

rm extra.pos.chr${CHR}.txt
rm int.txt

rm header.archaic.*
rm filtered.snps.*
rm 1.* 2.* 3.*
rm merged.chr${CHR}.archaic.vcf.*
rm temp.1000.chr${CHR}.*
rm samples.for.hmm.txt
rm archaic.merged.chr${CHR}.txt
rm 1000GP.pos.chr${CHR}.txt

