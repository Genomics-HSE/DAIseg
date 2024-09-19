#!/bin/bash

#change the following names and directories
CHR=$1
bed=$2
n1=$3
n2=$4
n3=$5

genome=${6}



#n1=/media/scglab/T7/Work/data/neand/33.19/chr${CHR}_mq25_mapab100.vcf.gz
#n2=/media/scglab/T7/Work/data/neand/altai/chr${CHR}_mq25_mapab100.vcf.gz
#n3=/media/scglab/T7/Work/data/neand/ChagyrskayaOkladnikov/split.${CHR}.vcf.gz


bcftools query -R ${bed} -f '%POS\n' ${n1} > 1.txt
bcftools query -R ${bed} -f '%POS\n' ${n2} > 2.txt
bcftools query -R ${bed} -f '%POS\n' ${n3} > 3.txt
cat 1.txt 2.txt 3.txt|sort -u > chr${CHR}.archaic.txt

rm 1.txt 2.txt 3.txt


python3 archaic.covering.py ${bed} chr${CHR}.archaic.txt ${CHR} ${genome}
rm chr${CHR}.archaic.txt
