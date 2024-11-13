#!/bin/bash

CHR=$1
panelfinal=$2


eu=$3
out=$4
arch=$5
aa=$6
bed=$7

outtxt=$8


bcftools query -f '%POS\t%REF\t%ALT\n' ${panelfinal}>positions.chr${CHR}.txt

bcftools query -S ${out} -f '[%GT]\n' ${panelfinal}|awk '{for(i=1;i<=length($0);i++){a[substr($0,i,1)]=1} for(i in a){printf("%s",i)} print "";delete a}'|sed 's/|//g' |sed 's|/||g' |sed -e 's/\.//g'|sed 's/^$/\./' >al.chr${CHR}.outgroup.txt

bcftools query -S ${arch} -f '[%GT]\n' ${panelfinal}|awk '{for(i=1;i<=length($0);i++){a[substr($0,i,1)]=1} for(i in a){printf("%s",i)} print "";delete a}'|sed 's|/||g' |sed 's/|//g'|sed -e 's/\.//g'|sed 's/^$/\./' >al.chr${CHR}.archaic.txt

bcftools query -S ${eu}  -f '[%GT ]\n'  ${panelfinal} |sed  's/|/ /g'|sed 's|/| |g' >  obs.chr${CHR}.ingroup.txt

paste  positions.chr${CHR}.txt al.chr${CHR}.outgroup.txt al.chr${CHR}.archaic.txt obs.chr${CHR}.ingroup.txt> 3.chr${CHR}.txt


printf '#POSITIONS\t#REF\t#ALT\t#OUTGROUP\t#ARCHAIC\t#OBSERVATIONS\n'>header.txt
cat header.txt 3.chr${CHR}.txt > temp.allels.ref.and.obs.chr${CHR}.txt

rm header.txt
paste positions.chr${CHR}.txt obs.chr${CHR}.ingroup.txt > obs.chr${CHR}.txt

rm al.*
rm positions.chr${CHR}.txt 
rm obs.chr${CHR}.ingroup.txt
rm 3.chr${CHR}.txt


python3 obs2.py ${CHR}  temp.allels.ref.and.obs.chr${CHR}.txt ${aa} ${bed} ${outtxt}

rm temp.allels.ref.and.obs.chr${CHR}.txt
rm obs.chr${CHR}.txt














