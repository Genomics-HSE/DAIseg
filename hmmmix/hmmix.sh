#!/bin/bash

CHR=2

A='ts.chr1.vcf'
A1='ancestral.chr1.fa'

for c in {2..2};
do
B="ts.chr${c}.vcf"
A="${A},${B}"
B1="ancestral.chr${c}.fa"
A1="${A1},${B1}"
done

prefix=./hmmix

hmmix create_outgroup -ind=individuals.json -vcf=${A} -out=outgroup.txt -ancestral ${A1}
 
hmmix mutation_rate -outgroup=outgroup.txt -window_size=10000000 -out mutrationrate.bed

hmmix create_ingroup -ind=individuals.json -vcf=${A} -out=obs -outgroup=outgroup.txt -ancestral ${A1} 


hmmix train -obs=obs.tsk_0.txt -mutrates=mutrationrate.bed -param=Initialguesses.json -out=trained.tsk_0.json -haploid
hmmix decode -obs=obs.tsk_0.txt -mutrates=mutrationrate.bed -param=trained.tsk_0.json -out=out.tsk_0  -haploid


#rm outgroup.txt
rm mutrationrate.bed
#rm obs.tsk_0.txt
rm trained.tsk_0.json

