#!/bin/bash

dir=$1
CHR=$2
bed=$3
n1=$4
n2=$5
n3=$6
GP1000=$7

obs=$8
out=$9
arch=${10}

outfilevcf=${11}
outtxt=${12}
aa=${13}



cd preparations


./archaic.covering.sh ${CHR} ${bed} ${n1} ${n2} ${n3} ${dir}



./new.panel.preparation.Linux.sh ${CHR} ${out} ${obs} ${bed} ${GP1000} ${n1} ${n2} ${n3} ${outfilevcf}


./new.make.obs.sh ${CHR} ${outfilevcf} ${obs} ${out} ${arch}  ${aa} ${bed} ${outtxt}



#./full.preparation.Linux.sh hg19 22 /media/scglab/T7/Work/18.09.eu/preparations/hg19/beds/chr22.hg19.bed /media/scglab/T7/Work/data/neand/33.19/22/chr22_mq25_mapab100.vcf.gz /media/scglab/T7/Work/data/neand/altai/22/chr22_mq25_mapab100.vcf.gz /media/scglab/T7/Work/data/neand/ChagyrskayaOkladnikov/split.22.vcf.gz /media/scglab/T7/Work/data/1000GP/22/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz /media/scglab/T7/Work/18.09.eu/samples/obs.samples.txt /media/scglab/T7/Work/18.09.eu/samples/Outgroup.txt /media/scglab/T7/Work/18.09.eu/samples/archaic.txt hg19.chr22.vcf.gz all.chr22.hg19.txt /media/scglab/T7/Work/data/hg19.Ancestral.Alleles/hg19.AA.chr22.txt

#./full.preparation.Linux.sh grch38 22 /media/scglab/T7/Work/18.09.eu/preparations/grch38/beds/chr22.grch38.bed /media/scglab/T7/Work/data/neand/33.19.grch38/33.19.chr22.grch38.vcf.gz /media/scglab/T7/Work/data/neand/altai.grch38/altai.chr22.grch38.vcf.gz /media/scglab/T7/Work/data/neand/ChagyrskayaOkladnikov.grch38/okl.grch38.22.vcf.gz /media/scglab/T7/Work/data/1000GP.grch38/22/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz /media/scglab/T7/Work/18.09.eu/samples/obs.samples.txt /media/scglab/T7/Work/18.09.eu/samples/Outgroup.txt /media/scglab/T7/Work/18.09.eu/samples/archaic.txt grch38.chr22.vcf.gz all.chr22.grch38.txt /media/scglab/T7/Work/data/grch38.Ancestral.Alleles/grch38.AA.chr22.txt

