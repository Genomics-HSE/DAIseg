#!/bin/bash

dir=$1
CHR=$2
bed=$3
n1=$4
n2=$5
n3=$6
GP1000=$7

out=$8
obs=$9
arch=${10}
outfile=${11}



cd $1


./archaic.covering.sh ${CHR} ${bed} ${n1} ${n2} ${n3}


./ancestral.alleles.sh ${CHR} ${GP1000}

./new.panel.preparation.Linux.sh ${CHR} ${out} ${obs} ${bed} ${GP1000} ${n1} ${n2} ${n3} ${outfile}


./new.make.obs.sh ${CHR} ${outfile} ${obs} ${out} ${arch}  ./Ancestral.Alleles/hg19.AA.chr${CHR}.txt ./regions/chr${CHR}.hg19.bed

cd ../


