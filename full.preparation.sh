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


./archaic.covering.sh ${CHR} ${bed} ${n1} ${n2} ${n3}



./new.panel.preparation.Linux.sh ${CHR} ${out} ${obs} ${bed} ${GP1000} ${n1} ${n2} ${n3} ${outfilevcf}


./new.make.obs.sh ${CHR} ${outfilevcf} ${obs} ${out} ${arch}  ${aa} ${bed} ${outtxt}




