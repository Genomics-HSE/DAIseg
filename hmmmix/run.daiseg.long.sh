#!/bin/bash

for c in {1..2};
do
python3 dai.seg.py --obs_samples samples.txt --bed test.bed --HMM_par par.file.txt --EM no --prepared_file obs.chr${c}.txt --o out.chr${c} --arch_cover arch.cover.txt --decoding viterbi
done


