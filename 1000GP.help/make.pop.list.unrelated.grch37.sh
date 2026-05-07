#!/usr/bin/env bash
set -euo pipefail

POP=$1

PED_FILE="/home/share/human.data/1000GP/pedigree.grch38/integrated_call_samples_v3.20200731.ALL.ped"
VCF="/home/share/human.data/1000GP/1000GP.grch38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz"
OUTDIR="/home/share/human.data/1000GP/pedigree.grch38"

cd "$OUTDIR"

bcftools query -l "$VCF" | sort > hc3202.samples.txt

if [ "$POP" = "IBS" ]; then
    awk -v p="$POP" '
        NR==1 {next}
        $7==p && ($8=="mother" || $8=="father") {print $2}
    ' "$PED_FILE" | sort > "${POP}.unrel.from_ped.txt"
else
    awk -v p="$POP" '
        NR==1 {next}
        $7==p && $8=="unrel" {print $2}
    ' "$PED_FILE" | sort > "${POP}.unrel.from_ped.txt"
fi

comm -12 "${POP}.unrel.from_ped.txt" hc3202.samples.txt > "${POP}.unrelated.hc3202.txt"

wc -l "${POP}.unrelated.hc3202.txt"
