#!/usr/bin/env bash

POP=$1

PED_FILE="20130606_g1k.ped"
VCF="/home/share/human.data/1000GP/1000GP.grch37/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

bcftools query -l "${VCF}" | sort > "grch37.samples.txt"

if [ "$POP" = "IBS" ]; then
    awk -v p="$POP" '$7==p && ($8=="mother" || $8=="father") {print $2}' "${PED_FILE}" \
    | sort \
    | comm -12 - "grch37.samples.txt" \
    > "${POP}.unrelated.grch37.ids"
else
    awk -v p="$POP" '$7==p && $8=="unrel" {print $2}' "${PED_FILE}" \
    | sort \
    | comm -12 - "grch37.samples.txt" \
    > "${POP}.unrelated.grch37.ids"
fi
