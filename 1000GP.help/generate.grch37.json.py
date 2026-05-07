#!/usr/bin/env python3

import json
import os
import sys

PED_DIR = "/home/share/human.data/1000GP/predigree.grch37"

def read_samples(pop):
    with open(os.path.join(PED_DIR, f"{pop}.unrelated.ids")) as f:
        samples = [line.strip() for line in f if line.strip()]
        samples.sort()  # чтобы порядок был стабильным (опционально)
        return samples

def get_base_config(pop_name, chrom, out_abs, yri_samples):
    return {
        "description": f"DAIseg.simple configuration to run for {pop_name}",
        "CHROM": str(chrom),
        "output": f"{pop_name}.YRI.grch37.chr{chrom}",
        "prefix": os.path.join(out_abs, f"{pop_name}.YRI"),
        "files": {
            "neand_files": {
                "Vindija33.19": {
                    "bed": f"/home/share/human.data/neand/33.19/bed/chr{chrom}_mask.bed.gz",
                    "vcf": f"/home/share/human.data/neand/33.19/chr{chrom}_mq25_mapab100.vcf.gz"
                },
                "Altai": {
                    "bed": f"/home/share/human.data/neand/altai/bed/chr{chrom}_mask.bed.gz",
                    "vcf": f"/home/share/human.data/neand/altai/chr{chrom}_mq25_mapab100.vcf.gz"
                },
                "Chagyrskaya-Phalanx": {
                    "bed": f"/home/share/human.data/neand/Chagyrskaya/bed/chr{chrom}_mask.bed.gz",
                    "vcf": f"/home/share/human.data/neand/Chagyrskaya/chr{chrom}.noRB.vcf.gz"
                }
            },
            "1000GP_files": {
                "bed": f"/home/share/human.data/1000GP/1000GP.grch37/bed/chr{chrom}.renamed.bed",
                "vcf": f"1kG_filtered.chr{chrom}.grch37.vcf.gz",
                "vcf_initial": f"/home/share/human.data/1000GP/1000GP.grch37/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
            },
            "ancestral": {"fasta": f"/home/share/human.data/Anc.fa/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_{chrom}.fa"},
            "reference": {"fasta": "/home/share/human.data/ref.fa/hs37.1000GP.fasta/hs37d5.fa.gz"},
            "chr_lengths": "/home/share/human.data/ref.fa/hg19.chr.lengths/hg19.chrom.len"
        },
        "samples": {
            "outgroup": yri_samples,
            "ingroup": [],
            "neand": ["Vindija33.19", "AltaiNeandertal", "Chagyrskaya-Phalanx"]
        },
        "parameters_initial": {
            "admixture_proportion": 0.02,
            "introgression_time": 55000,
            "rr": 1e-08,
            "mutation": 1.25e-08,
            "window_length": 1000,
            "generation_time": 29,
            "t_archaic_c": 550000,
            "t_split_c": 70000,
            "t_introgression_c": 55000,
            "t_introgression": 55000
        },
        "window_callability": {
            "Thousand_genomes": f"coverage_1kG.chr{chrom}.grch37.bed",
            "Nd_1k_genomes": f"coverage_1kG.nd.chr{chrom}.grch37.bed"
        },
        "data": f"prep.chr{chrom}.grch37.tsv",
        "gaps": "/home/share/human.data/ref.fa/gaps.grch37/gap.renamed.txt"
    }

def main():
    pop_name = sys.argv[1].upper()
    out = sys.argv[2]
    
    # Преобразуем в абсолютный путь
    out_abs = os.path.abspath(out)
    
    yri_samples = read_samples("YRI")
    ingroup_samples = read_samples(pop_name)
    
    output_dir = os.path.join(out_abs, f"{pop_name}.YRI.jsons")
    os.makedirs(output_dir, exist_ok=True)
    
    for chrom in range(1, 23):
        config = get_base_config(pop_name, chrom, out_abs, yri_samples)
        config["samples"]["ingroup"] = ingroup_samples
        
        out_file = os.path.join(output_dir, f"{pop_name}.YRI.grch37.chr{chrom}.json")
        with open(out_file, 'w') as f:
            json.dump(config, f, indent=2)

if __name__ == "__main__":
    main()
