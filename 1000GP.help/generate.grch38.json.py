#!/usr/bin/env python3

import json
import os
import sys

PED_DIR = "/home/share/human.data/1000GP/pedigree.grch38"

def read_samples(pop):
    """Читает образцы из файла {pop}.unrelated.2504.txt или {pop}.unrelated.hc3202.txt"""
    if pop == "YRI":
        file_path = os.path.join(PED_DIR, f"{pop}.unrelated.2504.txt")
    else:
        file_path = os.path.join(PED_DIR, f"{pop}.unrelated.hc3202.txt")
    
    with open(file_path) as f:
        return [line.strip() for line in f if line.strip()]

def get_base_config(pop_name, chrom, out_abs, yri_samples):
    return {
        "description": f"DAIseg.simple configuration to run for {pop_name} on GRCh38",
        "CHROM": f'chr{chrom}',
        "output": f"{pop_name}.YRI.grch38.chr{chrom}",
        "prefix": os.path.join(out_abs, f"{pop_name}.YRI"),
        "files": {
            "neand_files": {
                "Vindija33.19": {
                    "bed": f"/home/share/human.data/neand/33.19.grch38/bed/chr{chrom}_mask.bed.gz",
                    "vcf": f"/home/share/human.data/neand/33.19.grch38/chr{chrom}.vcf.gz"
                },
                "AltaiNeandertal": {
                    "bed": f"/home/share/human.data/neand/altai.grch38/bed/chr{chrom}_mask.bed.gz",
                    "vcf": f"/home/share/human.data/neand/altai.grch38/chr{chrom}.vcf.gz"
                },
                "Chagyrskaya-Phalanx": {
                    "bed": f"/home/share/human.data/neand/Chagyrskaya.grch38/bed/chr{chrom}_mask.bed.gz",
                    "vcf": f"/home/share/human.data/neand/Chagyrskaya.grch38/chr{chrom}.vcf.gz"
                }
            },
            "1000GP_files": {
                "bed": f"/home/share/human.data/1000GP/1000GP.grch38/bed/chr{chrom}.clean.bed",
                "vcf": f"1kG_filtered.chr{chrom}.grch38.vcf.gz",
                "vcf_initial": f"/home/share/human.data/1000GP/1000GP.grch38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"
            },
            "ancestral": {
                "fasta": f"/home/share/human.data/Anc.fa/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_{chrom}.fa"
            },
            "reference": {
                "fasta": "/home/share/human.data/ref.fa/grch38.fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa"
            },
            "chr_lengths": "/home/share/human.data/ref.fa/grch38.lengths/hg38.chrom.sizes"
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
            "Thousand_genomes": f"coverage_1kG.chr{chrom}.grch38.bed",
            "Nd_1k_genomes": f"coverage_1kG.nd.chr{chrom}.grch38.bed"
        },
        "data": f"prep.chr{chrom}.grch38.tsv",
        "gaps": "/home/share/human.data/ref.fa/gaps.grch38/gap.txt"
    }

def main():
    pop_name = sys.argv[1].upper()
    out = sys.argv[2]
    
    out_abs = os.path.abspath(out)
    
    yri_samples = read_samples("YRI")
    ingroup_samples = read_samples(pop_name)
    
    output_dir = os.path.join(out_abs, f"{pop_name}.YRI.jsons")
    os.makedirs(output_dir, exist_ok=True)
    
    for chrom in range(1, 23):
        config = get_base_config(pop_name, chrom, out_abs, yri_samples)
        config["samples"]["ingroup"] = ingroup_samples
        
        out_file = os.path.join(output_dir, f"{pop_name}.YRI.grch38.chr{chrom}.json")
        with open(out_file, 'w') as f:
            json.dump(config, f, indent=2)

if __name__ == "__main__":
    main()
