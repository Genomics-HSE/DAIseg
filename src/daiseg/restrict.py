import os
import json
from pathlib import Path
import subprocess as sbp


def run(json_file, threads):
    # check json exists
    if not Path(json_file).exists():
        raise ValueError(f"json file {json_file} not found")

    with open(json_file, "r") as f:
        config = json.load(f)

    prefix = Path(config["prefix"])
    os.makedirs(prefix, exist_ok=True)

    file_1kG = Path(config["files"]["1000GP_files"]["vcf_initial"])

    if not file_1kG.exists():
        raise ValueError(f"variant file {file_1kG} not found")

    chrom = config["CHROM"].removeprefix("chr")
    out_vcf = prefix / config["files"]["1000GP_files"]["vcf"]

    samples = config["samples"]["ingroup"] + config["samples"]["outgroup"]
    with open(prefix / "samples.txt", "w") as f:
        for s in samples:
            f.write(f"{s}\n")

    cmd = f"""
    bcftools view --threads {threads} -S {prefix / "samples.txt"} --force-samples --trim-alt-alleles -Ou {file_1kG} |
    bcftools norm --threads {threads} -m -any -Ou |
    bcftools view --threads {threads} -v snps -Ou |
    bcftools norm --threads {threads} -m +any -Ou |
    bcftools view -Wtbi --threads {threads} -m2 -M4 -Oz -o {out_vcf}
    """

    sbp.run(cmd, shell=True)

    (prefix / "samples.txt").unlink()
