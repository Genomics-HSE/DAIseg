import os
import json
from pathlib import Path
import subprocess as sbp


def run(json_file):
    if not Path(json_file).exists():
        raise ValueError(f"json file {json_file} not found")

    with open(json_file, "r") as f:
        config = json.load(f)

    # make sure prefix exists
    prefix = Path(config["prefix"])
    os.makedirs(prefix, exist_ok=True)

    # create merged mask
    strict_mask = Path(config["files"]["1000GP_files"]["bed"])
    tmp_strict_mask_sorted = prefix / "tmp_strict_mask.bed"
    strict_mask_cmd = f"bedtools sort -i {strict_mask} >{tmp_strict_mask_sorted}"
    sbp.run(strict_mask_cmd, shell=True)

    tmp_out = prefix / "tmp_intersect_mask.bed"
    merge_cmd = f"""
    jq -r '.files.neand_files[].bed' {json_file} | sed 's|^~|$HOME|' |
    xargs -I _ zcat -f _ | bedtools sort -i - |
    bedtools merge -i - |
    bedtools intersect -sorted -a {tmp_strict_mask_sorted}  -b - >{tmp_out}
    """

    sbp.run(merge_cmd, shell=True)

    # filter the sizes file
    sizes = Path(config["files"]["chr_lengths"])
    tmp_sizes = prefix / "temp.sizes"

    with open(sizes, "r") as f:
        target_chr_size = [x for x in f.readlines() if x.startswith(config["CHROM"])]

    assert len(target_chr_size) == 1
    with open(tmp_sizes, "w") as f:
        f.write(target_chr_size[0])

    # make_windows
    human_out = prefix / config["window_callability"]["Thousand_genomes"]
    nd_out = prefix / config["window_callability"]["Nd_1k_genomes"]

    nd_cmd = f"""
    bedtools makewindows -g {tmp_sizes} -w 1000 | bedtools coverage -a stdin -b {tmp_out} > {nd_out}
    """

    sbp.run(nd_cmd, shell=True)

    human_cmd = f"""
    bedtools makewindows -g {tmp_sizes} -w 1000 | bedtools coverage -a stdin -b {strict_mask} > {human_out}
    """

    sbp.run(human_cmd, shell=True)

    tmp_out.unlink()
    tmp_strict_mask_sorted.unlink()
    tmp_sizes.unlink()
