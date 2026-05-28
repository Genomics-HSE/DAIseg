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
    tmp_out = prefix / "tmp_intersect_mask.bed"

    merge_cmd = f"""
    jq -r '.files.neand_files[].bed' {json_file} | sed 's|^~|$HOME|' |
    xargs -I _ zcat -f _ | bedtools sort -i - |
    bedtools merge -i - |
    bedtools intersect -sorted -a <(bedtools sort -i {strict_mask})  -b - >{tmp_out}
    """

    sbp.run(merge_cmd, shell=True)

    # make_windows
    sizes = Path(config["files"]["chr_lengths"])
    human_out = prefix / config["window_callability"]["Thousand_genomes"]
    nd_out = prefix / config["window_callability"]["Nd_1k_genomes"]

    nd_cmd = f"""
    bedtools makewindows -g {sizes} -w 1000 | bedtools coverage -a stdin -b {tmp_out} > {nd_out}
    """

    sbp.run(nd_cmd, shell=True)

    human_cmd = f"""
    bedtools makewindows -g {sizes} -w 1000 | bedtools coverage -a stdin -b {strict_mask} > {human_out}
    """

    sbp.run(human_cmd, shell=True)

    tmp_out.unlink()
