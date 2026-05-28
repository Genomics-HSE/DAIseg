import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import daiseg.hmm as hmm
import daiseg.em_alg as em_alg
import daiseg.prep as prep
import daiseg.restrict as restrict
import daiseg.callability as callability


def main():
    script_description = """
    DAIseg: inference of introgressed archaic segments.
    """

    cmd = argparse.ArgumentParser(
        description=script_description.strip(),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    mode = cmd.add_subparsers(dest="mode")

    # 1. Standard Run
    cmd_run = mode.add_parser("run", help="Run HMM (Standard)")
    cmd_run.add_argument(
        "-json", help="File with parameters", type=str, nargs="+", required=True
    )

    # 2. Run with EM
    cmd_run_EM = mode.add_parser("run_EM", help="Run Global EM training and Inference")
    cmd_run_EM.add_argument(
        "-threads",
        help="Number of threads (optional)",
        type=int,
        required=False,
        default=4,
    )
    cmd_run_EM.add_argument(
        "-jsons", help="List of JSON files (e.g. sims/*.json)", nargs="+", required=True
    )
    cmd_run_EM.add_argument(
        "-out",
        help="Path to save merged results (e.g. all.results.tsv)",
        type=str,
        required=False,
    )

    # 2.1 Run with EM + Matrix Update
    cmd_run_trans = mode.add_parser(
        "run_EM_trans", help="Run EM optimizing Transitions "
    )
    cmd_run_trans.add_argument(
        "-threads", help="Threads", type=int, required=False, default=4
    )
    cmd_run_trans.add_argument(
        "-jsons", help="List of JSON files", nargs="+", required=True
    )
    cmd_run_trans.add_argument("-out", help="Path to save results", type=str)
    cmd_run_trans.add_argument("-iter", help="Max iterations", type=int, default=10)

    # 3. Data prep helpers
    cmd_prep = mode.add_parser("prep", help="Helper")
    cmd_prep.add_argument("-threads", type=int, required=False, default=1)
    cmd_prep.add_argument("-json", type=str, required=True)

    args = cmd.parse_args()

    if args.mode == "run":
        hmm.run_daiseg(args.json)

    elif args.mode == "run_EM":
        print(f"Starting Batch EM pipeline for {len(args.jsons)} config files...")
        em_alg.run_batch_em_pipeline(args.jsons, output_combined_file=args.out)

    elif args.mode == "run_EM_trans":
        print(f"Starting optimized EM pipeline for {len(args.jsons)} config files...")
        em_alg.run_batch_em_pipeline_v2(
            args.jsons,
            output_combined_file=args.out,
            threads=args.threads,
            max_iter=args.iter,
        )

    elif args.mode == "prep":
        restrict.run(args.json, args.threads)
        callability.run(args.json)
        prep.run(args.json)


if __name__ == "__main__":
    main()
