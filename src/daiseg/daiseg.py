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

    parser = argparse.ArgumentParser(
        description=script_description.strip(), 
        formatter_class=argparse.RawTextHelpFormatter
    )

    subparser = parser.add_subparsers(dest='mode')
 
    # 1. Standard Run
    decode_parser = subparser.add_parser('run', help='Run HMM (Standard)')
    decode_parser.add_argument("-json", help="File with parameters", type=str, nargs='+', required=True)


    # 2. Run with EM 
    decode_subparser = subparser.add_parser('run_EM', help='Run Global EM training and Inference')
    decode_subparser.add_argument("-threads", help="Number of threads (optional)", type=int, required=False, default=4)
    decode_subparser.add_argument("-jsons", help="List of JSON files (e.g. sims/*.json)", nargs='+', required=True)
    decode_subparser.add_argument("-out", help="Path to save merged results (e.g. all.results.tsv)", type=str, required=False)

    # 2.1 Run with EM + Matrix Update 
    decode_parser = subparser.add_parser('run_EM_trans', help='Run EM optimizing Transitions ')
    decode_parser.add_argument("-threads", help="Threads", type=int, required=False, default=4)
    decode_parser.add_argument("-jsons", help="List of JSON files", nargs='+', required=True)
    decode_parser.add_argument("-out", help="Path to save results", type=str)
    decode_parser.add_argument("-iter", help="Max iterations", type=int, default=10)

    # 3. Helpers
    decode_subparser = subparser.add_parser('restrict_1kG', help='Helper')
    decode_subparser.add_argument("-threads", type=int, required=False, default=1)
    decode_subparser.add_argument("-json", type=str, required=True)

    decode_subparser = subparser.add_parser('callability', help='Helper')
    decode_subparser.add_argument("-threads", type=int, required=False, default=1)
    decode_subparser.add_argument("-json", type=str, required=True)

    decode_subparser = subparser.add_parser('prep', help='Helper')
    decode_subparser.add_argument("-threads", type=int, required=False, default=1)
    decode_subparser.add_argument("-json", type=str, required=True)

    args = parser.parse_args()

    # get scripts dir as a subling to python source
    scripts_dir = Path(__file__).parent.parent / "scripts"

    if args.mode == 'run':
        hmm.run_daiseg(args.json)

    elif args.mode == 'run_EM':
        print(f"Starting Batch EM pipeline for {len(args.jsons)} config files...")
        em_alg.run_batch_em_pipeline(args.jsons, output_combined_file=args.out)

    elif args.mode == 'run_EM_trans':
        print(f"Starting optimized EM pipeline for {len(args.jsons)} config files...")
        em_alg.run_batch_em_pipeline_v2(
            args.jsons,
            output_combined_file=args.out,
            threads=args.threads,
            max_iter=args.iter
        )


    elif args.mode == 'restrict_1kG':
        restrict.run(args.json, args.threads)

        # result = subprocess.run(
        #     [str(scripts_dir / 'extract_samples.sh'), args.json, str(args.threads)],
        #     text=True,
        #     check=True
        # )

    elif args.mode == 'callability':
        callability.run(args.json)

        # result = subprocess.run(
        #     [str(scripts_dir / 'callability.sh'), args.json],
        #     text=True, check=True
        # )

    elif args.mode == 'prep':
        prep.run(args.json)

if __name__ == "__main__":
    main()

