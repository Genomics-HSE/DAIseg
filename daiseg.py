import argparse
import json
import os
import subprocess
import sys

# Imports of your modules
import hmm
import em_alg

def print_script_usage():
    prnt = f'''
    Script for identifying introgressed archaic segments
    '''
    return prnt


def main():

    parser = argparse.ArgumentParser(description=print_script_usage(), formatter_class=argparse.RawTextHelpFormatter)
    subparser = parser.add_subparsers(dest='mode')
    
    # 1. Standard Run
    decode_subparser = subparser.add_parser('run', help='Run HMM (Standard)')
    decode_subparser.add_argument("-json", help="File with parameters", type=str, nargs='+', required=True)
    
    # 2. Run with EM (Batch Mode)
    decode_subparser = subparser.add_parser('run.with.EM', help='Run Global EM training and Inference')
    decode_subparser.add_argument("-threads", help="Number of threads (optional)", type=int)
    decode_subparser.add_argument("-jsons", help="List of JSON files (e.g. sims/*.json)", nargs='+', required=True)
    decode_subparser.add_argument("-out", help="Path to save merged results (e.g. all.results.tsv)", type=str, required=False)
    
    # 3. Helpers
    decode_subparser = subparser.add_parser('restrict_1kG', help='Helper')
    decode_subparser.add_argument("-threads", type=int, required=True)
    decode_subparser.add_argument("-json", type=str, required=True)
 
    decode_subparser = subparser.add_parser('callability', help='Helper')
    decode_subparser.add_argument("-threads", type=int, required=True)
    decode_subparser.add_argument("-json", type=str, required=True)

    decode_subparser = subparser.add_parser('main.prep', help='Helper')
    decode_subparser.add_argument("-threads", type=int, required=True)
    decode_subparser.add_argument("-json", type=str, required=True)
    
    args = parser.parse_args()

    if args.mode == 'run':
        hmm.run_daiseg(args.json)
                
    elif args.mode == 'run.with.EM':
        print(f"Starting Batch EM pipeline on {len(args.jsons)} files...")
        em_alg.run_batch_em_pipeline(args.jsons, output_combined_file=args.out)

    elif args.mode == 'restrict_1kG':
        result = subprocess.run(
            ['bash', 'extract.samples.sh', args.json, str(args.threads)],
            capture_output=True,
            text=True,
            check=True
        )
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)

    elif args.mode == 'callability':
        result = subprocess.run(
            ['bash', 'callability.sh', args.json, str(args.threads)],
            capture_output=True,
            text=True
        )
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
        
    elif args.mode == 'main.prep':
        with open(args.json, 'r') as f:
            jsn = json.load(f)
            
        filename = jsn.get("data", f"prep.chr{jsn['CHROM']}.tsv")
        output_file = os.path.join(jsn["prefix"], filename)

        print(f"Running pipeline... Target Output: {output_file}")
        
        subprocess.run(
            ['python', '-u', 'main.prep.py', args.json, str(args.threads)],
            text=True,
            check=True
        )
        print("Done.")

if __name__ == "__main__":
    main()

