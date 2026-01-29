import os
import sys
import subprocess
import time

import pysam

import preprocessing as utils


def run_step(cmd, desc):
    """Run shell command, exit on error."""
    try:
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    except subprocess.CalledProcessError:
        sys.exit(f"[ERROR] Failed at step: {desc}")


def main():
    if len(sys.argv) < 2:
        sys.exit(" Usage: python main.prep.py <config.json>")

    start_time = time.time()

    # --- Load Config ---
    cfg = utils.load_config(sys.argv[1])
    chrom_raw = str(cfg["CHROM"])
    files = cfg["files"]
    prefix = cfg["prefix"]

    # Get output filename from JSON or use default
    output_filename = cfg.get("data", f"prep.chr{chrom_raw.lstrip('chr').lstrip('CHR')}.tsv")
    output_file = os.path.join(prefix, output_filename)

    print(f" [INFO] Output will be written to: {output_file}", file=sys.stderr)

    # Normalize chromosome name
    chrom_no_prefix = chrom_raw.lstrip("chr").lstrip("CHR")

    vcf_1kg = prefix + '/' + utils.expand_path(files["1000GP_files"]["vcf"])
    bed_strict = utils.expand_path(files["1000GP_files"]["bed"])
    anc_path = utils.expand_path(files["ancestral"]["fasta"])

    # Debug: Check FASTA file
    if not os.path.exists(anc_path):
        sys.exit(f"[ERROR] Ancestral FASTA not found: {anc_path}")

    # Try to read FASTA headers for debugging
    try:
        with open(anc_path, 'r') as f:
            headers = []
            for _ in range(10):
                line = f.readline()
                if not line:
                    break
                if line.startswith('>'):
                    headers.append(line.strip())

    except Exception as e:
        print(f" [DEBUG] Could not read FASTA headers: {e}", file=sys.stderr)

    temp_files = []

    print("Creating temporary files...", file=sys.stderr)

    # --- STEP 1: Filter 1000 Genomes ---
    tmp_1kg = f"{prefix}/temp_1kg_{chrom_no_prefix}.bcf"
    cmd_1kg = (
        f"bcftools view --threads 4 -r {chrom_raw} -R {bed_strict} "
        f"-O b -o {tmp_1kg} {vcf_1kg}"
    )
    run_step(cmd_1kg, "Filtering 1000 Genomes")
    run_step(f"bcftools index -f {tmp_1kg}", "Indexing 1000G")

    temp_files.append(tmp_1kg)
    temp_files.append(f"{tmp_1kg}.csi")

    # --- STEP 2: Filter Neanderthals (Parallel) ---
    print("Working with Neanderthals (Parallel)...", file=sys.stderr)
    neand_files = files.get("neand_files", {})
    nd_inputs = []
    running_procs = []

    for i, (name, paths) in enumerate(neand_files.items()):
        n_vcf = utils.expand_path(paths["vcf"])
        n_bed = utils.expand_path(paths["bed"])

        tmp_nd = f"{prefix}/temp_nd_{i}_{chrom_no_prefix}.bcf"

        # Check if BED file exists, if not try alternative naming
        bed_arg = ""
        if n_bed:
            if os.path.exists(n_bed):
                bed_arg = f"-T {n_bed}"
            else:
                # Try alternative naming
                bed_dir = os.path.dirname(n_bed)
                bed_base = os.path.basename(n_bed)

                # Try removing 'chr' prefix
                if bed_base.startswith("chr"):
                    alt_bed = os.path.join(bed_dir, bed_base[3:])
                    if os.path.exists(alt_bed):
                        bed_arg = f"-T {alt_bed}"
                        print(f"Using alternative BED: {alt_bed}", file=sys.stderr)

        # Run in background (threads=2 per job)
        cmd_str = (
            f"bcftools view --threads 2 -r {chrom_raw} {bed_arg} -O b -o {tmp_nd} {n_vcf} && "
            f"bcftools index -f {tmp_nd}"
        )

        print(f"Starting job for {name}...", file=sys.stderr)

        proc = subprocess.Popen(cmd_str, shell=True, executable='/bin/bash')
        running_procs.append(proc)

        temp_files.append(tmp_nd)
        temp_files.append(f"{tmp_nd}.csi")
        nd_inputs.append(tmp_nd)

    # Wait for completion
    for proc in running_procs:
        if proc.wait() != 0:
            for f in temp_files:
                if os.path.exists(f):
                    os.remove(f)
            sys.exit("[ERROR] One of the Neanderthal processing jobs failed.")

    # --- STEP 3: Merge Pipeline ---
    files_str = " ".join([tmp_1kg] + nd_inputs)

    # Use %TGT to extract alleles (e.g., A/G) directly
    pipeline_cmd = (
        f"bcftools merge --threads 4 --force-samples --merge all -O u {files_str} "
        f"| bcftools query -H -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%TGT]\\n'"
    )

    # --- STEP 4: Load Ancestral Genome ---
    try:
        af = pysam.FastaFile(anc_path)

        # Get all available chromosome names in FASTA
        available_chroms = af.references
        print(f"[DEBUG] Available chromosomes in ancestral FASTA: {available_chroms}", file=sys.stderr)
        print(f"[DEBUG] Looking for chromosome related to: '{chrom_raw}'", file=sys.stderr)

        # Smart chromosome finding for non-standard FASTA formats
        chrom_seq = None
        chrom_name_used = None

        # Case 1: If only one chromosome in FASTA, use it
        if len(available_chroms) == 1:
            chrom_name_used = available_chroms[0]
            chrom_seq = af.fetch(chrom_name_used)
            print(f"[DEBUG] Using only available chromosome: '{chrom_name_used}'", file=sys.stderr)

        # Case 2: Try to find chromosome by patterns
        else:
            # Create search patterns
            search_terms = [
                f":{chrom_no_prefix}:",           # Look for :21: in the name
                f"chromosome.*{chrom_no_prefix}",  # chromosome something 21
                f"chr{chrom_no_prefix}",          # chr21
                chrom_no_prefix,                  # 21
                f"CHR{chrom_no_prefix}",          # CHR21
            ]

            for chrom_name in available_chroms:
                for term in search_terms:
                    if term in chrom_name:
                        chrom_name_used = chrom_name
                        chrom_seq = af.fetch(chrom_name_used)
                        print(f"[DEBUG] Found chromosome using pattern '{term}': '{chrom_name_used}'", file=sys.stderr)
                        break
                if chrom_seq:
                    break

        # Case 3: If still not found, try exact match or partial match
        if chrom_seq is None:
            for chrom_name in available_chroms:
                # Try exact match first
                if chrom_name == chrom_raw or chrom_name == f"chr{chrom_no_prefix}" or chrom_name == chrom_no_prefix:
                    chrom_name_used = chrom_name
                    chrom_seq = af.fetch(chrom_name_used)
                    print(f"[DEBUG] Found exact match: '{chrom_name_used}'", file=sys.stderr)
                    break

            # Try partial match
            if chrom_seq is None:
                for chrom_name in available_chroms:
                    if chrom_no_prefix in chrom_name or chrom_raw in chrom_name:
                        chrom_name_used = chrom_name
                        chrom_seq = af.fetch(chrom_name_used)
                        print(f"[DEBUG] Found partial match: '{chrom_name_used}'", file=sys.stderr)
                        break

        # Case 4: Final fallback - if nothing found
        if chrom_seq is None:
            error_msg = f"Chromosome '{chrom_raw}' not found in Ancestral FASTA.\n"
            error_msg += f"Available chromosomes: {available_chroms}\n"
            error_msg += f"Searching for patterns containing: '{chrom_no_prefix}'\n"
            error_msg += "Please check your FASTA file format and chromosome naming."
            raise ValueError(error_msg)

        chrom_len = len(chrom_seq)
        print(f"🔍 [DEBUG] Chromosome '{chrom_name_used}' length: {chrom_len} bp", file=sys.stderr)
        af.close()
    except Exception as e:
        for f in temp_files:
            if os.path.exists(f):
                os.remove(f)
        sys.exit(f"[ERROR] Ancestral fasta issue: {e}")

    # --- STEP 5: Run Pipeline and Write to File ---
    process = subprocess.Popen(pipeline_cmd, shell=True, stdout=subprocess.PIPE,
                               text=True, executable='/bin/bash')

    try:
        # Open output file for writing
        with open(output_file, 'w') as out_f:
            header_line = process.stdout.readline().strip()
            if not header_line.startswith("#"):
                raise RuntimeError("Pipeline returned no header.")

            # Clean header (# [1]CHROM -> CHROM)
            raw_headers = header_line.lstrip("#").split('\t')
            headers = []
            for h in raw_headers:
                if "]" in h:
                    h = h.split("]")[-1]
                if ":" in h:
                    h = h.split(":")[0]
                headers.append(h.strip())

            # Map columns (Strict JSON order)
            idx_chrom, idx_pos, idx_ref, idx_alt, cols_yri, cols_ibs, cols_nd = utils.map_columns(headers, cfg["samples"])

            # === Create Output Header ===
            # Split Ingroup into Sample_1 and Sample_2
            ibs_split_headers = []
            for i in cols_ibs:
                name = headers[i]
                ibs_split_headers.append(f"{name}_1")
                ibs_split_headers.append(f"{name}_2")

            ibs_header_str = "\t".join(ibs_split_headers)

            # Write header to file
            out_f.write(f"#CHROM\tPOS\tREF\tALT\tAncestral\tOutgroup\tNeand\t{ibs_header_str}\n")

            row_count = 0
            # === Process Rows ===
            for line in process.stdout:
                parts = line.strip().split('\t')
                if len(parts) != len(headers):
                    continue

                try:
                    ref = parts[idx_ref]
                    alt = parts[idx_alt]

                    # Filter: Remove Indels
                    if len(ref) > 1:
                        continue
                    if any(len(a) > 1 for a in alt.split(',')):
                        continue

                    # Ancestral Allele
                    pos = int(parts[idx_pos])
                    anc = chrom_seq[pos - 1] if (pos - 1) < chrom_len else "."

                    # Helper: Get unique alleles set
                    def get_alleles_set(cols):
                        s = set()
                        for i in cols:
                            gt = parts[i]
                            if gt in [".", "./.", ".|."]:
                                continue
                            for b in gt.replace('|', '/').split('/'):
                                if b != '.' and len(b) == 1:
                                    s.add(b)
                        return s

                    s1 = get_alleles_set(cols_yri)

                    # --- INGROUP PROCESSING (Split Haplotypes) ---
                    s3 = set()
                    ibs_row_values = []

                    for i in cols_ibs:
                        gt = parts[i]
                        hap1, hap2 = ".", "."

                        if gt not in [".", "./.", ".|."]:
                            alleles = gt.replace('|', '/').split('/')
                            if len(alleles) >= 2:
                                hap1, hap2 = alleles[0], alleles[1]
                            elif len(alleles) == 1:
                                hap1 = alleles[0]

                        # Append both haplotypes
                        ibs_row_values.append(hap1)
                        ibs_row_values.append(hap2)

                        # Update Set for filtering
                        if hap1 != "." and len(hap1) == 1:
                            s3.add(hap1)
                        if hap2 != "." and len(hap2) == 1:
                            s3.add(hap2)

                    if not s3:
                        continue

                    s2 = get_alleles_set(cols_nd)

                    # --- Site Filtering ---
                    # Keep if Ingroup differs from Outgroup OR Neanderthals
                    diff_yri = s3 - s1
                    diff_nd = s3 - s2

                    keep_site = False
                    if not s2:
                        if diff_yri:
                            keep_site = True
                    else:
                        if diff_yri or diff_nd:
                            keep_site = True

                    if not keep_site:
                        continue

                    # Output formatting
                    s1_str = "{" + ",".join(sorted(s1)) + "}"
                    s2_str = "{" + ",".join(sorted(s2)) + "}"
                    s3_cols_str = "\t".join(ibs_row_values)

                    # Output chromosome as normalized format (chr21)
                    output_chrom = f"{chrom_no_prefix}"

                    # Write to file
                    out_f.write(f"{output_chrom}\t{parts[idx_pos]}\t{parts[idx_ref]}\t{parts[idx_alt]}\t{anc}\t{s1_str}\t{s2_str}\t{s3_cols_str}\n")

                    row_count += 1
                    if row_count % 10000 == 0:
                        print(f"[DEBUG] Processed {row_count} rows...", file=sys.stderr)

                except Exception:
                    continue

            print(f"[INFO] Total rows written: {row_count}", file=sys.stderr)

            if row_count == 0:
                print("[WARNING] No variants passed filters! Output file is empty.", file=sys.stderr)

    except Exception as e:
        sys.stderr.write(f"[ERROR] Stream processing failed: {e}\n")
        sys.exit(1)
    finally:
        process.wait()
        for f in temp_files:
            if os.path.exists(f):
                os.remove(f)

    # Check if output file was created
    if os.path.exists(output_file):
        file_size = os.path.getsize(output_file)
        print(f"[INFO] Output file created: {output_file} ({file_size} bytes)", file=sys.stderr)
        if file_size > 0:
            with open(output_file, 'r') as f:
                line_count = sum(1 for _ in f)
            print(f"[INFO] Output file contains {line_count} lines", file=sys.stderr)
    else:
        print(f"[ERROR] Output file was not created: {output_file}", file=sys.stderr)

    elapsed = time.time() - start_time
    print(f"[INFO] Pipeline finished in {elapsed:.2f} seconds", file=sys.stderr)


if __name__ == "__main__":
    main()

