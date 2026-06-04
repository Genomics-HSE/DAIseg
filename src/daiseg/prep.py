"""
Prepare merged SNP table for DAIseg.

1. 1000G is filtered by the 1000G strict mask, as before.
2. Each Neanderthal VCF is filtered by its own mask intersected with the 1000G strict mask.
3. Original 1000G exact positions and SV intervals are loaded in one pass from the original VCF.
4. Sites are kept when IBS has an allele absent from YRI or absent from Neanderthals.
5. If IBS/YRI are missing after merge, they are treated as REF/REF only when the exact site was absent from the original 1000G VCF and the site does not overlap an original 1000G SV interval.
"""

import os
import sys
import subprocess
import time
import signal
from bisect import bisect_right
from collections import defaultdict
import pysam
import daiseg.utils as utils


child_procs = []
temp_files = []


def run_step(cmd, desc):
    try:
        subprocess.check_call(cmd, shell=True, preexec_fn=os.setsid)
    except subprocess.CalledProcessError:
        sys.exit(f"[ERROR] Failed at step: {desc}")


def terminate_children():
    for proc in child_procs:
        try:
            if proc.poll() is None:
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        except Exception:
            pass

    time.sleep(1)

    for proc in child_procs:
        try:
            if proc.poll() is None:
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except Exception:
            pass


def cleanup_temp_files():
    for f in temp_files:
        if os.path.exists(f):
            try:
                os.remove(f)
            except Exception:
                pass


def handle_stop_signal(signum, frame):
    print("[INFO] Interrupted. Terminating child processes...", file=sys.stderr)
    terminate_children()
    cleanup_temp_files()
    sys.exit("[INFO] Stopped by user.")


def run(json_path, threads=16):

    signal.signal(signal.SIGINT, handle_stop_signal)
    signal.signal(signal.SIGTERM, handle_stop_signal)

    start_time = time.time()

    cfg = utils.load_config(json_path)
    chrom_raw = str(cfg["CHROM"])
    files = cfg["files"]
    prefix = cfg["prefix"]

    threads_1kg = max(1, threads // 2)
    threads_nd = max(1, threads // 4)
    threads_merge = max(1, threads)

    output_filename = cfg.get(
        "data", f"prep.chr{chrom_raw.lstrip('chr').lstrip('CHR')}.tsv"
    )
    output_file = os.path.join(prefix, output_filename)

    print(f" [INFO] Output will be written to: {output_file}", file=sys.stderr)

    chrom_no_prefix = chrom_raw.lstrip("chr").lstrip("CHR")

    vcf_1kg = prefix + "/" + utils.expand_path(files["1000GP_files"]["vcf"])
    vcf_1kg_initial = utils.expand_path(files["1000GP_files"].get("vcf_initial", ""))
    bed_strict = utils.expand_path(files["1000GP_files"]["bed"])
    anc_path = utils.expand_path(files["ancestral"]["fasta"])

    if not os.path.exists(anc_path):
        sys.exit(f"[ERROR] Ancestral FASTA not found: {anc_path}")

    original_1kg_positions = set()
    original_1kg_sv_intervals = []

    if vcf_1kg_initial and os.path.exists(vcf_1kg_initial):
        print(" [INFO] Loading original 1000G positions and SV intervals...", file=sys.stderr)
        cmd_original = (
            f"bcftools query -r {chrom_raw} "
            f"-f '%CHROM\\t%POS\\t%INFO/END\\t%ALT\\n' "
            f"{vcf_1kg_initial}"
        )
        proc_original = subprocess.Popen(
            cmd_original,
            shell=True,
            stdout=subprocess.PIPE,
            text=True,
            preexec_fn=os.setsid,
        )
        child_procs.append(proc_original)

        for line in proc_original.stdout:
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < 4:
                continue

            chrom, pos, end, alt = fields[:4]
            original_1kg_positions.add((chrom, pos))

            if "<" in alt and ">" in alt:
                start_pos = int(pos)
                end_pos = int(end) if end not in ["", "."] else start_pos
                original_1kg_sv_intervals.append((chrom, start_pos, end_pos))

        if proc_original.wait() != 0:
            cleanup_temp_files()
            sys.exit("[ERROR] Failed to load original 1000G positions/SV intervals.")

        print(
            f" [INFO] Original 1000G positions loaded: {len(original_1kg_positions)}",
            file=sys.stderr,
        )
        print(
            f" [INFO] Original 1000G SV intervals loaded: {len(original_1kg_sv_intervals)}",
            file=sys.stderr,
        )
    else:
        print(
            " [WARNING] Original 1000G VCF not found; REF/REF restoration will be disabled.",
            file=sys.stderr,
        )

    def build_sv_index(intervals):
        by_chrom = defaultdict(list)
        for chrom, start, end in intervals:
            by_chrom[chrom].append((start, end))

        index = {}
        for chrom, chrom_intervals in by_chrom.items():
            chrom_intervals.sort()
            merged = []
            for start, end in chrom_intervals:
                if not merged or start > merged[-1][1] + 1:
                    merged.append([start, end])
                elif end > merged[-1][1]:
                    merged[-1][1] = end

            starts = [x[0] for x in merged]
            ends = [x[1] for x in merged]
            index[chrom] = (starts, ends)
        return index

    original_1kg_sv_index = build_sv_index(original_1kg_sv_intervals)
    print(
        f" [INFO] Original 1000G merged SV intervals indexed: {sum(len(v[0]) for v in original_1kg_sv_index.values())}",
        file=sys.stderr,
    )

    def in_original_1kg_sv(chrom, pos):
        starts, ends = original_1kg_sv_index.get(chrom, ([], []))
        if not starts:
            return False
        pos = int(pos)
        i = bisect_right(starts, pos) - 1
        return i >= 0 and pos <= ends[i]

    print(f" [INFO] Creating temporary files...", file=sys.stderr)

    tmp_1kg = f"{prefix}/temp_1kg_{chrom_no_prefix}.bcf"
    cmd_1kg = (
        f"bcftools view --threads {threads_1kg} -r {chrom_raw} -R {bed_strict} "
        f"-O b -o {tmp_1kg} {vcf_1kg}"
    )
    run_step(cmd_1kg, "Filtering 1000 Genomes")
    run_step(f"bcftools index -f {tmp_1kg}", "Indexing 1000G")

    temp_files.append(tmp_1kg)
    temp_files.append(f"{tmp_1kg}.csi")

    print(f" [INFO] Working with Neanderthals (Parallel)...", file=sys.stderr)
    neand_files = files.get("neand_files", {})
    nd_inputs = []
    running_procs = []

    for i, (name, paths) in enumerate(neand_files.items()):
        n_vcf = utils.expand_path(paths["vcf"])
        n_bed = utils.expand_path(paths["bed"])

        tmp_nd = f"{prefix}/temp_nd_{i}_{chrom_no_prefix}.bcf"

        bed_arg = ""
        if n_bed:
            n_bed_used = ""
            if os.path.exists(n_bed):
                n_bed_used = n_bed
            else:
                bed_dir = os.path.dirname(n_bed)
                bed_base = os.path.basename(n_bed)

                if bed_base.startswith("chr"):
                    alt_bed = os.path.join(bed_dir, bed_base[3:])
                    if os.path.exists(alt_bed):
                        n_bed_used = alt_bed
                        print(f" [INFO] Using alternative BED: {alt_bed}", file=sys.stderr)

            if n_bed_used:
                tmp_bed = f"{prefix}/temp_nd_{i}_{chrom_no_prefix}.bed"
                run_step(
                    f"zcat -f {n_bed_used} | bedtools intersect -a - -b {bed_strict} > {tmp_bed}",
                    f"Intersecting Neanderthal BED with 1000G mask for {name}",
                )
                bed_arg = f"-T {tmp_bed}"
                temp_files.append(tmp_bed)

        cmd_str = (
            f"bcftools view --threads {threads_nd} -r {chrom_raw} {bed_arg} -O b -o {tmp_nd} {n_vcf} && "
            f"bcftools index -f {tmp_nd}"
        )

        print(f" [INFO] Starting job for {name}...", file=sys.stderr)

        proc = subprocess.Popen(cmd_str, shell=True, preexec_fn=os.setsid)
        running_procs.append(proc)
        child_procs.append(proc)

        temp_files.append(tmp_nd)
        temp_files.append(f"{tmp_nd}.csi")
        nd_inputs.append(tmp_nd)

    for proc in running_procs:
        if proc.wait() != 0:
            cleanup_temp_files()
            sys.exit("[ERROR] One of the Neanderthal processing jobs failed.")

    files_str = " ".join([tmp_1kg] + nd_inputs)

    pipeline_cmd = (
        f"bcftools merge --threads {threads_merge} --force-samples --merge all -O u {files_str} "
        f"| bcftools query -H -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%TGT]\\n'"
    )

    try:
        af = pysam.FastaFile(anc_path)

        available_chroms = af.references
        print(
            f"[INFO] Available chromosomes in ancestral FASTA: {available_chroms}",
            file=sys.stderr,
        )
        print(
            f"[INFO] Looking for chromosome related to: '{chrom_raw}'", file=sys.stderr
        )

        chrom_seq = None
        chrom_name_used = None

        if len(available_chroms) == 1:
            chrom_name_used = available_chroms[0]
            chrom_seq = af.fetch(chrom_name_used)
            print(
                f"[INFO] Using only available chromosome: '{chrom_name_used}'",
                file=sys.stderr,
            )
        else:
            search_terms = [
                f":{chrom_no_prefix}:",
                f"chromosome.*{chrom_no_prefix}",
                f"chr{chrom_no_prefix}",
                chrom_no_prefix,
                f"CHR{chrom_no_prefix}",
            ]

            for chrom_name in available_chroms:
                for term in search_terms:
                    if term in chrom_name:
                        chrom_name_used = chrom_name
                        chrom_seq = af.fetch(chrom_name_used)
                        print(
                            f"[INFO] Found chromosome using pattern '{term}': '{chrom_name_used}'",
                            file=sys.stderr,
                        )
                        break
                if chrom_seq:
                    break

        if chrom_seq is None:
            for chrom_name in available_chroms:
                if (
                    chrom_name == chrom_raw
                    or chrom_name == f"chr{chrom_no_prefix}"
                    or chrom_name == chrom_no_prefix
                ):
                    chrom_name_used = chrom_name
                    chrom_seq = af.fetch(chrom_name_used)
                    print(
                        f"[INFO] Found exact match: '{chrom_name_used}'",
                        file=sys.stderr,
                    )
                    break

            if chrom_seq is None:
                for chrom_name in available_chroms:
                    if chrom_no_prefix in chrom_name or chrom_raw in chrom_name:
                        chrom_name_used = chrom_name
                        chrom_seq = af.fetch(chrom_name_used)
                        print(
                            f"[INFO] Found partial match: '{chrom_name_used}'",
                            file=sys.stderr,
                        )
                        break

        if chrom_seq is None:
            error_msg = f"Chromosome '{chrom_raw}' not found in Ancestral FASTA.\n"
            error_msg += f"Available chromosomes: {available_chroms}\n"
            error_msg += f"Searching for patterns containing: '{chrom_no_prefix}'\n"
            error_msg += "Please check your FASTA file format and chromosome naming."
            raise ValueError(error_msg)

        chrom_len = len(chrom_seq)
        print(
            f"[INFO] Chromosome '{chrom_name_used}' length: {chrom_len} bp",
            file=sys.stderr,
        )
        af.close()
    except Exception as e:
        cleanup_temp_files()
        sys.exit(f"[ERROR] Ancestral fasta issue: {e}")

    process = subprocess.Popen(
        pipeline_cmd,
        shell=True,
        stdout=subprocess.PIPE,
        text=True,
        preexec_fn=os.setsid,
    )
    child_procs.append(process)

    try:
        with open(output_file, "w") as out_f:
            header_line = process.stdout.readline().strip()
            if not header_line.startswith("#"):
                raise RuntimeError("Pipeline returned no header.")

            raw_headers = header_line.lstrip("#").split("\t")
            headers = []
            for h in raw_headers:
                if "]" in h:
                    h = h.split("]")[-1]
                if ":" in h:
                    h = h.split(":")[0]
                headers.append(h.strip())

            idx_chrom, idx_pos, idx_ref, idx_alt, cols_yri, cols_ibs, cols_nd = (
                utils.map_columns(headers, cfg["samples"])
            )

            ibs_split_headers = []
            for i in cols_ibs:
                name = headers[i]
                ibs_split_headers.append(f"{name}_1")
                ibs_split_headers.append(f"{name}_2")

            ibs_header_str = "\t".join(ibs_split_headers)

            out_f.write(
                f"#CHROM\tPOS\tREF\tALT\tAncestral\tOutgroup\tNeand\t{ibs_header_str}\n"
            )

            row_count = 0

            for line in process.stdout:
                parts = line.strip().split("\t")
                if len(parts) != len(headers):
                    continue

                try:
                    ref = parts[idx_ref]
                    alt = parts[idx_alt]

                    if len(ref) > 1:
                        continue
                    if any(len(a) > 1 for a in alt.split(",")):
                        continue

                    pos = int(parts[idx_pos])
                    anc = chrom_seq[pos - 1] if (pos - 1) < chrom_len else "."

                    def get_alleles_set(cols):
                        s = set()
                        for i in cols:
                            gt = parts[i]
                            if gt in [".", "./.", ".|."]:
                                continue
                            for b in gt.replace("|", "/").split("/"):
                                if b != "." and len(b) == 1:
                                    s.add(b)
                        return s

                    s1 = get_alleles_set(cols_yri)

                    s3 = set()
                    ibs_row_values = []

                    for i in cols_ibs:
                        gt = parts[i]
                        hap1, hap2 = ".", "."

                        if gt not in [".", "./.", ".|."]:
                            alleles = gt.replace("|", "/").split("/")
                            if len(alleles) >= 2:
                                hap1, hap2 = alleles[0], alleles[1]
                            elif len(alleles) == 1:
                                hap1 = alleles[0]

                        ibs_row_values.append(hap1)
                        ibs_row_values.append(hap2)

                        if hap1 != "." and len(hap1) == 1:
                            s3.add(hap1)
                        if hap2 != "." and len(hap2) == 1:
                            s3.add(hap2)

                    s2 = get_alleles_set(cols_nd)

                    pos_key = (parts[idx_chrom], parts[idx_pos])
                    can_restore_ref = (
                        bool(original_1kg_positions)
                        and pos_key not in original_1kg_positions
                        and not in_original_1kg_sv(parts[idx_chrom], parts[idx_pos])
                    )

                    if not s3 and can_restore_ref:
                        s3 = {ref}
                        if not s1:
                            s1 = {ref}
                        ibs_row_values = []
                        for _ in cols_ibs:
                            ibs_row_values.append(ref)
                            ibs_row_values.append(ref)

                    if not s3:
                        continue

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

                    s1_str = "{" + ",".join(sorted(s1)) + "}"
                    s2_str = "{" + ",".join(sorted(s2)) + "}"
                    s3_cols_str = "\t".join(ibs_row_values)

                    output_chrom = f"{chrom_no_prefix}"

                    out_f.write(
                        f"{output_chrom}\t{parts[idx_pos]}\t{parts[idx_ref]}\t{parts[idx_alt]}\t{anc}\t{s1_str}\t{s2_str}\t{s3_cols_str}\n"
                    )

                    row_count += 1
                    if row_count % 10000 == 0:
                        print(f"[INFO] Processed {row_count} rows...", file=sys.stderr)

                except Exception:
                    continue

            print(f"[INFO] Total rows written: {row_count}", file=sys.stderr)

            if row_count == 0:
                print(
                    "[WARNING] No variants passed filters! Output file is empty.",
                    file=sys.stderr,
                )

    except Exception as e:
        sys.stderr.write(f"[ERROR] Stream processing failed: {e}\n")
        terminate_children()
        sys.exit(1)
    finally:
        if process.poll() is None:
            try:
                process.wait()
            except Exception:
                pass
        cleanup_temp_files()

    if os.path.exists(output_file):
        file_size = os.path.getsize(output_file)
        print(
            f"[INFO] Output file created: {output_file} ({file_size} bytes)",
            file=sys.stderr,
        )
        if file_size > 0:
            with open(output_file, "r") as f:
                line_count = sum(1 for _ in f)
            print(f"[INFO] Output file contains {line_count} lines", file=sys.stderr)
    else:
        print(f"[ERROR] Output file was not created: {output_file}", file=sys.stderr)

    elapsed = time.time() - start_time
    print(f"[INFO] Pipeline finished in {elapsed:.2f} seconds", file=sys.stderr)
