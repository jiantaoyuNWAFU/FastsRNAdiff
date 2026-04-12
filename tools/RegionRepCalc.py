#!/usr/bin/env python3
import sys
import time
import pysam
import os

def main():
    if len(sys.argv) != 4:
        print("Usage: python RegionRepCalc.py <bam file> <sRNA cluster file> <output file>", file=sys.stderr)
        print("Function: Count Reads, Rep-total, UniqueReads, DicerCall without any tags like XX.")
        sys.exit(1)

    bam_file = sys.argv[1]
    region_file = sys.argv[2]
    output_file = sys.argv[3]
    total_output_file = output_file + ".total"

    print("Start calculating sRNA cluster expression...")
    start_time = time.time()

    global_qname_set = set()

    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except FileNotFoundError:
        print(f"Error: {bam_file} does not exist!", file=sys.stderr)
        sys.exit(1)

    bai_file = bam_file + ".bai"
    if not os.path.exists(bai_file):
        print(f"Warning: Creating index for {bam_file}...", file=sys.stderr)
        try:
            pysam.index(bam_file)
        except Exception as e:
            print(f"Error: Failed to create index: {e}", file=sys.stderr)
            sys.exit(1)

    with open(output_file, "w") as f:
        f.write("#Locus\tReads\tRep-total\tUniqueReads\tDicerCall\n")

    total_count = 0
    success_count = 0

    print("Processing sRNA cluster file...")
    with open(region_file, "r") as f:
        header = f.readline().strip().split()
        locus_idx = header.index("#Locus")
        dicer_idx = header.index("DicerCall") if 'DicerCall' in header else -1

        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            total_count += 1
            if total_count % 1000 == 0:
                print(f"Finished {total_count} regions...")

            parts = line.split("\t")
            locus = parts[locus_idx]

            if dicer_idx != -1:
                dicer = parts[dicer_idx]
            else:
                dicer = 24

            try:
                chrom, pos = locus.split(":")
                start_str, end_str = pos.split("-")
                start = int(start_str)
                end = int(end_str)
            except:
                with open(output_file, "a") as out_f:
                    out_f.write(f"{locus}\t0\t0\t0\tNotValid\n")
                continue

            if start > end:
                with open(output_file, "a") as out_f:
                    out_f.write(f"{locus}\t0\t0\t0\tNotValid\n")
                continue

            read_counts = {}
            total_reads = 0
            unique_reads = 0

            try:
                for read in bam.fetch(chrom, start - 1, end):
                    qname = read.qname
                    total_reads += 1

                    if qname in read_counts:
                        read_counts[qname] += 1
                    else:
                        read_counts[qname] = 1

                    global_qname_set.add(qname)

            except Exception as e:
                print(f"Warning: Error processing {locus}: {e}", file=sys.stderr)
                with open(output_file, "a") as out_f:
                    out_f.write(f"{locus}\t0\t0.000\t0\t{dicer}\n")
                continue

            rep_total = 0.0
            for cnt in read_counts.values():
                rep_total += 1.0 / cnt

            unique_reads = sum(1 for cnt in read_counts.values() if cnt == 1)

            with open(output_file, "a") as out_f:
                out_f.write(f"{locus}\t{total_reads}\t{rep_total:.3f}\t{unique_reads}\t{dicer}\n")

            success_count += 1

    rep_total_global = len(global_qname_set)
    with open(total_output_file, "w") as f:
        f.write(f"Repeat-normalized Total Mapped Reads: {rep_total_global}\n")

    print("\n=============================================")
    print("Processing completed!")
    print(f"Total time: {int(time.time() - start_time)}s")
    print(f"Total clusters: {total_count}")
    print(f"Success: {success_count}")
    print(f"Total Mapped Reads: {rep_total_global}")
    print(f"Result saved to: {output_file}")
    print(f"Total count saved to: {total_output_file}")
    print("=============================================")

if __name__ == "__main__":
    main()
