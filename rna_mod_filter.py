#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pysam
from pybedtools import BedTool

# Map modification type to the expected value in column 4 of the BED file
mod_dic = {
    "m6a": "a",
    "m5c": "m",
    "psi": "17802"
}

def filter_bed(bed_file, fasta_file, mod_type, output_file, base):
    mod_code = mod_dic.get(mod_type)
    if mod_code is None:
        raise ValueError(f"Unsupported mod_type: {mod_type}")

    # Load the reference genome
    fasta = pysam.FastaFile(fasta_file)
    filtered_lines = []

    for interval in BedTool(bed_file):
        chrom = interval.chrom
        pos = int(interval.start)
        strand = interval.strand

        # Fetch the reference base from the genome
        ref_base = fasta.fetch(chrom, pos, pos + 1).strip().upper()

        # If it's on the negative strand, get the reverse complement
        if strand == "-":
            complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
            ref_base = complement.get(ref_base, "N")

        # Parse fields from the BED line
        fields = interval.fields
        mod = str(fields[3]).strip()
        valid_cov = int(fields[4])
        mod_rate = float(fields[10])

        # Apply all filters:
        # 1. correct modification type
        # 2. reference base matches target (default: T)
        # 3. coverage >= 20
        # 4. modification rate >= 5%
        if ref_base == base and mod == mod_code and valid_cov >= 20 and mod_rate >= 5:
            filtered_lines.append(str(interval))

    # Write filtered lines to output
    with open(output_file, "w") as f:
        f.writelines(line + "\n" for line in filtered_lines)

def main():
    parser = argparse.ArgumentParser(description="Filter modkit BED entries by reference base and thresholds")
    parser.add_argument("--bed", required=True, help="Input BED file from modkit")
    parser.add_argument("--fasta", required=True, help="Reference FASTA file (must be indexed)")
    parser.add_argument("--mod", required=True, choices=mod_dic.keys(), help="Modification type: m6a, m5c, psi")
    parser.add_argument("--base", default="T", help="Reference base to match (default: T)")
    parser.add_argument("--out", required=True, help="Output BED file with filtered entries")
    args = parser.parse_args()

    filter_bed(args.bed, args.fasta, args.mod, args.out, args.base.upper())

if __name__ == "__main__":
    main()