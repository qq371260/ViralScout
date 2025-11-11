#!/usr/bin/env python3
"""
Generate sRNAs feature vectors with 13, 17, 26, 52, 104 dimensions
"""
import pysam
from collections import defaultdict
import csv
import argparse
import math
import os
import multiprocessing


def process_chunk_optimized(chunk_args):
    """Processing function to count all combinations based on selected dimensions"""
    ref_names, bam_path, min_len, max_len, dimension = chunk_args
    combined_counts = defaultdict(int)

    with pysam.AlignmentFile(bam_path, "rb") as _bam_file:
        for ref_name in ref_names:
            try:
                for read in _bam_file.fetch(ref_name):
                    if read.is_unmapped or read.is_secondary:
                        continue
                    length = read.query_length
                    if min_len <= length <= max_len:
                        if dimension == 13:  # size
                            key = (ref_name, length)
                            combined_counts[key] += 1
                        elif dimension == 17:  # size + 5'nt
                            key1 = (ref_name, length)
                            key2 = (ref_name, 'total')
                            first_base = read.query_sequence[0].upper()
                            if first_base in ['A', 'T', 'C', 'G']:
                                key3 = (ref_name, 'base', first_base)
                                combined_counts[key1] += 1
                                combined_counts[key2] += 1
                                combined_counts[key3] += 1
                        elif dimension == 26:  # size × 2 strands
                            strand = '-' if read.is_reverse else '+'
                            key = (ref_name, length, strand)
                            combined_counts[key] += 1
                        elif dimension == 52:  # size × 5'nt
                            first_base = read.query_sequence[0].upper()
                            if first_base in ['A', 'T', 'C', 'G']:
                                key = (ref_name, length, first_base)
                                combined_counts[key] += 1
                        elif dimension == 104:  # size × 5'nt × 2 strands
                            first_base = read.query_sequence[0].upper()
                            if first_base in ['A', 'T', 'C', 'G']:
                                strand = '-' if read.is_reverse else '+'
                                key = (ref_name, length, first_base, strand)
                                combined_counts[key] += 1
            except ValueError:
                continue

    return dict(combined_counts)


def main():
    parser = argparse.ArgumentParser(
        description='Process BAM file to count sRNA reads with selectable dimensions.')
    parser.add_argument('-i', '--input', required=True,
                        help='Input BAM file path')
    parser.add_argument('-o', '--output', default='sRNA_features.csv',
                        help='Output CSV file path')
    parser.add_argument('-c', '--cores', type=int, default=os.cpu_count(),
                        help=f"Number of CPU cores to use (default: {os.cpu_count()})")
    parser.add_argument('--min_length', type=int, default=18,
                        help='Minimum read length (default: 18)')
    parser.add_argument('--max_length', type=int, default=30,
                        help='Maximum read length (default: 30)')
    parser.add_argument('-d', '--dimension', type=int, choices=[13, 17, 26, 52, 104], default=52,
                        help='Feature vector dimension: 13 (length only), 17 (length + total 5\' nt), '
                             '26 (length × strand), 52 (length × 5\' nt), 104 (length × 5\' nt × strand) (default: 52)')
    args = parser.parse_args()

    with pysam.AlignmentFile(args.input, "rb") as bamfile:
        references = bamfile.references
        if not references:
            raise ValueError("BAM file contains no reference sequences")

    print(f"Found {len(references)} reference sequences")
    print(f"Counting length range: {args.min_length}-{args.max_length} nt")
    print(f"Generating {args.dimension}-dimensional feature vector")

    # Chunk strategy
    chunk_size = max(1, math.ceil(len(references) / args.cores))
    args_list = [
        (references[i:i + chunk_size], args.input, args.min_length, args.max_length, args.dimension)
        for i in range(0, len(references), chunk_size)
    ]

    print(f"Using {len(args_list)} processes for processing...")

    # Multiprocessing
    with multiprocessing.Pool(processes=args.cores) as pool:
        results = pool.map(process_chunk_optimized, args_list)

    # Results merged
    merged_counts = defaultdict(int)
    for chunk_result in results:
        for key, count in chunk_result.items():
            merged_counts[key] += count

    # Feature vectors
    length_range = range(args.min_length, args.max_length + 1)
    nucleotides = ['A', 'T', 'C', 'G']
    strands = ['+', '-']

    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        headers = ['Reference']

        if args.dimension == 13:
            # 13D: size
            for l in length_range:
                headers.append(f'{l}nt')

        elif args.dimension == 17:
            # 17D: size + 5'nt
            for l in length_range:
                headers.append(f'{l}nt')
            for nt in nucleotides:
                headers.append(f'Total_{nt}')

        elif args.dimension == 26:
            # 26D: size × 2 strands
            for l in length_range:
                for strand in strands:
                    strand_name = 'plus' if strand == '+' else 'minus'
                    headers.append(f'{l}nt_{strand_name}')

        elif args.dimension == 52:
            # 52D: size × 5'nt
            for l in length_range:
                for nt in nucleotides:
                    headers.append(f'{l}nt_{nt}')

        elif args.dimension == 104:
            # 104D: size × 5'nt × 2 strands
            for l in length_range:
                for nt in nucleotides:
                    for strand in strands:
                        strand_name = 'plus' if strand == '+' else 'minus'
                        headers.append(f'{l}nt_{nt}_{strand_name}')

        writer.writerow(headers)

        # Feature vectors
        for ref_name in references:
            row = [ref_name]

            if args.dimension == 13:
                # size
                for l in length_range:
                    count = merged_counts.get((ref_name, l), 0)
                    row.append(count)

            elif args.dimension == 17:
                # size + 5'nt
                for l in length_range:
                    count = merged_counts.get((ref_name, l), 0)
                    row.append(count)
                for nt in nucleotides:
                    count = merged_counts.get((ref_name, 'base', nt), 0)
                    row.append(count)

            elif args.dimension == 26:
                # size × 2 strands
                for l in length_range:
                    for strand in strands:
                        count = merged_counts.get((ref_name, l, strand), 0)
                        row.append(count)

            elif args.dimension == 52:
                # size × 5'nt
                for l in length_range:
                    for nt in nucleotides:
                        count = merged_counts.get((ref_name, l, nt), 0)
                        row.append(count)

            elif args.dimension == 104:
                # size × 5'nt × 2 strands
                for l in length_range:
                    for nt in nucleotides:
                        for strand in strands:
                            count = merged_counts.get((ref_name, l, nt, strand), 0)
                            row.append(count)

            writer.writerow(row)

    print(f"Processing completed! Results saved to {args.output}")
    print(f"Generated feature vector dimension: {args.dimension} dimensions")

    if args.dimension == 13:
        print(f"Feature order: Length only ({args.min_length}-{args.max_length}nt)")
    elif args.dimension == 17:
        print(f"Feature order: Length ({args.min_length}-{args.max_length}nt) + Total 5' nucleotides (A,T,C,G)")
    elif args.dimension == 26:
        print(f"Feature order: All combinations of length ({args.min_length}-{args.max_length}nt) and strands (plus, minus)")
    elif args.dimension == 52:
        print(
            f"Feature order: All combinations of length ({args.min_length}-{args.max_length}nt) and 5' nucleotides (A,T,C,G)")
    elif args.dimension == 104:
        print(
            f"Feature order: All combinations of length ({args.min_length}-{args.max_length}nt), 5' nucleotides (A,T,C,G), and strands (+, -)")


if __name__ == '__main__':
    main()
