#!/usr/bin/env python3
"""
Generate sRNAs feature vectors with dynamic dimensions based on length range and feature combinations
Modes: size, size_P_5nt, sizeXstr, sizeX5nt, sizeX5ntXstr
"""

import pysam
from collections import defaultdict
import csv
import argparse
import math
import os
import multiprocessing


def process_chunk_optimized(chunk_args):
    """
    Process a chunk of references to count sRNA reads based on selected mode

    Args:
        chunk_args: Tuple containing (reference_names, bam_path, min_len, max_len, mode)

    Returns:
        Dictionary with counts for each feature combination
    """
    ref_names, bam_path, min_len, max_len, mode = chunk_args
    combined_counts = defaultdict(int)

    with pysam.AlignmentFile(bam_path, "rb") as _bam_file:
        for ref_name in ref_names:
            try:
                for read in _bam_file.fetch(ref_name):
                    # Skip unmapped and secondary alignments
                    if read.is_unmapped or read.is_secondary:
                        continue

                    length = read.query_length
                    if min_len <= length <= max_len:
                        if mode == 'size':  # Length only
                            key = (ref_name, length)
                            combined_counts[key] += 1

                        elif mode == 'size_P_5nt':  # Length + 5' nucleotides
                            key1 = (ref_name, length)
                            first_base = read.query_sequence[0].upper()
                            if first_base in ['A', 'T', 'C', 'G']:
                                key2 = (ref_name, 'base', first_base)
                                combined_counts[key1] += 1
                                combined_counts[key2] += 1

                        elif mode == 'sizeXstr':  # Length × strand
                            strand = '-' if read.is_reverse else '+'
                            key = (ref_name, length, strand)
                            combined_counts[key] += 1

                        elif mode == 'sizeX5nt':  # Length × 5' nucleotides
                            first_base = read.query_sequence[0].upper()
                            if first_base in ['A', 'T', 'C', 'G']:
                                key = (ref_name, length, first_base)
                                combined_counts[key] += 1

                        elif mode == 'sizeX5ntXstr':  # Length × 5' nucleotides × strand
                            first_base = read.query_sequence[0].upper()
                            if first_base in ['A', 'T', 'C', 'G']:
                                strand = '-' if read.is_reverse else '+'
                                key = (ref_name, length, first_base, strand)
                                combined_counts[key] += 1
            except ValueError:
                # Skip references that cannot be found in BAM file
                continue

    return dict(combined_counts)


def main():
    parser = argparse.ArgumentParser(
        description='Process BAM file to count sRNA reads with selectable feature combinations.')
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
    parser.add_argument('-m', '--mode',
                        choices=['size', 'size_P_5nt', 'sizeXstr', 'sizeX5nt', 'sizeX5ntXstr'],
                        default='sizeX5nt',
                        help='Feature combination mode: size (length only), size_P_5nt (length + 5\' nucleotides), '
                             'sizeXstr (length × strand), sizeX5nt (length × 5\' nucleotides), '
                             'sizeX5ntXstr (length × 5\' nucleotides × strand) (default: sizeX5nt)')
    args = parser.parse_args()

    # Validate length range
    if args.min_length > args.max_length:
        raise ValueError("Minimum length cannot be greater than maximum length")

    if args.min_length < 1:
        raise ValueError("Minimum length must be at least 1")

    # Open BAM file to get references
    with pysam.AlignmentFile(args.input, "rb") as bamfile:
        references = bamfile.references
        if not references:
            raise ValueError("BAM file contains no reference sequences")

    # Calculate actual dimensions
    length_count = args.max_length - args.min_length + 1

    if args.mode == 'size':
        dimensions = length_count
    elif args.mode == 'size_P_5nt':
        dimensions = length_count + 4  # +4 for A,T,C,G
    elif args.mode == 'sizeXstr':
        dimensions = length_count * 2  # ×2 for strands
    elif args.mode == 'sizeX5nt':
        dimensions = length_count * 4  # ×4 for nucleotides
    elif args.mode == 'sizeX5ntXstr':
        dimensions = length_count * 8  # ×4 for nucleotides ×2 for strands

    print(f"Found {len(references)} reference sequences")
    print(f"Counting length range: {args.min_length}-{args.max_length} nt ({length_count} lengths)")
    print(f"Selected mode: {args.mode}")
    print(f"Generated feature vector dimension: {dimensions} dimensions")

    # Split references into chunks for parallel processing
    chunk_size = max(1, math.ceil(len(references) / args.cores))
    args_list = [
        (references[i:i + chunk_size], args.input, args.min_length, args.max_length, args.mode)
        for i in range(0, len(references), chunk_size)
    ]

    print(f"Using {len(args_list)} processes for processing...")

    # Process chunks in parallel
    with multiprocessing.Pool(processes=args.cores) as pool:
        results = pool.map(process_chunk_optimized, args_list)

    # Merge results from all chunks
    merged_counts = defaultdict(int)
    for chunk_result in results:
        for key, count in chunk_result.items():
            merged_counts[key] += count

    # Define feature components
    length_range = range(args.min_length, args.max_length + 1)
    nucleotides = ['A', 'T', 'C', 'G']
    strands = ['+', '-']

    # Write results to CSV
    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Create header based on mode
        headers = ['Reference']

        if args.mode == 'size':
            # Length only
            for l in length_range:
                headers.append(f'{l}nt')

        elif args.mode == 'size_P_5nt':
            # Length + 5' nucleotides
            for l in length_range:
                headers.append(f'{l}nt')
            for nt in nucleotides:
                headers.append(f'Total_{nt}')

        elif args.mode == 'sizeXstr':
            # Length × strand
            for l in length_range:
                for strand in strands:
                    strand_name = 'plus' if strand == '+' else 'minus'
                    headers.append(f'{l}nt_{strand_name}')

        elif args.mode == 'sizeX5nt':
            # Length × 5' nucleotides
            for l in length_range:
                for nt in nucleotides:
                    headers.append(f'{l}nt_{nt}')

        elif args.mode == 'sizeX5ntXstr':
            # Length × 5' nucleotides × strand
            for l in length_range:
                for nt in nucleotides:
                    for strand in strands:
                        strand_name = 'plus' if strand == '+' else 'minus'
                        headers.append(f'{l}nt_{nt}_{strand_name}')

        writer.writerow(headers)

        # Generate feature vectors for each reference and filter all-zero rows
        print("Generating feature vectors and filtering all-zero rows...")
        non_zero_rows = []

        for ref_name in references:
            row = [ref_name]

            if args.mode == 'size':
                # Length only
                for l in length_range:
                    count = merged_counts.get((ref_name, l), 0)
                    row.append(count)

            elif args.mode == 'size_P_5nt':
                # Length + 5' nucleotides
                for l in length_range:
                    count = merged_counts.get((ref_name, l), 0)
                    row.append(count)
                for nt in nucleotides:
                    count = merged_counts.get((ref_name, 'base', nt), 0)
                    row.append(count)

            elif args.mode == 'sizeXstr':
                # Length × strand
                for l in length_range:
                    for strand in strands:
                        count = merged_counts.get((ref_name, l, strand), 0)
                        row.append(count)

            elif args.mode == 'sizeX5nt':
                # Length × 5' nucleotides
                for l in length_range:
                    for nt in nucleotides:
                        count = merged_counts.get((ref_name, l, nt), 0)
                        row.append(count)

            elif args.mode == 'sizeX5ntXstr':
                # Length × 5' nucleotides × strand
                for l in length_range:
                    for nt in nucleotides:
                        for strand in strands:
                            count = merged_counts.get((ref_name, l, nt, strand), 0)
                            row.append(count)
            if sum(row[1:]) > 0:
                non_zero_rows.append(row)
                
        print(f"Writing {len(non_zero_rows)} non-zero rows to CSV (filtered out {len(references) - len(non_zero_rows)} all-zero rows)")
        writer.writerow(headers)
        for row in non_zero_rows:
            writer.writerow(row)

    print(f"Processing completed! Results saved to {args.output}")
    print(f"Final feature vector: {dimensions} dimensions")
    print(f"Feature combination: {args.mode}")
    print(f"Length range: {args.min_length}-{args.max_length} nt")


if __name__ == '__main__':
    main()
