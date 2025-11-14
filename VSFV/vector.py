#!/usr/bin/env python3
"""
Generate sRNAs feature vectors for specific contigs with dynamic dimensions based on length range and feature combinations
Modes: size, size_P_5nt, sizeXstr, sizeX5nt, sizeX5ntXstr
For sizeXstr and sizeX5ntXstr modes, generate both plus and minus strand versions for all vectors
"""

import pysam
from collections import defaultdict
import csv
import argparse
import math
import os
import multiprocessing
import re
import pandas as pd


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Process BAM file to count sRNA reads with selectable feature combinations for specific contigs.')
    parser.add_argument('-i', '--input', required=True,
                        help='Input BAM file path')
    parser.add_argument('-o', '--output', default='sRNA_features.csv',
                        help='Output CSV file path')
    parser.add_argument('-f', '--filtered_contigs', type=str, default=None,
                        help='File containing filtered contig names (txt or csv)')
    parser.add_argument('-b', '--benchmarks', type=str, default=None,
                        help='Benchmark contig txt or csv file path. If absent, use default vsiRNA simulant reference')
    parser.add_argument('-c', '--cores', type=int, default=os.cpu_count(),
                        help='Number of CPU cores to use')
    parser.add_argument('--min_length', type=int, default=18,
                        help='Minimum sRNA size')
    parser.add_argument('--max_length', type=int, default=30,
                        help='Maximum sRNA size')
    parser.add_argument('-m', '--mode',
                        choices=['size', 'size_P_5nt', 'sizeXstr', 'sizeX5nt', 'sizeX5ntXstr'], default='sizeX5nt',
                        help='Five feature vector modes')
    return parser.parse_args()


def extract_first_word(name):
    first_word = re.split(r'[\s\t,;]', name.strip())[0]
    return first_word


def load_contig_list(contig_file):
    """Load contig names from file"""
    print(f"Loading contig list from: {contig_file}")

    if contig_file.endswith('.csv'):
        df = pd.read_csv(contig_file)
        contig_names = df.iloc[:, 0].tolist()
    else:
        # Assume txt file with one contig name per line
        with open(contig_file, 'r') as f:
            contig_names = [line.strip() for line in f if line.strip()]

    # Extract first word from each contig name
    processed_names = [extract_first_word(name) for name in contig_names]
    print(f"Loaded {len(processed_names)} contigs")
    return processed_names


def load_benchmarks(benchmarks_file, mode='sizeX5nt'):
    """Load benchmark contigs from file or use default vsiRNA simulant"""

    if benchmarks_file is None:
        # Use default vsiRNA simulant file
        default_filename = f"vsiRNA_simulant_{mode}.csv"
        possible_locations = [
            default_filename,
            f"../{default_filename}",
            f"./{default_filename}"
        ]

        for ref_file in possible_locations:
            if os.path.isfile(ref_file):
                print(f"Using default vsiRNA simulant file: {ref_file}")
                ref_df = pd.read_csv(ref_file)
                # Return the DataFrame and a flag indicating it's a default vector file
                return None, ref_df

        print(f"Error: No benchmarks provided and default vsiRNA simulant files not found")
        print(f"Expected files: {default_filename}")
        return None, None

    else:
        # User provided benchmark file with contig names
        print(f"Loading benchmark contigs from: {benchmarks_file}")

        if benchmarks_file.endswith('.csv'):
            ref_df = pd.read_csv(benchmarks_file)
            benchmark_names = ref_df.iloc[:, 0].tolist()
        else:
            # Assume txt file with one contig name per line
            with open(benchmarks_file, 'r') as f:
                benchmark_names = [line.strip() for line in f if line.strip()]

        # Extract first word
        benchmark_names = [extract_first_word(name) for name in benchmark_names]
        print(f"Loaded {len(benchmark_names)} benchmark contigs")
        return benchmark_names, None


def process_chunk_optimized(chunk_args):
    """Process a chunk of references to count sRNA reads based on selected mode"""

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


def swap_strand_counts(row, mode, length_range, nucleotides):
    """Create a minus strand version by swapping plus and minus strand counts"""

    if mode == 'sizeXstr':
        # For sizeXstr: swap plus and minus for each length
        swapped_row = [row[0] + '_minus']  # Reference name with _minus suffix
        for i in range(1, len(row), 2):  # Process in pairs: plus, minus
            if i + 1 < len(row):
                swapped_row.append(row[i + 1])  # minus becomes first
                swapped_row.append(row[i])  # plus becomes second
        return swapped_row

    elif mode == 'sizeX5ntXstr':
        # For sizeX5ntXstr: swap plus and minus for each length and nucleotide
        swapped_row = [row[0] + '_minus']  # Reference name with _minus suffix
        # Each length has 8 columns: A_plus, A_minus, T_plus, T_minus, C_plus, C_minus, G_plus, G_minus
        for i in range(1, len(row), 8):
            for nt_idx in range(4):  # For each nucleotide (A, T, C, G)
                plus_idx = i + (nt_idx * 2)
                minus_idx = plus_idx + 1
                if minus_idx < len(row):
                    swapped_row.append(row[minus_idx])  # minus becomes first
                    swapped_row.append(row[plus_idx])  # plus becomes second
        return swapped_row

    else:
        return None


def process_external_vectors(external_vectors_df, mode, length_range, nucleotides):
    """Process external vectors (vsiRNA simulant) to generate plus-minus strand versions for strand-specific modes"""
    print(f"Processing external vectors for {mode} mode...")

    external_rows = []

    # Iterate through each row in the external vectors DataFrame
    for idx, row in external_vectors_df.iterrows():
        # Convert the row to a list
        row_list = [row.iloc[0]] + row.iloc[1:].tolist()

        # Add the original row
        external_rows.append(row_list)

        # For strand-specific modes, also create minus strand version
        if mode in ['sizeXstr', 'sizeX5ntXstr']:
            minus_row = swap_strand_counts(row_list, mode, length_range, nucleotides)
            if minus_row:
                external_rows.append(minus_row)

    return external_rows


def main():
    # Parse arguments
    args = parse_arguments()

    # Validate length range
    if args.min_length > args.max_length:
        raise ValueError("Minimum length cannot be greater than maximum length")

    if args.min_length < 1:
        raise ValueError("Minimum length must be at least 1")

    # Load contig lists
    bam_contigs = []  # Contigs that exist in BAM file
    external_vectors_df = None  # External vectors (like vsiRNA simulants)

    # Load filtered contigs if provided
    if args.filtered_contigs:
        filtered_contigs = load_contig_list(args.filtered_contigs)
        bam_contigs.extend(filtered_contigs)

    # Load benchmarks
    benchmark_names, benchmark_vectors_df = load_benchmarks(args.benchmarks, args.mode)

    if benchmark_vectors_df is not None:
        # We have default vsiRNA vectors (external vectors)
        print("Using default vsiRNA simulant vectors as external references")
        external_vectors_df = benchmark_vectors_df
    elif benchmark_names:
        # We have benchmark contig names that should exist in BAM file
        bam_contigs.extend(benchmark_names)

    # Remove duplicates
    bam_contigs = list(set(bam_contigs))

    if not bam_contigs and external_vectors_df is None:
        raise ValueError("No contigs provided. Please specify either filtered contigs or benchmarks.")

    # Open BAM file to get all references
    with pysam.AlignmentFile(args.input, "rb") as bamfile:
        bam_references = bamfile.references
        if not bam_references:
            raise ValueError("BAM file contains no reference sequences")

    # Filter references to only include our target contigs
    target_references = []
    for ref in bam_references:
        processed_name = extract_first_word(ref)
        if processed_name in bam_contigs:
            target_references.append(ref)

    print(f"Found {len(target_references)} matching references in BAM file")
    if external_vectors_df is not None:
        print(f"Using {len(external_vectors_df)} external reference vectors")

    # Calculate actual dimensions
    length_count = args.max_length - args.min_length + 1
    nucleotides = ['A', 'T', 'C', 'G']
    strands = ['+', '-']
    length_range = range(args.min_length, args.max_length + 1)

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

    print(f"Processing {len(target_references)} target references from BAM file")
    print(f"Counting length range: {args.min_length}-{args.max_length} nt ({length_count} lengths)")
    print(f"Selected mode: {args.mode}")
    print(f"Generated feature vector dimension: {dimensions} dimensions")

    # For strand-specific modes, we'll generate twice as many rows for both BAM contigs and external vectors
    if args.mode in ['sizeXstr', 'sizeX5ntXstr']:
        print(
            f"Note: For {args.mode} mode, each contig (both BAM and external) will generate plus-minus strand versions")

    # Process BAM contigs in parallel
    all_rows = []
    if target_references:
        # Split references into chunks for parallel processing
        chunk_size = max(1, math.ceil(len(target_references) / args.cores))
        args_list = [
            (target_references[i:i + chunk_size], args.input, args.min_length, args.max_length, args.mode)
            for i in range(0, len(target_references), chunk_size)
        ]

        print(f"Using {len(args_list)} processes for processing BAM contigs...")

        # Process chunks in parallel
        with multiprocessing.Pool(processes=args.cores) as pool:
            results = pool.map(process_chunk_optimized, args_list)

        # Merge results from all chunks
        merged_counts = defaultdict(int)
        for chunk_result in results:
            for key, count in chunk_result.items():
                merged_counts[key] += count

        # Generate feature vectors for each reference and filter all-zero rows
        print("Generating feature vectors for BAM contigs...")

        for ref_name in target_references:
            processed_name = extract_first_word(ref_name)
            row = [processed_name]

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

            # Only add if row has non-zero counts
            if sum(row[1:]) > 0:
                all_rows.append(row)

                # For strand-specific modes, also create minus strand version
                if args.mode in ['sizeXstr', 'sizeX5ntXstr']:
                    minus_row = swap_strand_counts(row, args.mode, length_range, nucleotides)
                    if minus_row and sum(minus_row[1:]) > 0:
                        all_rows.append(minus_row)

    # Process external vectors
    if external_vectors_df is not None:
        external_rows = process_external_vectors(external_vectors_df, args.mode, length_range, nucleotides)
        all_rows.extend(external_rows)
        print(f"Added {len(external_rows)} external vector rows")

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

        print(f"Writing {len(all_rows)} total rows to CSV")
        writer.writerow(headers)
        for row in all_rows:
            writer.writerow(row)

    print(f"Processing completed! Results saved to {args.output}")
    print(f"Final feature vector: {dimensions} dimensions")
    print(f"Feature combination: {args.mode}")
    print(f"Length range: {args.min_length}-{args.max_length} nt")

    if args.mode in ['sizeXstr', 'sizeX5ntXstr']:
        print(f"Note: Output contains both plus and minus strand versions for all contigs")


if __name__ == '__main__':
    main()