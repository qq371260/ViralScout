#!/usr/bin/env python3
"""
Filter contigs using correlation analysis of sRNA size vector
"""

import pysam
from collections import defaultdict
import csv
import argparse
import os
import multiprocessing
import math
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from scipy.stats import spearmanr, pearsonr
import warnings
from scipy.stats import combine_pvalues
import re

warnings.filterwarnings('ignore')


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Generate size vectors to find highly correlated contigs')
    parser.add_argument('-i', '--input_bam', type=str, required=True,
                        help='Input BAM file path')
    parser.add_argument('-b', '--benchmarks', type=str, default=None,
                        help='Benchmark contig txt or csv file path. If absent, use default vsiRNA simulant reference')
    parser.add_argument('--correlation_method', type=str, default='spearman',
                        choices=['spearman', 'pearson'],
                        help='Correlation method: spearman or pearson (default: spearman)')
    parser.add_argument('-o', '--output', type=str, default='correlation_results',
                        help='Output file prefix (default: correlation_results)')
    parser.add_argument('--min_length', type=int, default=18,
                        help='Minimum sRNA size')
    parser.add_argument('--max_length', type=int, default=30,
                        help='Maximum sRNA size')
    parser.add_argument('--threads', type=int, default=-1,
                        help='Number of parallel threads for all processing steps (default to use all CPU cores)')
    parser.add_argument('--mean_r', type=float, default=0.8,
                        help='Mean correlation threshold')
    parser.add_argument('--p_value', type=float, default=0.05,
                        help='P-value threshold')
    parser.add_argument('-n', '--number_21_nt', type=int, default=10,
                        help='21nt threshold for filtering')
    parser.add_argument('-m', '--mode', type=str, default='sizeX5nt',
                        choices=['size', 'size_P_5nt', 'sizeXstr', 'sizeX5nt', 'sizeX5ntXstr'],
                        help='Five feature vector modes')
    return parser.parse_args()


def extract_first_word(name):
    first_word = re.split(r'[\s\t,;]', name.strip())[0]
    return first_word


def process_chunk_optimized(chunk_args):
    """ Process a chunk of references to count sRNA reads for size mode only """

    ref_names, bam_path, min_len, max_len = chunk_args
    combined_counts = defaultdict(int)

    with pysam.AlignmentFile(bam_path, "rb") as _bam_file:
        for ref_name in ref_names:
            try:
                # Extract first word for consistent naming
                processed_name = extract_first_word(ref_name)

                for read in _bam_file.fetch(ref_name):
                    # Skip unmapped and secondary alignments
                    if read.is_unmapped or read.is_secondary:
                        continue

                    length = read.query_length
                    if min_len <= length <= max_len:
                        # Only count size mode: length distribution
                        key = (processed_name, length)
                        combined_counts[key] += 1
            except ValueError:
                # Skip references that cannot be found in BAM file
                continue

    return dict(combined_counts)


def generate_feature_vectors(bam_file, output_csv, min_length=18, max_length=30, threads=-1):
    """Generate feature vector data from BAM file - always in size mode"""
    if threads == -1:
        threads = os.cpu_count()
    elif threads is None:
        threads = os.cpu_count()

    print("=" * 60)
    print("Step 1: Generating feature vector data (size mode only)")
    print("=" * 60)

    # Open BAM file to get references
    with pysam.AlignmentFile(bam_file, "rb") as bamfile:
        references = bamfile.references
        if not references:
            raise ValueError("BAM file contains no reference sequences")

    # Calculate actual dimensions
    length_count = max_length - min_length + 1
    dimensions = length_count

    print(f"Found {len(references)} reference sequences")
    print(f"Counting length range: {min_length}-{max_length} nt ({length_count} lengths)")
    print(f"Generated feature vector dimension: {dimensions} dimensions")
    print(f"Using {threads} threads for processing")

    # Split references into chunks for parallel processing
    chunk_size = max(1, math.ceil(len(references) / threads))
    args_list = [
        (references[i:i + chunk_size], bam_file, min_length, max_length)
        for i in range(0, len(references), chunk_size)
    ]

    print(f"Using {len(args_list)} processes for processing...")

    # Process chunks in parallel
    with multiprocessing.Pool(processes=threads) as pool:
        results = pool.map(process_chunk_optimized, args_list)

    # Merge results from all chunks
    merged_counts = defaultdict(int)
    for chunk_result in results:
        for key, count in chunk_result.items():
            merged_counts[key] += count

    # Define feature components
    length_range = range(min_length, max_length + 1)

    # Create a mapping from processed names to original names
    name_mapping = {}
    for ref_name in references:
        processed_name = extract_first_word(ref_name)
        if processed_name not in name_mapping:
            name_mapping[processed_name] = ref_name

    # Write results to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Create header for size mode
        headers = ['Reference']
        for l in length_range:
            headers.append(f'{l}nt')

        # Generate feature vectors for each processed reference name and filter all-zero rows
        print("Generating feature vectors and filtering all-zero rows...")
        non_zero_rows = []

        # Use processed names for the feature matrix
        processed_names = list(name_mapping.keys())

        for processed_name in processed_names:
            row = [processed_name]  # Use processed name as reference

            # Length only (size mode)
            for l in length_range:
                count = merged_counts.get((processed_name, l), 0)
                row.append(count)

            if sum(row[1:]) > 0:
                non_zero_rows.append(row)

        print(
            f"Writing {len(non_zero_rows)} non-zero rows to CSV (filtered out {len(processed_names) - len(non_zero_rows)} all-zero rows)")
        writer.writerow(headers)
        for row in non_zero_rows:
            writer.writerow(row)

        # Save name mapping for reference
        mapping_file = output_csv.replace('.csv', '_name_mapping.csv')
        with open(mapping_file, 'w', newline='') as mapfile:
            map_writer = csv.writer(mapfile)
            map_writer.writerow(['Processed_Name', 'Original_Name'])
            for proc_name, orig_name in name_mapping.items():
                map_writer.writerow([proc_name, orig_name])
        print(f"Name mapping saved to: {mapping_file}")

    print(f"Feature vector data saved to: {output_csv}")
    return output_csv


def load_benchmarks(benchmarks_file, mode='sizeX5nt'):
    """Load benchmark contigs from file or use default vsiRNA simulants"""

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
                ref_df = pd.read_csv(ref_file, header=0)
                # Return the whole DataFrame as we need vector data
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


def convert_vectors_to_size_mode(vectors_df, original_mode, min_length=18, max_length=30):
    """Convert vectors from any mode to size mode"""
    print(f"Converting vectors from {original_mode} mode to size mode...")

    length_range = range(min_length, max_length + 1)
    nucleotides = ['A', 'T', 'C', 'G']
    strands = ['+', '-']

    # Create new DataFrame for size mode vectors
    size_vectors = []

    # Handle both DataFrame and Series input
    if isinstance(vectors_df, pd.DataFrame):
        rows_to_process = vectors_df.iterrows()
        index_col = vectors_df.columns[0]
    else:
        # Should not happen in our usage
        return None

    for idx, row in rows_to_process:
        contig_name = row[index_col]
        size_row = [contig_name]

        if original_mode == 'size':
            # Already in size mode, just copy
            for l in length_range:
                col_name = f'{l}nt'
                size_row.append(row[col_name])

        elif original_mode == 'size_P_5nt':
            # Take only the length columns, ignore base columns
            for l in length_range:
                col_name = f'{l}nt'
                size_row.append(row[col_name])

        elif original_mode == 'sizeXstr':
            # Sum plus and minus strands for each length
            for l in length_range:
                plus_col = f'{l}nt_plus'
                minus_col = f'{l}nt_minus'
                total = row.get(plus_col, 0) + row.get(minus_col, 0)
                size_row.append(total)

        elif original_mode == 'sizeX5nt':
            # Sum all nucleotides for each length
            for l in length_range:
                total = 0
                for nt in nucleotides:
                    col_name = f'{l}nt_{nt}'
                    total += row.get(col_name, 0)
                size_row.append(total)

        elif original_mode == 'sizeX5ntXstr':
            # Sum all nucleotides and strands for each length
            for l in length_range:
                total = 0
                for nt in nucleotides:
                    for strand in strands:
                        strand_name = 'plus' if strand == '+' else 'minus'
                        col_name = f'{l}nt_{nt}_{strand_name}'
                        total += row.get(col_name, 0)
                size_row.append(total)

        size_vectors.append(size_row)

    # Create header for size vectors
    headers = ['Reference']
    for l in length_range:
        headers.append(f'{l}nt')

    size_df = pd.DataFrame(size_vectors, columns=headers)
    size_df = size_df.set_index('Reference')

    return size_df


def calculate_correlation_and_pvalue(target_row, other_row, df, numeric_cols, corr_method='spearman'):
    """Calculate correlation coefficient and p-value between two rows"""
    try:
        y = df.loc[target_row, numeric_cols].values.astype(float)
        x = df.loc[other_row, numeric_cols].values.astype(float)

        valid_mask = ~(np.isnan(x) | np.isnan(y))
        x_valid = x[valid_mask]
        y_valid = y[valid_mask]

        if len(x_valid) <= 3:
            return 0.0, 1.0

        # Use selected correlation method
        if corr_method == 'pearson':
            corr, p_value = pearsonr(x_valid, y_valid)
        else:  # spearman
            corr, p_value = spearmanr(x_valid, y_valid)

        if np.isnan(corr):
            return 0.0, 1.0

        return corr, p_value

    except Exception as e:
        print(f"Error calculating correlation between {target_row} and {other_row}: {e}")
        return 0.0, 1.0


def get_top_k_with_p_value(target_row, k, df, numeric_cols, corr_method='spearman'):
    """Get top k rows with highest correlation to target row and their p-values"""
    correlations = {}
    p_values = {}

    for other_row in df.index:
        if other_row == target_row:
            continue

        corr, p_value = calculate_correlation_and_pvalue(target_row, other_row, df, numeric_cols, corr_method)
        if corr > 0:
            correlations[other_row] = corr
            p_values[other_row] = p_value

    sorted_rows = sorted(correlations.items(), key=lambda x: x[1], reverse=True)[:k]
    return [(row, correlations[row], p_values[row]) for row, _ in sorted_rows]


def filter_by_21nt_threshold(df, number_21_nt, benchmark_set=None):
    """Filter dataframe based on 21nt total threshold, but always keep benchmark contigs"""
    print(f"Checking and filtering by 21nt threshold: {number_21_nt}")

    # If no benchmark set provided, create empty set
    if benchmark_set is None:
        benchmark_set = set()

    print(f"Preserving {len(benchmark_set)} benchmark contigs from filtering")

    # Check if 21nt column exists
    if '21nt' in df.columns:
        # Directly use 21nt column for filtering
        initial_count = len(df)

        # Create filter condition: keep if 21nt >= threshold OR if it's a benchmark contig
        condition = (df['21nt'] >= number_21_nt) | (df.index.isin(benchmark_set))
        df_filtered = df[condition]

        filtered_count = len(df_filtered)
        benchmark_preserved = sum(df_filtered.index.isin(benchmark_set))

        print(f"Filtered by 21nt threshold {number_21_nt}: {initial_count} -> {filtered_count} rows")
        print(f"Preserved {benchmark_preserved} benchmark contigs")
        return df_filtered
    else:
        # Check for 21nt-related columns
        twentyone_nt_cols = [col for col in df.columns if '21' in str(col)]
        if len(twentyone_nt_cols) > 0:
            print(f"Found 21nt-related columns: {twentyone_nt_cols}")
            # Calculate sum of 21nt-related columns
            df['21nt_total'] = df[twentyone_nt_cols].sum(axis=1)
            initial_count = len(df)

            # Create filter condition: keep if 21nt_total >= threshold OR if it's a benchmark contig
            condition = (df['21nt_total'] >= number_21_nt) | (df.index.isin(benchmark_set))
            df_filtered = df[condition]

            filtered_count = len(df_filtered)
            benchmark_preserved = sum(df_filtered.index.isin(benchmark_set))

            print(f"Filtered by 21nt total threshold {number_21_nt}: {initial_count} -> {filtered_count} rows")
            print(f"Preserved {benchmark_preserved} benchmark contigs")
            return df_filtered
        else:
            print(f"Warning: No 21nt column found, skipping 21nt filtering")
            return df


def find_highly_correlated_contigs(vector_csv, benchmarks_file, output_prefix,
                                   mean_r=0.8, p_value=0.05, number_21_nt=10,
                                   threads=-1, mode='sizeX5nt', min_length=18, max_length=30,
                                   correlation_method='spearman'):
    """Find highly correlated contigs using selected correlation method"""

    print("\n" + "=" * 60)
    print(f"Step 2: Finding highly correlated contigs using {correlation_method} correlation")
    print("=" * 60)

    # Set fixed parameters
    specific_contigs = 1000  # K = 1000
    common_contigs = 999  # L = 999 (K > L)

    print(f"Using fixed parameters: specific_contigs(K)={specific_contigs}, common_contigs(L)={common_contigs}")
    print(f"Using {threads} threads for correlation analysis")
    print(f"Correlation method: {correlation_method}")

    # Reading feature vector data (size mode)
    print("Reading feature vector data...")
    df = pd.read_csv(vector_csv)
    df = df.set_index(df.columns[0])

    # Loading benchmarks
    benchmark_names, benchmark_vectors_df = load_benchmarks(benchmarks_file, mode)

    if benchmark_names is None and benchmark_vectors_df is None:
        print("Error: Cannot load benchmark contigs")
        return

    # Initialize benchmark_set to track which contigs are benchmarks
    benchmark_set = set()

    # Process benchmark data
    if benchmark_vectors_df is not None:
        # Using default vsiRNA file, need to convert to size mode
        print("Using default vsiRNA simulant vectors")
        benchmark_vectors_size = convert_vectors_to_size_mode(
            benchmark_vectors_df, mode,
            min_length=min_length, max_length=max_length
        )
        # Merge benchmark vectors with main data
        df = pd.concat([df, benchmark_vectors_size])
        target_rows = benchmark_vectors_size.index.tolist()
        benchmark_set = set(target_rows)
    else:
        # Using user-provided benchmark names
        target_rows = [row for row in benchmark_names if row in df.index]
        benchmark_set = set(target_rows)

    # Validate target row availability
    valid_target_rows = [row for row in target_rows if row in df.index]
    print(f"Valid benchmark contigs: {len(valid_target_rows)}/{len(target_rows)}")
    if len(valid_target_rows) < len(target_rows):
        missing = set(target_rows) - set(valid_target_rows)
        print(f"Missing benchmark contigs: {missing}")

    if not valid_target_rows:
        print("Error: No valid benchmark contigs found in the vector data")
        return

    # Filter by 21nt threshold, but preserve benchmark contigs
    df = filter_by_21nt_threshold(df, number_21_nt, benchmark_set)

    # Check if benchmark contigs are still present after filtering
    remaining_benchmarks = [row for row in valid_target_rows if row in df.index]
    if len(remaining_benchmarks) < len(valid_target_rows):
        lost_benchmarks = set(valid_target_rows) - set(remaining_benchmarks)
        print(f"WARNING: {len(lost_benchmarks)} benchmark contigs lost during filtering: {lost_benchmarks}")
        valid_target_rows = remaining_benchmarks

    if not valid_target_rows:
        print("Error: All benchmark contigs were lost during filtering!")
        return

    # Get all numeric columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    print(f"Using {len(numeric_cols)} feature dimensions")

    # Parallel processing for each target row
    print(f"Calculating {correlation_method} correlation for feature vectors...")
    top_k_with_p_values = Parallel(n_jobs=threads)(
        delayed(get_top_k_with_p_value)(target_row, specific_contigs, df, numeric_cols, correlation_method)
        for target_row in valid_target_rows
    )

    # Organize results as dictionary
    row_rankings_with_p_values = {
        target_row: {row: (corr, p_value) for row, corr, p_value in items}
        for target_row, items in zip(valid_target_rows, top_k_with_p_values)
    }

    # Check if we have valid results
    if not row_rankings_with_p_values:
        print("Error: No valid correlation results found")
        return

    # Calculate common rows
    common_rows = set.intersection(*[set(d.keys()) for d in row_rankings_with_p_values.values()])
    print(f"Found {len(common_rows)} common highly correlated contigs")

    # Calculate composite score and combined p-value
    scoring = defaultdict(lambda: {'score': 0.0, 'p_value': 1.0})
    for row in common_rows:
        p_values = []
        for target_row in valid_target_rows:
            corr, individual_p_value = row_rankings_with_p_values[target_row][row]
            scoring[row]['score'] += corr
            p_values.append(individual_p_value)

        scoring[row]['score'] /= len(valid_target_rows)

        if p_values:
            clean_p_values = [max(min(p, 0.999), 0.001) for p in p_values]
            try:
                _, combined_p = combine_pvalues(clean_p_values, method='fisher')
                scoring[row]['p_value'] = combined_p
            except:
                scoring[row]['p_value'] = min(1.0, len(p_values) * min(p_values))

    # Take top L results, sorted by composite score in descending order
    final_rows = sorted(scoring.items(), key=lambda x: x[1]['score'], reverse=True)[:common_contigs]

    # Apply threshold filtering
    filtered_rows = [(row, values['score'], values['p_value'])
                     for row, values in final_rows
                     if values['score'] >= mean_r and values['p_value'] < p_value]

    print(
        f"After threshold filtering (mean correlation>={mean_r}, p-value<{p_value}), retained {len(filtered_rows)} contigs")

    # Check if any contigs passed filtering
    if len(filtered_rows) == 0:
        print("\n" + "=" * 60)
        print("WARNING: No contigs passed the threshold filtering!")
        print("=" * 60)
        return

    # Save filtered contig names
    output_file = f"{output_prefix}_highly_correlated_contigs.txt"
    with open(output_file, 'w') as f:
        for row, score, p_val in filtered_rows:
            f.write(f"{row}\n")

    # Save detailed results with appropriate column names
    detailed_file = f"{output_prefix}_detailed_results.csv"

    # Use appropriate column name based on correlation method
    if correlation_method == 'pearson':
        corr_column_name = "Mean_Pearson_Correlation"
    else:
        corr_column_name = "Mean_Spearman_Correlation"

    df_detailed = pd.DataFrame(filtered_rows, columns=["Contig_Name", corr_column_name, "P_Value"])
    df_detailed.to_csv(detailed_file, index=False)

    print(f"\nResults saved:")
    print(f"  - {output_file}: Filtered contig names")
    print(f"  - {detailed_file}: Detailed results with correlation scores and p-values")

    # Print top results
    print(f"\nTop 10 highly correlated contigs:")
    for i, (row, score, p_val) in enumerate(filtered_rows[:10], 1):
        print(f"  {i}. {row} (mean correlation: {score:.4f}, p-value: {p_val:.4e})")


def main():
    # Parse arguments
    args = parse_arguments()

    # Step 1: Generate feature vectors (always use size mode for main analysis)
    vector_csv = f"{args.output}_vectors.csv"
    generate_feature_vectors(
        bam_file=args.input_bam,
        output_csv=vector_csv,
        min_length=args.min_length,
        max_length=args.max_length,
        threads=args.threads
    )

    # Step 2: Find highly correlated contigs
    find_highly_correlated_contigs(
        vector_csv=vector_csv,
        benchmarks_file=args.benchmarks,
        output_prefix=args.output,
        mean_r=args.mean_r,
        p_value=args.p_value,
        number_21_nt=args.number_21_nt,
        threads=args.threads,
        mode=args.mode,
        min_length=args.min_length,
        max_length=args.max_length,
        correlation_method=args.correlation_method
    )

    print("\n" + "=" * 60)
    print("Analysis completed successfully!")
    print("=" * 60)


if __name__ == '__main__':
    main()
