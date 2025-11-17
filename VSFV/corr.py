#!/usr/bin/env python3
"""
Correlation analysis of sRNAs feature vectors of viral and other contigs
"""

import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from collections import defaultdict
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage
import warnings
import argparse
import os
from scipy.spatial.distance import squareform
import re
import math
import matplotlib.patches as patches
from matplotlib import font_manager

warnings.filterwarnings('ignore')
import matplotlib.colors as mcolors
from scipy.stats import combine_pvalues


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Correlation-based feature vector clustering analysis',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Input CSV file path containing feature vector data')
    parser.add_argument('-b', '--benchmarks', type=str, default=None,
                        help='Benchmark contig txt or csv file path. If absent, use default vsiRNA simulant reference')
    parser.add_argument('--correlation_method', type=str, default='spearman',
                        choices=['spearman', 'pearson'],
                        help='Correlation method: spearman or pearson')
    parser.add_argument('-m', '--mean_r', type=float, default=0.8,
                        help='Mean correlation threshold')
    parser.add_argument('-p', '--p_value', type=float, default=0.05,
                        help='P-value threshold')
    parser.add_argument('-o', '--output', type=str, default='',
                        help='Output file prefix (default: use vector mode as prefix)')
    parser.add_argument('-t', '--threads', type=int, default=-1,
                        help='Number of parallel threads (default to use all CPU cores)')
    parser.add_argument('--min_length', type=int, default=18,
                        help='Minimum sRNA size for dimension calculation')
    parser.add_argument('--max_length', type=int, default=30,
                        help='Maximum sRNA size for dimension calculation')
    return parser.parse_args()


def load_benchmarks(benchmarks_arg, vector_mode=None):
    """Load benchmark contigs from -b parameter or a default reference file with multi-vector-modes
    """

    def extract_first_word(text):
        match = re.search(r'[^\s\t,;]+', text.strip())
        return match.group() if match else ""

    if benchmarks_arg is None:
        # Use default vsiRNA simulant reference with dynamic mode
        default_filename = f"vsiRNA_simulant_{vector_mode}.csv"
        possible_locations = [
            default_filename,
            f"../{default_filename}",
            f"./{default_filename}"
        ]

        for ref_file in possible_locations:
            if os.path.isfile(ref_file):
                print(f"Using default vsiRNA simulant file: {ref_file}")
                ref_df = pd.read_csv(ref_file)
                raw_benchmarks = ref_df.iloc[:, 0].tolist()  # First column contains names
                benchmarks = [extract_first_word(name) for name in raw_benchmarks]
                print(f"Loaded {len(benchmarks)} benchmark contigs from default vsiRNA simulant file")
                return benchmarks

        print(f"Error: No benchmarks provided and default vsiRNA simulant files not found")
        print(f"Expected files: {default_filename}")
        print("Checked locations:")
        for loc in possible_locations:
            print(f"  - {loc}")
        return []

    # Check if it's a file path
    possible_paths = [
        benchmarks_arg,
        f"../{benchmarks_arg}",
        f"./{benchmarks_arg}"
    ]

    actual_path = None
    for path in possible_paths:
        if os.path.isfile(path):
            actual_path = path
            break

    if actual_path:
        # Check file extension
        if actual_path.endswith('.csv'):
            print(f"Loading benchmark contigs from CSV file: {actual_path}")
            ref_df = pd.read_csv(actual_path)
            raw_benchmarks = ref_df.iloc[:, 0].tolist()  # First column contains names
            benchmarks = [extract_first_word(name) for name in raw_benchmarks]
            print(f"Loaded {len(benchmarks)} benchmark contigs from CSV file")
            return benchmarks
        else:
            # Assume txt file with one contig name per line
            print(f"Loading benchmark contigs from text file: {actual_path}")
            with open(actual_path, 'r') as f:
                raw_benchmarks = [line.strip() for line in f if line.strip()]
                benchmarks = [extract_first_word(name) for name in raw_benchmarks]
            print(f"Loaded {len(benchmarks)} benchmark contigs from text file")
            return benchmarks
    else:
        # Direct input of contig names (space separated) or file not found
        if ' ' in benchmarks_arg or '/' in benchmarks_arg or '\\' in benchmarks_arg:
            # Looks like file path but not found
            print(f"Error: Benchmark file not found: {benchmarks_arg}")
            print("Checked locations:")
            for path in possible_paths:
                print(f"  - {path}")
            return []
        else:
            # Direct input of contig names
            raw_benchmarks = benchmarks_arg.split()
            benchmarks = [extract_first_word(name) for name in raw_benchmarks]
            print(f"Direct input of {len(benchmarks)} benchmark contigs")
            return benchmarks


def detect_vector_mode(numeric_cols, min_length, max_length):
    """
    Detect vector mode based on number of columns and length range

    Modes:
    - size: length only (L dimensions)
    - size_P_5nt: length + 5' nucleotides (L + 4 dimensions)
    - sizeXstr: length × strand (L × 2 dimensions)
    - sizeX5nt: length × 5' nucleotides (L × 4 dimensions)
    - sizeX5ntXstr: length × 5' nucleotides × strand (L × 8 dimensions)
    """
    num_cols = len(numeric_cols)
    length_count = max_length - min_length + 1

    # Calculate expected dimensions for each mode
    expected_dims = {
        'size': length_count,
        'size_P_5nt': length_count + 4,
        'sizeXstr': length_count * 2,
        'sizeX5nt': length_count * 4,
        'sizeX5ntXstr': length_count * 8
    }

    # Find exact match
    for mode, expected_dim in expected_dims.items():
        if num_cols == expected_dim:
            return mode

    # If no exact match, find the closest
    differences = {mode: abs(num_cols - dim) for mode, dim in expected_dims.items()}
    best_match = min(differences, key=differences.get)

    print(f"Warning: Expected {expected_dims[best_match]} columns for {best_match} mode, but found {num_cols}")
    print(f"Using closest match: {best_match}")
    return best_match


def calculate_correlation(corr_method, x, y):
    """Calculate correlation using specified method"""
    if corr_method == 'spearman':
        return spearmanr(x, y)
    else:  # pearson
        return pearsonr(x, y)


def process_strand_data(df, valid_target_rows, numeric_cols, corr_method, args):
    """
    Process strand-specific data with unified polarity
    Returns processed dataframe with only selected strands
    """
    print("Processing strand-specific data with unified polarity...")

    # Create a copy to avoid modifying original
    df_processed = df.copy()

    # Extract base names (remove _minus suffix if present)
    base_names = {}
    for contig in df_processed.index:
        if contig.endswith('_minus'):
            base_name = contig[:-6]
        else:
            base_name = contig
        if base_name not in base_names:
            base_names[base_name] = []
        base_names[base_name].append(contig)

    # Step 1: Select reference benchmark (first one, prefer positive strand)
    if not valid_target_rows:
        print("No valid target rows found")
        return df_processed, valid_target_rows

    selected_benchmarks = []

    # For first benchmark, prefer positive strand
    first_benchmark = valid_target_rows[0]
    first_benchmark_base = first_benchmark
    if first_benchmark_base.endswith('_minus'):
        first_benchmark_base = first_benchmark_base[:-6]

    if first_benchmark_base in base_names:
        available_strands = base_names[first_benchmark_base]
        # Prefer positive strand (without _minus)
        positive_strand = None
        minus_strand = None
        for strand in available_strands:
            if strand.endswith('_minus'):
                minus_strand = strand
            else:
                positive_strand = strand

        if positive_strand and positive_strand in df_processed.index:
            selected_benchmarks.append(positive_strand)
            print(f"First benchmark: selected positive strand '{positive_strand}'")
        elif minus_strand and minus_strand in df_processed.index:
            selected_benchmarks.append(minus_strand)
            print(f"First benchmark: using minus strand '{minus_strand}' (no positive available)")
        else:
            # Fallback: use original if available
            if first_benchmark in df_processed.index:
                selected_benchmarks.append(first_benchmark)
                print(f"First benchmark: using original '{first_benchmark}' (fallback)")
    else:
        # Fallback: use original if available
        if first_benchmark in df_processed.index:
            selected_benchmarks.append(first_benchmark)
            print(f"First benchmark: using original '{first_benchmark}' (fallback)")

    # Step 2: For other benchmarks, select strand with highest correlation to first benchmark
    for i, benchmark in enumerate(valid_target_rows[1:], 1):
        benchmark_base = benchmark
        if benchmark_base.endswith('_minus'):
            benchmark_base = benchmark_base[:-6]

        if benchmark_base in base_names:
            available_strands = base_names[benchmark_base]
            strand_scores = {}

            for strand in available_strands:
                if strand in df_processed.index and selected_benchmarks:  # Ensure we have reference
                    try:
                        # Calculate correlation with first benchmark
                        first_bench_vector = df_processed.loc[selected_benchmarks[0], numeric_cols].values.astype(float)
                        strand_vector = df_processed.loc[strand, numeric_cols].values.astype(float)

                        valid_mask = ~(np.isnan(first_bench_vector) | np.isnan(strand_vector))
                        first_valid = first_bench_vector[valid_mask]
                        strand_valid = strand_vector[valid_mask]

                        if len(first_valid) > 3:
                            corr, _ = calculate_correlation(corr_method, first_valid, strand_valid)
                            if not np.isnan(corr):
                                strand_scores[strand] = corr
                    except Exception as e:
                        print(f"Error calculating correlation for {strand}: {e}")

            if strand_scores:
                # Select strand with highest correlation to first benchmark
                best_strand = max(strand_scores.items(), key=lambda x: x[1])
                selected_benchmarks.append(best_strand[0])
                print(f"Benchmark {i}: selected '{best_strand[0]}' (corr: {best_strand[1]:.4f})")
            else:
                # Fallback: use positive strand if available
                positive_strand = None
                for strand in available_strands:
                    if not strand.endswith('_minus') and strand in df_processed.index:
                        positive_strand = strand
                        break
                if positive_strand:
                    selected_benchmarks.append(positive_strand)
                    print(f"Benchmark {i}: using positive strand '{positive_strand}' (fallback)")
                elif available_strands and available_strands[0] in df_processed.index:
                    selected_benchmarks.append(available_strands[0])
                    print(f"Benchmark {i}: using '{available_strands[0]}' (fallback)")
        else:
            # Fallback: use original if available
            if benchmark in df_processed.index:
                selected_benchmarks.append(benchmark)
                print(f"Benchmark {i}: using original '{benchmark}' (fallback)")

    print(f"Selected {len(selected_benchmarks)} benchmark strands (reduced from {len(valid_target_rows)})")

    # Step 3: Select strands for non-benchmark contigs based on correlation with selected benchmarks
    print("Selecting strands for non-benchmark contigs...")

    # Create a set of benchmark base names for quick lookup
    benchmark_bases = set()
    for benchmark in selected_benchmarks:
        base_name = benchmark.rstrip('_minus')
        benchmark_bases.add(base_name)

    selected_contigs = []

    for base_contig, strands in base_names.items():
        # If it's a benchmark, directly use the strand selected in steps 1 and 2
        if base_contig in benchmark_bases:
            # Find the corresponding benchmark strand
            for strand in strands:
                if strand in selected_benchmarks:
                    selected_contigs.append(strand)
                    break
            else:
                # If not found, use the first strand
                selected_contigs.append(strands[0])
                print(f"  Warning: Benchmark {base_contig} strand not found in selected_benchmarks")
        else:
            # Non-benchmark contigs: select strand based on correlation
            if len(strands) == 1:
                selected_contigs.append(strands[0])
            else:
                strand_scores = {}

                for strand in strands:
                    if strand in df_processed.index:
                        try:
                            strand_vector = df_processed.loc[strand, numeric_cols].values.astype(float)
                            total_corr = 0.0
                            valid_count = 0

                            for benchmark in selected_benchmarks:
                                if benchmark in df_processed.index:
                                    bench_vector = df_processed.loc[benchmark, numeric_cols].values.astype(float)

                                    valid_mask = ~(np.isnan(strand_vector) | np.isnan(bench_vector))
                                    strand_valid = strand_vector[valid_mask]
                                    bench_valid = bench_vector[valid_mask]

                                    if len(strand_valid) > 3:
                                        corr, _ = calculate_correlation(corr_method, strand_valid, bench_valid)
                                        if not np.isnan(corr):
                                            total_corr += corr
                                            valid_count += 1

                            if valid_count > 0:
                                avg_corr = total_corr / valid_count
                                strand_scores[strand] = avg_corr
                            else:
                                strand_scores[strand] = 0.0

                        except Exception as e:
                            print(f"Error processing {strand}: {e}")
                            strand_scores[strand] = 0.0

                # Select the strand with the highest correlation
                if strand_scores:
                    best_strand = max(strand_scores.items(), key=lambda x: x[1])
                    selected_contigs.append(best_strand[0])
                else:
                    # Fallback: use positive strand
                    positive_strand = None
                    for strand in strands:
                        if not strand.endswith('_minus'):
                            positive_strand = strand
                            break
                    if positive_strand:
                        selected_contigs.append(positive_strand)
                    else:
                        selected_contigs.append(strands[0])

    # Create new dataframe with selected contigs only
    df_final = df_processed.loc[selected_contigs].copy()

    # Update target rows to only include selected benchmark strands
    final_target_rows = [benchmark for benchmark in selected_benchmarks if benchmark in df_final.index]

    print(f"Strand processing complete:")
    print(f"  - Original contigs: {len(df_processed)}")
    print(f"  - Final contigs: {len(df_final)}")
    print(f"  - Original benchmarks: {len(valid_target_rows)}")
    print(f"  - Final benchmarks: {len(final_target_rows)}")
    print(f"  - Reduction factor: {len(df_final) / len(df_processed):.2f}")

    return df_final, final_target_rows


def main():
    # Parse arguments
    args = parse_arguments()

    # Reading and detect vector mode
    print("Reading feature vector data...")
    df = pd.read_csv(args.input)
    df = df.set_index(df.columns[0])

    # Get all numeric columns
    numeric_cols = [col for col in df.select_dtypes(include=[np.number]).columns.tolist()]

    # Auto-detect vector mode
    vector_mode = detect_vector_mode(numeric_cols, args.min_length, args.max_length)

    # Loading benchmarks using vector_mode
    target_rows = load_benchmarks(args.benchmarks, vector_mode)
    if not target_rows:
        print("Error: Cannot load benchmark contigs")
        return

    # Output name prefix
    output_prefix = args.output if args.output else vector_mode

    print("=" * 60)
    print(f"{args.correlation_method.capitalize()} Correlation Analysis")
    print("=" * 60)
    print(f"Input file: {args.input}")
    print(f"Detected vector mode: {vector_mode}")
    print(f"Correlation method: {args.correlation_method}")
    print(f"Length range: {args.min_length}-{args.max_length} nt")
    print(f"Actual dimension: {len(numeric_cols)}")
    print(f"Benchmark contigs: {len(target_rows)}")
    print(f"Mean correlation threshold: {args.mean_r}")
    print(f"P-value threshold: {args.p_value}")
    print(f"Output prefix: {output_prefix}")
    print(f"Threads: {args.threads}")
    print("=" * 60)

    # Read vector data (if not already read)
    if 'df' not in locals():
        df = pd.read_csv(args.input)
        df = df.set_index(df.columns[0])

    # Validate target row availability
    valid_target_rows = [row for row in target_rows if row in df.index]
    print(f"Valid benchmark contigs: {len(valid_target_rows)}/{len(target_rows)}")
    if len(valid_target_rows) < len(target_rows):
        missing = set(target_rows) - set(valid_target_rows)
        print(f"Missing benchmark contigs: {missing}")

    # Get all numeric columns
    numeric_cols = [col for col in df.select_dtypes(include=[np.number]).columns.tolist()]

    print(f"Using vector mode: {vector_mode}")
    print(f"Feature vector analysis - using {len(numeric_cols)} feature dimensions")

    # Feature description based on mode
    length_count = args.max_length - args.min_length + 1
    if vector_mode == 'size':
        print(
            f"Feature description: {length_count}-dimensional length distribution ({args.min_length}-{args.max_length}nt)")
    elif vector_mode == 'size_P_5nt':
        print(f"Feature description: {length_count + 4}-dimensional (length distribution + A/T/C/G base composition)")
    elif vector_mode == 'sizeXstr':
        print(f"Feature description: {length_count * 2}-dimensional (length × strand specificity)")
    elif vector_mode == 'sizeX5nt':
        print(f"Feature description: {length_count * 4}-dimensional (length × A/T/C/G base specificity)")
    elif vector_mode == 'sizeX5ntXstr':
        print(f"Feature description: {length_count * 8}-dimensional (length × A/T/C/G × strand specificity)")

    # Initialize analysis variables
    strand_specific_mode = vector_mode in ['sizeXstr', 'sizeX5ntXstr']

    # For strand-specific modes, process strand data
    if strand_specific_mode:
        df_processed, final_target_rows = process_strand_data(df, valid_target_rows, numeric_cols,
                                                              args.correlation_method, args)

        analysis_df = df_processed
        final_target_indices = [row for row in final_target_rows if row in analysis_df.index]

        print(f"Strand mode filtering: {len(analysis_df)} non-redundant contigs")
    else:
        print("Non-strand mode, using original data")
        analysis_df = df
        final_target_indices = valid_target_rows


    def calculate_correlation_and_pvalue(target_row, other_row, df, numeric_cols, corr_method):
        """Calculate correlation coefficient and p-value between two rows"""
        try:
            y = df.loc[target_row, numeric_cols].values.astype(float)
            x = df.loc[other_row, numeric_cols].values.astype(float)

            valid_mask = ~(np.isnan(x) | np.isnan(y))
            x_valid = x[valid_mask]
            y_valid = y[valid_mask]

            if len(x_valid) <= 3:
                return 0.0, 1.0

            corr, p_value = calculate_correlation(corr_method, x_valid, y_valid)

            if np.isnan(corr):
                return 0.0, 1.0

            return corr, p_value

        except Exception as e:
            print(f"Error calculating correlation between {target_row} and {other_row}: {e}")
            return 0.0, 1.0

    def calculate_all_correlations(target_row, df, numeric_cols, corr_method):
        """Calculate correlations between target row and all other rows"""
        correlations = {}
        p_values = {}

        for other_row in df.index:
            if other_row == target_row:
                continue

            corr, p_value = calculate_correlation_and_pvalue(target_row, other_row, df, numeric_cols, corr_method)
            correlations[other_row] = corr
            p_values[other_row] = p_value

        return correlations, p_values

    # Calculate correlations for all benchmarks
    print(f"Calculating {args.correlation_method} correlations for all contigs...")
    all_correlations = Parallel(n_jobs=args.threads)(
        delayed(calculate_all_correlations)(target_row, analysis_df, numeric_cols, args.correlation_method)
        for target_row in final_target_indices
    )

    # Combine results from all benchmarks
    combined_scores = defaultdict(lambda: {'correlations': [], 'p_values': []})

    for i, (correlations, p_values) in enumerate(all_correlations):
        target_row = final_target_indices[i]
        for contig, corr in correlations.items():
            combined_scores[contig]['correlations'].append(corr)
            combined_scores[contig]['p_values'].append(p_values[contig])

    # Calculate mean correlation and combined p-value for each contig
    filtered_rows = []
    for contig, scores in combined_scores.items():
        if scores['correlations']:
            mean_corr = np.mean(scores['correlations'])

            # Combine p-values using Fisher's method
            clean_p_values = [max(min(p, 0.999), 0.001) for p in scores['p_values']]
            try:
                _, combined_p = combine_pvalues(clean_p_values, method='fisher')
            except:
                combined_p = min(1.0, len(clean_p_values) * min(clean_p_values))

            # Apply threshold filtering
            if mean_corr >= args.mean_r and combined_p < args.p_value:
                filtered_rows.append((contig, mean_corr, combined_p))

    # Sort by mean correlation in descending order
    filtered_rows.sort(key=lambda x: x[1], reverse=True)

    print(
        f"After threshold filtering (mean correlation>={args.mean_r}, p-value<{args.p_value}), retained {len(filtered_rows)} contigs")

    # Quit if no contigs retained
    if len(filtered_rows) == 0:
        print("\n" + "=" * 60)
        print("WARNING: No contigs passed the threshold filtering!")
        print("=" * 60)
        print(f"Analysis completed but no contigs met the criteria:")
        print(f"  - Mean correlation threshold: {args.mean_r}")
        print(f"  - P-value threshold: {args.p_value}")

        try:
            with open(f'{output_prefix}_analysis_summary.txt', 'w') as f:
                f.write(f"{vector_mode} {args.correlation_method.capitalize()} Correlation Analysis Summary\n")
                f.write("=" * 50 + "\n")
                f.write("WARNING: No contigs passed the threshold filtering!\n\n")
                f.write("Filtering Criteria:\n")
                f.write(f"  Mean correlation threshold: {args.mean_r}\n")
                f.write(f"  P-value threshold: {args.p_value}\n\n")
                f.write("Analysis Parameters:\n")
                f.write(f"  Input file: {args.input}\n")
                f.write(f"  Detected vector mode: {vector_mode}\n")
                f.write(f"  Correlation method: {args.correlation_method}\n")
                f.write(f"  Benchmark contigs: {len(target_rows)}\n")
                f.write(f"  Valid benchmark contigs: {len(final_target_indices)}\n")

            print(f"Analysis summary saved as '{output_prefix}_analysis_summary.txt'")
        except Exception as e:
            print(f"Error saving analysis summary: {e}")

        print("\nProgram terminated: No valid contigs found after filtering.")
        return

    # Extract filtered row names (excluding benchmarks)
    filtered_row_indices = [row for row, score, p_value in filtered_rows if row not in final_target_indices]
    print(f"After excluding benchmarks, {len(filtered_row_indices)} contigs remain for output")

    # Combine target rows and filtered rows for clustering
    combined_row_indices = filtered_row_indices + final_target_indices
    combined_row_indices = list(set(combined_row_indices))

    # Ensure all rows exist in dataframe
    available_indices = [idx for idx in combined_row_indices if idx in analysis_df.index]
    print(f"Total rows for heatmap analysis: {len(available_indices)}")

    available_names = available_indices
    final_target_names = final_target_indices

    # Strict correlation-based clustering
    print(f"Starting strict {args.correlation_method} correlation-based clustering analysis...")

    # Extract feature data for filtered rows
    feature_data = analysis_df.loc[available_indices, numeric_cols].copy()

    # Normalize each row (value/total, i.e., calculate proportion)
    print("Normalizing row data (calculating proportions)...")
    row_sums = feature_data.sum(axis=1)
    normalized_data = feature_data.div(row_sums, axis=0)

    print(f"Normalized data shape: {normalized_data.shape}")
    print(f"Normalized data range: min={normalized_data.values.min():.6f}, max={normalized_data.values.max():.6f}")

    # Calculate distance matrix (1 - correlation)
    print(f"Calculating {args.correlation_method} distance matrix...")
    n_samples = len(available_indices)
    distance_matrix = np.zeros((n_samples, n_samples))

    for i in range(n_samples):
        for j in range(i, n_samples):
            if i == j:
                distance_matrix[i, j] = 0.0
            else:
                idx1 = available_indices[i]
                idx2 = available_indices[j]
                data1 = normalized_data.iloc[i].values
                data2 = normalized_data.iloc[j].values

                # Create valid data mask
                valid_mask = ~(np.isnan(data1) | np.isnan(data2))
                data1_valid = data1[valid_mask]
                data2_valid = data2[valid_mask]

                if len(data1_valid) > 3:
                    corr, _ = calculate_correlation(args.correlation_method, data1_valid, data2_valid)
                    # Convert correlation to distance: distance = 1 - correlation
                    distance = 1 - corr
                    if np.isnan(distance):
                        distance = 2.0
                else:
                    distance = 2.0  # Maximum distance

                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance

    print(f"Distance matrix range: {distance_matrix.min():.3f} to {distance_matrix.max():.3f}")

    # Hierarchical clustering using distance matrix
    print("Performing hierarchical clustering...")

    # distance_matrix is n x n symmetric
    condensed_dist = squareform(distance_matrix, checks=False)
    row_linkage = linkage(condensed_dist, method='average')

    # Visualization section
    print("Plotting heatmap and dendrogram using clustermap...")

    # Figure size (automatically adjusted based on number of rows and columns)
    fig_size_width = max(10, len(numeric_cols) * 0.1)
    fig_size_height = max(8, len(available_indices) * 0.1)
    print(f"Setting figure size: {fig_size_width:.1f} × {fig_size_height:.1f} inches")

    # Custom color map
    colors = ['#F0F0F0', '#E0F0E0', '#C8E6C9', '#A5D6A7', '#81C784',
              '#66BB6A', '#4CAF50', '#388E3C', '#2E7D32', '#1B5E20']
    cmap_custom = mcolors.LinearSegmentedColormap.from_list("custom_gray_green", colors, N=100)

    # Set font to Times New Roman if available
    font_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "times.ttf")
    if os.path.exists(font_path):
        font_manager.fontManager.addfont(font_path)
        # Set global font to Times New Roman
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman']
        plt.rcParams['mathtext.fontset'] = 'stix'
        print(f"Using Times New Roman font from: {font_path}")
    else:
        print("Times New Roman font file not found, using default font")

    # Create a DataFrame with correct names for visualization
    normalized_data_named = normalized_data.copy()
    normalized_data_named.index = available_names

    # Plot clustered heatmap
    g = sns.clustermap(
        normalized_data_named,  # Use named version for visualization
        row_linkage=row_linkage,  # Precomputed hierarchical clustering for rows
        col_cluster=False,  # Disable column clustering (keep original order)
        cmap=cmap_custom,  # Custom color map for heatmap visualization
        vmin=0,  # Minimum value for color scale
        vmax=0.1,  # Maximum value for color scale
        figsize=(fig_size_width, fig_size_height),  # Dynamic figure dimensions
        dendrogram_ratio=0.2,  # Proportion of figure dedicated to dendrogram
        cbar_pos=(0.02, 0.8, 0.03, 0.15),  # Color bar position (left, bottom, width, height)
        yticklabels=False,  # Disable default y-axis labels (we'll add custom ones)
        xticklabels=False,  # Disable default x-axis labels (add custom ones)
        cbar_kws={"label": ""}  # Color bar settings (empty label)
    )

    # Calculate dynamic font sizes
    def calculate_optimal_fontsize(n_items, max_font=10, min_font=4):
        """Calculate optimal font size based on number of items"""
        import math
        if n_items <= 10:
            return max_font
        elif n_items <= 50:
            return max_font - (n_items - 10) * 0.1
        else:
            # For large datasets, use logarithmic scaling
            return max(min_font, max_font - math.log(n_items) * 1.5)

    y_fontsize = calculate_optimal_fontsize(len(available_indices))
    x_fontsize = calculate_optimal_fontsize(len(numeric_cols))
    print(f"Optimal font sizes - y: {y_fontsize:.1f}, x: {x_fontsize:.1f}")

    # Set y-axis labels at the center of each row with correct names
    if len(available_indices) > 0:
        # Get the reordered indices from clustering
        reordered_indices = g.dendrogram_row.reordered_ind
        reordered_names = [available_names[i] for i in reordered_indices]

        # Set ticks at the center of each row (0.5, 1.5, 2.5, ...)
        y_ticks = [i + 0.5 for i in range(len(available_indices))]
        g.ax_heatmap.set_yticks(y_ticks)

        # Set the labels to the reordered names (clustered order)
        g.ax_heatmap.set_yticklabels(reordered_names, fontsize=y_fontsize, va='center')

        print(f"Set {len(reordered_names)} y-axis labels in clustered order")

    # Set x-axis labels at the center of each column
    if len(numeric_cols) > 0:
        # Set ticks at the center of each column (0.5, 1.5, 2.5, ...)
        x_ticks = [i + 0.5 for i in range(len(numeric_cols))]
        g.ax_heatmap.set_xticks(x_ticks)
        g.ax_heatmap.set_xticklabels(numeric_cols, fontsize=x_fontsize, rotation=60, ha='center')

    # Get clustering order and draw enclosing rectangles for FINAL target rows
    reordered_indices = g.dendrogram_row.reordered_ind
    reordered_names = [available_names[i] for i in reordered_indices]

    # Draw rectangles around FINAL target rows (after deduplication)
    for i, row_name in enumerate(reordered_names):
        if row_name in final_target_names:
            # Add enclosing rectangle - now using the center-aligned positions
            rect = patches.Rectangle(
                (0, i),  # Bottom-left coordinates (start of row)
                len(numeric_cols),  # Width = number of columns
                1.0,  # Height = single row
                linewidth=0.8,  # Border line width
                edgecolor='red',  # Border color
                facecolor='none',  # No fill
                linestyle='solid',  # Solid line
                joinstyle='miter',  # Corner style
                zorder=10,  # Stacking order
                clip_on=False  # Not clipped by heatmap edges
            )
            g.ax_heatmap.add_patch(rect)

    # Layout optimization and saving
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    plt.tight_layout()
    heatmap_filename = f"{output_prefix}_{args.correlation_method}_clustering_heatmap.svg"
    plt.savefig(heatmap_filename, bbox_inches="tight", dpi=300)
    plt.show()
    print(f"Heatmap saved as: {heatmap_filename}")

    # Save results
    print("\nSaving processed data...")
    try:
        filtered_output = [(row, score, p_value)
                           for row, score, p_value in filtered_rows
                           if row not in final_target_indices]

        # Save filtered results (NON-BENCHMARKS ONLY)
        if args.correlation_method == 'pearson':
            df_filtered = pd.DataFrame(
                filtered_output,
                columns=["Contig_Name", "Mean_Pearson_Correlation", "P_Value"]
            )
        else:
            df_filtered = pd.DataFrame(
                filtered_output,
                columns=["Contig_Name", "Mean_Spearman_Correlation", "P_Value"]
            )

        df_filtered.to_csv(f"{output_prefix}_filtered_contigs.csv", index=False)
        print(f"Filtered results (non-benchmarks only) saved as '{output_prefix}_filtered_contigs.csv'")

        # Save normalized data with original names
        normalized_data_with_names = normalized_data.copy()
        normalized_data_with_names.index = available_names
        normalized_data_with_names.to_csv(f'{output_prefix}_normalized_feature_vectors.csv')
        print(f"Normalized feature data saved as '{output_prefix}_normalized_feature_vectors.csv'")

        # Save distance matrix with original names
        distance_df = pd.DataFrame(distance_matrix, index=available_names, columns=available_names)
        distance_df.to_csv(f'{output_prefix}_{args.correlation_method}_distance_matrix.csv')
        print(
            f"{args.correlation_method.capitalize()} distance matrix saved as '{output_prefix}_{args.correlation_method}_distance_matrix.csv'")

        # Save clustering order
        reordered_indices = g.dendrogram_row.reordered_ind
        reordered_names = [available_names[i] for i in reordered_indices]

        clustering_info = pd.DataFrame({
            'Original_Name': available_names,
            'Cluster_Order': [list(reordered_names).index(name) if name in reordered_names else -1
                              for name in available_names]
        })
        clustering_info.to_csv(f'{output_prefix}_clustering_order.csv', index=False)
        print(f"Clustering order saved as '{output_prefix}_clustering_order.csv'")

        # Save detailed information
        with open(f'{output_prefix}_analysis_summary.txt', 'w') as f:
            f.write(f"{vector_mode} {args.correlation_method.capitalize()} Correlation Analysis Summary\n")
            f.write("=" * 50 + "\n")
            f.write(f"Input file: {args.input}\n")
            f.write(f"Detected vector mode: {vector_mode}\n")
            f.write(f"Correlation method: {args.correlation_method}\n")
            f.write(f"Length range: {args.min_length}-{args.max_length} nt\n")
            f.write(f"Actual dimension: {len(numeric_cols)}\n")
            f.write(f"Benchmark contigs: {len(target_rows)}\n")
            f.write(f"Valid benchmark contigs: {len(final_target_indices)}\n")
            f.write(f"Mean correlation threshold: {args.mean_r}\n")
            f.write(f"P-value threshold: {args.p_value}\n")
            f.write(f"Heatmap data: Normalized feature values (0-1 proportion)\n")
            f.write(
                f"Clustering method: Hierarchical clustering based on {args.correlation_method} correlation (distance = 1 - correlation)\n")
            f.write(f"Data normalization: Each feature value divided by row total (proportion normalization)\n")

            # Feature description based on mode
            length_count = args.max_length - args.min_length + 1
            if vector_mode == 'size':
                f.write(
                    f"Feature description: {length_count}-dimensional length distribution ({args.min_length}-{args.max_length}nt)\n")
            elif vector_mode == 'size_P_5nt':
                f.write(
                    f"Feature description: {length_count + 4}-dimensional (length distribution + A/T/C/G base composition)\n")
            elif vector_mode == 'sizeXstr':
                f.write(f"Feature description: {length_count * 2}-dimensional (length × strand specificity)\n")
            elif vector_mode == 'sizeX5nt':
                f.write(f"Feature description: {length_count * 4}-dimensional (length × A/T/C/G base specificity)\n")
            elif vector_mode == 'sizeX5ntXstr':
                f.write(
                    f"Feature description: {length_count * 8}-dimensional (length × A/T/C/G × strand specificity)\n")

            f.write(f"Contigs after threshold filtering: {len(filtered_rows)}\n")
            f.write(f"Non-benchmark contigs for output: {len(filtered_output)}\n")
            f.write(f"Heatmap data dimension: {normalized_data.shape[0]} rows × {normalized_data.shape[1]} columns\n")
            f.write(f"Feature value range: {normalized_data.values.min():.6f} to {normalized_data.values.max():.6f}\n")
            f.write(f"Distance matrix range: {distance_matrix.min():.3f} to {distance_matrix.max():.3f}\n")
            f.write("\nBenchmark contigs (after strand selection, marked with red borders):\n")
            for name in final_target_names:
                f.write(f"  - {name}\n")
            f.write("\nFiltered highly correlated contigs (non-benchmarks, top 10):\n")
            non_benchmark_top10 = filtered_output[:10]
            for row, score, p_value in non_benchmark_top10:
                f.write(f"  - {row} (mean correlation: {score:.4f}, p-value: {p_value:.4e})\n")

        print("All data files saved successfully:")
        print(f"  - {output_prefix}_filtered_contigs.csv: Filtered results (non-benchmarks only)")
        print(f"  - {output_prefix}_normalized_feature_vectors.csv: Normalized feature data")
        print(
            f"  - {output_prefix}_{args.correlation_method}_distance_matrix.csv: {args.correlation_method.capitalize()} distance matrix")
        print(f"  - {output_prefix}_clustering_order.csv: Clustering order")
        print(f"  - {output_prefix}_analysis_summary.txt: Analysis summary")
        print(f"  - {output_prefix}_{args.correlation_method}_clustering_heatmap.svg: Clustering heatmap")

    except Exception as e:
        print(f"Error saving files: {e}")

    print("\n" + "=" * 60)
    print(f"{vector_mode} {args.correlation_method.capitalize()} Analysis Completion Summary:")
    print(f"Input file: {args.input}")
    print(f"Detected vector mode: {vector_mode}")
    print(f"Correlation method: {args.correlation_method}")
    print(f"Length range: {args.min_length}-{args.max_length} nt")
    print(f"Actual dimension: {len(numeric_cols)}")
    print(f"Benchmark contigs: {len(final_target_indices)}")
    print(f"Filtered contigs: {len(filtered_rows)}")
    print(f"Non-benchmark contigs for output: {len(filtered_output)}")
    print(f"Heatmap data dimension: {normalized_data.shape[0]} rows × {normalized_data.shape[1]} columns")
    print("=" * 60)

    # Output detailed information of filtered results (NON-BENCHMARKS ONLY)
    if filtered_output:
        print(f"\nTop 10 filtered results (NON-BENCHMARKS ONLY, sorted by mean correlation):")
        for i, (row, score, p_value) in enumerate(filtered_output[:10], 1):
            print(f"{i}. {row} (mean correlation: {score:.4f}, p-value: {p_value:.4e})")


if __name__ == "__main__":
    main()