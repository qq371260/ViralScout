#!/usr/bin/env python3
"""
Spearman correlation of sRNAs feature vectors of viral and other contigs
"""

import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from collections import defaultdict
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage
import warnings
import matplotlib as mpl
import argparse
import os
from scipy.spatial.distance import squareform
import re

warnings.filterwarnings('ignore')
import matplotlib.colors as mcolors
from scipy.stats import combine_pvalues


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Spearman correlation-based feature vector clustering analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Required parameters
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Input CSV file path containing feature vector data')

    # Benchmark parameters (optional)
    parser.add_argument('-b', '--benchmarks', type=str, default=None,
                        help='Benchmark contigs file path (txt format, one contig name per line). If not provided, use default vsiRNA simulant reference')

    # Other parameters
    parser.add_argument('-s', '--specific_contigs', type=int, default=1000,
                        help='Top K highly correlated contigs for each benchmark (default: 1000)')
    parser.add_argument('-c', '--common_contigs', type=int, default=999,
                        help='Final number of common top L contigs (default: 999)')
    parser.add_argument('-n', '--number_21_nt', type=int, default=10,
                        help='21nt four-column sum threshold (default: 10)')
    parser.add_argument('-m', '--mean_r', type=float, default=0.8,
                        help='Composite score threshold (default: 0.8)')
    parser.add_argument('-p', '--p_value', type=float, default=0.05,
                        help='P-value threshold (default: 0.05)')
    parser.add_argument('-o', '--output', type=str, default='',
                        help='Output file prefix (default: use vector mode as prefix)')
    parser.add_argument('-t', '--threads', type=int, default=-1,
                        help='Number of parallel threads (default: -1, use all CPU cores)')
    parser.add_argument('--cell_height', type=float, default=0.1,
                        help='Heatmap cell height (default: 0.1)')
    parser.add_argument('--cell_width', type=float, default=0.1,
                        help='Heatmap cell width (default: 0.1)')
    parser.add_argument('--no_plot', action='store_true',
                        help='Do not generate heatmap')
    parser.add_argument('--min_length', type=int, default=18,
                        help='Minimum read length for dimension calculation (default: 18)')
    parser.add_argument('--max_length', type=int, default=30,
                        help='Maximum read length for dimension calculation (default: 30)')

    return parser.parse_args()


def load_benchmarks(benchmarks_arg, vector_mode=None):
    """
    Load benchmark contigs
    Support file path, direct contig names, or use default reference
    vector_mode: detected vector mode for dynamic filename
    """

    def extract_first_word(name):
        first_word = re.split(r'[\s\t,;]', name.strip())[0]
        return first_word

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


def main():
    # Parse arguments
    args = parse_arguments()

    # Parameter validation
    assert args.specific_contigs > args.common_contigs > 1, "Must satisfy: specific_contigs > common_contigs > 1"

    # Reading and detect vector mode
    print("Reading feature vector data...")
    df = pd.read_csv(args.input)
    df = df.set_index(df.columns[0])

    # Get all numeric columns
    numeric_cols = [col for col in df.select_dtypes(include=[np.number]).columns.tolist() if col != '21nt_total']

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
    print("Spearman Correlation Clustering Analysis")
    print("=" * 60)
    print(f"Input file: {args.input}")
    print(f"Detected vector mode: {vector_mode}")
    print(f"Length range: {args.min_length}-{args.max_length} nt")
    print(f"Actual dimension: {len(numeric_cols)}")
    print(f"Benchmark contigs: {len(target_rows)}")
    print(f"Specific contigs (K): {args.specific_contigs}")
    print(f"Common contigs (L): {args.common_contigs}")
    print(f"21nt threshold: {args.number_21_nt}")
    print(f"Mean correlation threshold: {args.mean_r}")
    print(f"P-value threshold: {args.p_value}")
    print(f"Output prefix: {output_prefix}")
    print(f"Threads: {args.threads}")
    print("=" * 60)

    # Read vector data (if not already read)
    if 'df' not in locals():
        df = pd.read_csv(args.input)
        df = df.set_index(df.columns[0])

    # Check and extract 21nt related columns
    twentyone_nt_cols = [col for col in df.columns if '21' in col]
    if len(twentyone_nt_cols) != 4:
        print(f"Warning: Found {len(twentyone_nt_cols)} 21nt-related columns, expected 4")
        print(f"Found columns: {twentyone_nt_cols}")
    else:
        print(f"Found 21nt-related columns: {twentyone_nt_cols}")

    # Calculate 21nt four-column sum and add to dataframe
    df['21nt_total'] = df[twentyone_nt_cols].sum(axis=1)
    print(
        f"21nt total statistics: min={df['21nt_total'].min()}, max={df['21nt_total'].max()}, mean={df['21nt_total'].mean():.2f}")

    # Filter based on 21nt total
    initial_count = len(df)
    df = df[df['21nt_total'] >= args.number_21_nt]
    filtered_count = len(df)
    print(f"Filtered by 21nt total threshold {args.number_21_nt}: {initial_count} -> {filtered_count} rows")

    # Validate target row availability
    valid_target_rows = [row for row in target_rows if row in df.index]
    print(f"Valid benchmark contigs: {len(valid_target_rows)}/{len(target_rows)}")
    if len(valid_target_rows) < len(target_rows):
        missing = set(target_rows) - set(valid_target_rows)
        print(f"Missing benchmark contigs: {missing}")

    # Get all numeric columns
    numeric_cols = [col for col in df.select_dtypes(include=[np.number]).columns.tolist() if col != '21nt_total']

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

    def calculate_correlation_and_pvalue(target_row, other_row, df, numeric_cols):
        """Calculate Spearman correlation coefficient and p-value between two rows"""
        try:
            y = df.loc[target_row, numeric_cols].values.astype(float)
            x = df.loc[other_row, numeric_cols].values.astype(float)

            valid_mask = ~(np.isnan(x) | np.isnan(y))
            x_valid = x[valid_mask]
            y_valid = y[valid_mask]

            if len(x_valid) <= 3:
                return 0.0, 1.0

            corr, p_value = spearmanr(x_valid, y_valid)

            if np.isnan(corr):
                return 0.0, 1.0

            return corr, p_value

        except Exception as e:
            print(f"Error calculating correlation between {target_row} and {other_row}: {e}")
            return 0.0, 1.0

    def get_top_k_with_p_value(target_row, k, df, numeric_cols):
        """Get top k rows with highest correlation to target row and their p-values"""
        correlations = {}
        p_values = {}

        for other_row in df.index:
            if other_row == target_row:
                continue

            corr, p_value = calculate_correlation_and_pvalue(target_row, other_row, df, numeric_cols)
            if corr > 0:
                correlations[other_row] = corr
                p_values[other_row] = p_value

        sorted_rows = sorted(correlations.items(), key=lambda x: x[1], reverse=True)[:k]
        return [(row, correlations[row], p_values[row]) for row, _ in sorted_rows]

    # Parallel processing for each target row
    print("Calculating Spearman correlation for feature vectors...")
    top_k_with_p_values = Parallel(n_jobs=args.threads)(
        delayed(get_top_k_with_p_value)(target_row, args.specific_contigs, df, numeric_cols)
        for target_row in valid_target_rows
    )

    # Organize results as dictionary
    row_rankings_with_p_values = {
        target_row: {row: (corr, p_value) for row, corr, p_value in items}
        for target_row, items in zip(valid_target_rows, top_k_with_p_values)
    }

    # Fix: Check if we have valid results before calculating intersection
    if not row_rankings_with_p_values:
        print("Error: No valid correlation results found")
        return

    # Fix: Proper set intersection calculation
    common_rows = set.intersection(*[set(d.keys()) for d in row_rankings_with_p_values.values()])
    print(f"Found {len(common_rows)} common highly correlated contigs")

    # Calculate composite score and combined p-value
    scoring = defaultdict(lambda: {'score': 0.0, 'p_value': 1.0})
    for row in common_rows:
        p_values = []
        for target_row in valid_target_rows:
            corr, p_value = row_rankings_with_p_values[target_row][row]
            scoring[row]['score'] += corr
            p_values.append(p_value)

        scoring[row]['score'] /= len(valid_target_rows)

        if p_values:
            clean_p_values = [max(min(p, 0.999), 0.001) for p in p_values]
            try:
                _, combined_p = combine_pvalues(clean_p_values, method='fisher')
                scoring[row]['p_value'] = combined_p
            except:
                scoring[row]['p_value'] = min(1.0, len(p_values) * min(p_values))

    # Take top L results, sorted by composite score in descending order
    final_rows = sorted(scoring.items(), key=lambda x: x[1]['score'], reverse=True)[:args.common_contigs]

    # Apply threshold filtering
    filtered_rows = [(row, values['score'], values['p_value'])
                     for row, values in final_rows
                     if values['score'] >= args.mean_r and values['p_value'] < args.p_value]

    print(
        f"After threshold filtering (mean correlation>={args.mean_r}, p-value<{args.p_value}), retained {len(filtered_rows)} contigs")

    # Extract filtered row names
    filtered_row_names = [row for row, score, p_value in filtered_rows]

    # Combine target rows and filtered rows for clustering
    combined_row_names = filtered_row_names + valid_target_rows
    combined_row_names = list(set(combined_row_names))

    # Ensure all rows exist in dataframe
    available_rows = [row for row in combined_row_names if row in df.index]
    print(f"Total rows for heatmap analysis: {len(available_rows)}")

    # Strict Spearman correlation-based clustering
    print("Starting strict Spearman correlation-based clustering analysis...")

    # Extract feature data for filtered rows
    feature_data = df.loc[available_rows, numeric_cols].copy()

    # Normalize each row (value/total, i.e., calculate proportion)
    print("Normalizing row data (calculating proportions)...")
    row_sums = feature_data.sum(axis=1)
    normalized_data = feature_data.div(row_sums, axis=0)

    print(f"Normalized data shape: {normalized_data.shape}")
    print(f"Normalized data range: min={normalized_data.values.min():.6f}, max={normalized_data.values.max():.6f}")

    # Calculate Spearman distance matrix (1 - Spearman correlation)
    print("Calculating Spearman distance matrix...")
    n_samples = len(available_rows)
    distance_matrix = np.zeros((n_samples, n_samples))

    for i in range(n_samples):
        for j in range(i, n_samples):
            if i == j:
                distance_matrix[i, j] = 0.0
            else:
                row1 = available_rows[i]
                row2 = available_rows[j]
                data1 = normalized_data.loc[row1].values
                data2 = normalized_data.loc[row2].values

                # Create valid data mask
                valid_mask = ~(np.isnan(data1) | np.isnan(data2))
                data1_valid = data1[valid_mask]
                data2_valid = data2[valid_mask]

                if len(data1_valid) > 3:
                    corr, _ = spearmanr(data1_valid, data2_valid)
                    # Convert correlation to distance: distance = 1 - correlation
                    distance = 1 - corr
                    if np.isnan(distance):
                        distance = 2.0
                else:
                    distance = 2.0  # Maximum distance

                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance

    print(f"Distance matrix range: {distance_matrix.min():.3f} to {distance_matrix.max():.3f}")

    # Hierarchical clustering using Spearman distance matrix
    print("Performing hierarchical clustering...")

    # distance_matrix is n x n symmetric
    condensed_dist = squareform(distance_matrix, checks=False)
    row_linkage = linkage(condensed_dist, method='average')

    # Visualization section
    if not args.no_plot:
        print("Plotting strictly corresponding heatmap and dendrogram using clustermap...")

        import matplotlib.patches as patches

        # Custom figure size (automatically adjusted based on number of rows and columns)
        min_width = 10
        min_height = 8
        fig_size_width = max(min_width, len(numeric_cols) * args.cell_width)
        fig_size_height = max(min_height, len(available_rows) * args.cell_height)
        print(f"Setting figure size: {fig_size_width:.1f} × {fig_size_height:.1f} inches")

        # Custom color map
        colors = ['#F0F0F0', '#E0F0E0', '#C8E6C9', '#A5D6A7', '#81C784',
                  '#66BB6A', '#4CAF50', '#388E3C', '#2E7D32', '#1B5E20']
        cmap_custom = mcolors.LinearSegmentedColormap.from_list("custom_gray_green", colors, N=100)

        # Plot clustered heatmap
        g = sns.clustermap(
            normalized_data,
            row_linkage=row_linkage,
            col_cluster=False,
            cmap=cmap_custom,
            vmin=0,
            vmax=0.1,
            figsize=(fig_size_width, fig_size_height),
            dendrogram_ratio=0.15,
            cbar_pos=(0.02, 0.8, 0.03, 0.15),
            yticklabels=False,
            xticklabels=False,
            cbar_kws={"label": ""}
        )

        # Get clustering order and draw enclosing rectangles for target rows
        reordered_indices = g.dendrogram_row.reordered_ind
        reordered_rows = [available_rows[i] for i in reordered_indices]

        for i, row_name in enumerate(reordered_rows):
            if row_name in valid_target_rows:
                # Add enclosing rectangle (adjusted to ensure all sides are fully displayed)
                rect = patches.Rectangle(
                    (-0.01, i - 0.01),  # Bottom-left coordinates (slightly offset inward)
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
        heatmap_filename = f"{output_prefix}_spearman_clustering_heatmap.svg"
        plt.savefig(heatmap_filename, bbox_inches="tight", dpi=300)
        plt.show()
        print(f"Heatmap saved as: {heatmap_filename}")

    # Save results
    print("\nSaving processed data...")
    try:
        # Save filtered results
        df_filtered = pd.DataFrame(filtered_rows, columns=["Contig_Name", "Mean_Spearman_Correlation", "P_Value"])
        df_filtered.to_csv(f"{output_prefix}_filtered_contigs.csv", index=False)
        print(f"Filtered results saved as '{output_prefix}_filtered_contigs.csv'")

        # Save normalized data
        normalized_data.to_csv(f'{output_prefix}_normalized_features.csv')
        print(f"Normalized feature data saved as '{output_prefix}_normalized_features.csv'")

        # Save distance matrix
        distance_df = pd.DataFrame(distance_matrix, index=available_rows, columns=available_rows)
        distance_df.to_csv(f'{output_prefix}_spearman_distance_matrix.csv')
        print(f"Spearman distance matrix saved as '{output_prefix}_spearman_distance_matrix.csv'")

        # Save clustering order
        clustering_info = pd.DataFrame({
            'Original_Index': range(len(available_rows)),
            'Contig_Name': available_rows,
            'Cluster_Order': [list(reordered_rows).index(row) if row in reordered_rows else -1
                              for row in available_rows]
        })
        clustering_info.to_csv(f'{output_prefix}_clustering_order.csv', index=False)
        print(f"Clustering order saved as '{output_prefix}_clustering_order.csv'")

        # Save detailed information
        with open(f'{output_prefix}_analysis_summary.txt', 'w') as f:
            f.write(f"{vector_mode} Spearman Correlation Analysis Summary\n")
            f.write("=" * 50 + "\n")
            f.write(f"Input file: {args.input}\n")
            f.write(f"Detected vector mode: {vector_mode}\n")
            f.write(f"Length range: {args.min_length}-{args.max_length} nt\n")
            f.write(f"Actual dimension: {len(numeric_cols)}\n")
            f.write(f"Benchmark contigs: {len(target_rows)}\n")
            f.write(f"Valid benchmark contigs: {len(valid_target_rows)}\n")
            f.write(f"Specific contigs (K): {args.specific_contigs}\n")
            f.write(f"Common contigs (L): {args.common_contigs}\n")
            f.write(f"21nt threshold: {args.number_21_nt}\n")
            f.write(f"Mean correlation threshold: {args.mean_r}\n")
            f.write(f"P-value threshold: {args.p_value}\n")
            f.write(f"Heatmap data: Normalized feature values (0-1 proportion)\n")
            f.write(
                f"Clustering method: Hierarchical clustering based on Spearman correlation (distance = 1 - Spearman correlation)\n")
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

            f.write(f"21nt total filter threshold: {args.number_21_nt}\n")
            f.write(f"Common highly correlated contigs: {len(common_rows)}\n")
            f.write(f"Contigs after threshold filtering: {len(filtered_rows)}\n")
            f.write(f"Heatmap data dimension: {normalized_data.shape[0]} rows × {normalized_data.shape[1]} columns\n")
            f.write(f"Feature value range: {normalized_data.values.min():.6f} to {normalized_data.values.max():.6f}\n")
            f.write(f"Distance matrix range: {distance_matrix.min():.3f} to {distance_matrix.max():.3f}\n")
            f.write("\nBenchmark contigs (marked with red borders):\n")
            for row in valid_target_rows:
                f.write(f"  - {row}\n")
            f.write("\nFiltered highly correlated contigs (top 10):\n")
            for row, score, p_value in filtered_rows[:10]:
                f.write(f"  - {row} (mean correlation: {score:.4f}, p-value: {p_value:.4e})\n")

        print("All data files saved successfully:")
        print(f"  - {output_prefix}_filtered_contigs.csv: Filtered results")
        print(f"  - {output_prefix}_normalized_features.csv: Normalized feature data")
        print(f"  - {output_prefix}_spearman_distance_matrix.csv: Spearman distance matrix")
        print(f"  - {output_prefix}_clustering_order.csv: Clustering order")
        print(f"  - {output_prefix}_analysis_summary.txt: Analysis summary")
        if not args.no_plot:
            print(f"  - {output_prefix}_spearman_clustering_heatmap.svg: Clustering heatmap")

    except Exception as e:
        print(f"Error saving files: {e}")

    print("\n" + "=" * 60)
    print(f"{vector_mode} Spearman Clustering Analysis Completion Summary:")
    print(f"Input file: {args.input}")
    print(f"Detected vector mode: {vector_mode}")
    print(f"Length range: {args.min_length}-{args.max_length} nt")
    print(f"Actual dimension: {len(numeric_cols)}")
    print(f"Benchmark contigs: {len(valid_target_rows)}")
    print(f"Filtered contigs: {len(filtered_rows)}")
    print(f"Heatmap data dimension: {normalized_data.shape[0]} rows × {normalized_data.shape[1]} columns")
    print("=" * 60)

    # Output detailed information of filtered results
    if filtered_rows:
        print(f"\nTop 10 filtered results (sorted by mean correlation in descending order):")
        for i, (row, score, p_value) in enumerate(filtered_rows[:10], 1):
            print(f"{i}. {row} (mean correlation: {score:.4f}, p-value: {p_value:.4e})")


if __name__ == "__main__":
    main()