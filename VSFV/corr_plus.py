#!/usr/bin/env python3
"""
Enhanced correlation analysis with additional genomic features
Adds AC_normalized, GC_normalized, and 21-22_nt_normalized to feature vectors
with length-based column merging
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
    parser = argparse.ArgumentParser(
        description='Enhanced correlation analysis with genomic features and length-based column merging',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='Input CSV file path containing normalized feature vector data for ALL contigs')
    parser.add_argument('-f', '--features', type=str, required=True,
                        help='CSV file with genomic features for NON-BENCHMARK contigs (Avg_Coverage, GC_Content, sRNA percentages)')
    parser.add_argument('-b', '--benchmark_features', type=str, required=True,
                        help='CSV file with genomic features for BENCHMARK contigs (Avg_Coverage, GC_Content, sRNA percentages)')
    parser.add_argument('-o', '--output', type=str, default='enhanced_analysis',
                        help='Output file prefix')
    parser.add_argument('-d', '--output_dir', type=str, default='enhanced_corr_results',
                        help='Output directory for all files')
    parser.add_argument('-t', '--threads', type=int, default=max(1, (os.cpu_count() or 1) - 1),
                        help='Number of parallel threads')
    parser.add_argument('--min_length', type=int, default=18,
                        help='Minimum sRNA size for column merging')
    parser.add_argument('--max_length', type=int, default=30,
                        help='Maximum sRNA size for column merging')
    return parser.parse_args()


def ensure_directory(directory):
    """Create directory if it doesn't exist"""
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Created output directory: {directory}")
    return directory


def extract_length_from_column(column_name):
    """
    Extract length number from column name.
    Returns None if no length number is found.
    """
    # Look for numbers in the column name
    matches = re.findall(r'\b(\d+)\b', str(column_name))
    if matches:
        return int(matches[0])
    return None


def merge_length_columns(df_norm, min_length, max_length):
    """
    Merge columns with the same length number within specified range.
    Only processes columns that contain length numbers.
    """
    print(f"\nMerging columns for lengths {min_length} to {max_length}...")

    # Create a dictionary to group columns by length
    length_groups = {}

    # Process each column
    for col in df_norm.columns:
        length = extract_length_from_column(col)
        if length is not None and min_length <= length <= max_length:
            if length not in length_groups:
                length_groups[length] = []
            length_groups[length].append(col)

    print(f"Found {len(length_groups)} unique lengths in range")

    # Create new merged dataframe
    merged_data = {}

    # For each length, sum all corresponding columns
    for length in sorted(length_groups.keys()):
        cols = length_groups[length]
        if len(cols) > 1:
            print(f"  Merging {len(cols)} columns for length {length}: {', '.join(cols[:3])}...")
        merged_data[f"{length}-nt"] = df_norm[cols].sum(axis=1)

    # Create new dataframe with merged columns
    df_merged = pd.DataFrame(merged_data, index=df_norm.index)

    print(f"Original columns: {len(df_norm.columns)}")
    print(f"Merged columns: {len(df_merged.columns)}")
    print(f"Final dimension: {len(df_merged.columns)} length columns")

    return df_merged


def load_and_enhance_features(normalized_file, non_benchmark_features_file, benchmark_features_file, min_length,
                              max_length):
    """
    Load normalized feature vectors, merge length-based columns, and enhance with genomic features
    """
    print("Loading and enhancing feature vectors...")

    # Load normalized feature vectors (all contigs)
    df_norm = pd.read_csv(normalized_file)
    df_norm = df_norm.set_index(df_norm.columns[0])
    print(f"Loaded normalized features for {len(df_norm)} contigs")
    print(f"Original feature columns: {len(df_norm.columns)}")

    # STEP 1: Merge columns by length
    df_merged = merge_length_columns(df_norm, min_length, max_length)

    # Load genomic features for non-benchmark contigs
    df_non_benchmark_features = pd.read_csv(non_benchmark_features_file)
    print(f"Loaded genomic features for {len(df_non_benchmark_features)} non-benchmark contigs")

    # Load genomic features for benchmark contigs
    df_benchmark_features = pd.read_csv(benchmark_features_file)
    print(f"Loaded genomic features for {len(df_benchmark_features)} benchmark contigs")

    # Extract benchmark contig names
    benchmark_contigs = df_benchmark_features['Contig_Name'].tolist()
    print(f"Identified {len(benchmark_contigs)} benchmark contigs")

    # Simple and robust name matching solution
    def remove_minus_suffix(name):
        """Remove _minus suffix from contig names"""
        return name[:-6] if name.endswith('_minus') else name

    # Create name mapping for vector file (base name -> all possible original names)
    vector_name_map = {}
    for orig_name in df_merged.index:
        base_name = remove_minus_suffix(orig_name)
        if base_name not in vector_name_map:
            vector_name_map[base_name] = []
        vector_name_map[base_name].append(orig_name)

    # Remap benchmark contigs to actual names in vector file
    mapped_benchmarks = []
    for bench_name in benchmark_contigs:
        base_bench = remove_minus_suffix(bench_name)
        if base_bench in vector_name_map:
            # Prefer exact match with benchmark name, otherwise use first available
            if bench_name in vector_name_map[base_bench]:
                selected_name = bench_name
            else:
                selected_name = vector_name_map[base_bench][0]
            mapped_benchmarks.append(selected_name)
            if bench_name != selected_name:
                print(f"  Benchmark mapped: '{bench_name}' -> '{selected_name}'")
        else:
            print(f"Warning: Benchmark '{bench_name}' not found in vector file")

    print(f"Successfully mapped {len(mapped_benchmarks)}/{len(benchmark_contigs)} benchmark contigs")

    # Update benchmark contigs list with mapped names
    benchmark_contigs = mapped_benchmarks

    # Combine all genomic features
    df_features = pd.concat([df_non_benchmark_features, df_benchmark_features], ignore_index=True)

    # Create name mapping for genomic features
    feature_name_map = {}
    for orig_name in df_features['Contig_Name']:
        base_name = remove_minus_suffix(orig_name)
        feature_name_map[base_name] = orig_name  # Store original name

    # STEP 2: Calculate total dimensions
    length_dimensions = max_length - min_length + 1
    total_dimensions = length_dimensions + 3

    print(f"\nDimension calculation:")
    print(f"  Length columns: {length_dimensions} ({min_length}-{max_length} nt)")
    print(f"  Extra features: 3")
    print(f"  Total dimensions: {total_dimensions}")

    # Create enhanced feature dataframe
    enhanced_features = []

    for contig in df_merged.index:
        # Get existing merged length features
        length_features = df_merged.loc[contig].to_dict()

        # Find genomic features using base name matching
        base_contig = remove_minus_suffix(contig)
        if base_contig in feature_name_map:
            feature_contig_name = feature_name_map[base_contig]
            genomic_data = df_features[df_features['Contig_Name'] == feature_contig_name]

            if len(genomic_data) > 0:
                row = genomic_data.iloc[0]
                try:
                    AC_normalized = np.log10(float(row['Avg_Coverage']) + 1) / 10
                    GC_normalized = float(row['GC_Content_percent']) / 100
                    sRNA_21_22_nt_normalized = float(row['sRNA_21_22nt_Percent']) / 100

                    # Create enhanced feature vector
                    enhanced_row = {
                        'Contig_Name': contig,  # Keep original name from vector file
                        'AC_normalized': AC_normalized,
                        'GC_content': GC_normalized,
                        '21-22-nt': sRNA_21_22_nt_normalized
                    }

                    # Add all merged length features
                    for key, value in length_features.items():
                        enhanced_row[key] = value

                    enhanced_features.append(enhanced_row)

                except (ValueError, KeyError) as e:
                    print(f"Warning: Could not process genomic features for {contig}: {e}")
            else:
                print(f"Warning: No genomic features found for {contig} (mapped from '{feature_contig_name}')")
        else:
            print(f"Warning: No genomic features mapping found for {contig}")

    # Create enhanced dataframe
    df_enhanced = pd.DataFrame(enhanced_features)
    df_enhanced = df_enhanced.set_index('Contig_Name')

    print(f"\nEnhanced feature matrix: {df_enhanced.shape[0]} contigs × {df_enhanced.shape[1]} features")
    print(f"Length features: {length_dimensions} columns ({min_length}-nt to {max_length}-nt)")
    print(f"Additional features: AC_normalized, GC_content, 21-22-nt")
    print(f"Benchmark contigs: {len(benchmark_contigs)}")

    # Reorder columns to have additional features first, then length columns
    additional_cols = ['AC_normalized', 'GC_content', '21-22-nt']
    length_cols = [col for col in df_enhanced.columns if col not in additional_cols]
    df_enhanced = df_enhanced[additional_cols + sorted(length_cols, key=lambda x: int(x.split('-')[0]))]

    return df_enhanced, benchmark_contigs


def calculate_correlation_and_pvalue(target_row, other_row, df, numeric_cols):
    """Calculate pearson correlation coefficient and p-value between two rows"""
    try:
        y = df.loc[target_row, numeric_cols].values.astype(float)
        x = df.loc[other_row, numeric_cols].values.astype(float)

        valid_mask = ~(np.isnan(x) | np.isnan(y))
        x_valid = x[valid_mask]
        y_valid = y[valid_mask]

        if len(x_valid) <= 3:
            return 0.0, 1.0

        corr, p_value = pearsonr(x_valid, y_valid)

        if np.isnan(corr):
            return 0.0, 1.0

        return corr, p_value

    except Exception as e:
        print(f"Error calculating correlation between {target_row} and {other_row}: {e}")
        return 0.0, 1.0


def calculate_all_correlations(target_row, df, numeric_cols):
    """Calculate correlations between target row and all other rows"""
    correlations = {}
    p_values = {}

    for other_row in df.index:
        if other_row == target_row:
            continue

        corr, p_value = calculate_correlation_and_pvalue(target_row, other_row, df, numeric_cols)
        correlations[other_row] = corr
        p_values[other_row] = p_value

    return correlations, p_values


def create_enhanced_heatmap(df_enhanced, benchmark_contigs, output_prefix, output_dir, min_length, max_length):
    """Create enhanced clustering heatmap"""
    print("Creating enhanced clustering heatmap...")

    # Get all numeric columns
    numeric_cols = [col for col in df_enhanced.select_dtypes(include=[np.number]).columns.tolist()]
    normalized_data = df_enhanced[numeric_cols].copy()

    print(f"Normalized data shape: {normalized_data.shape}")
    print(f"Normalized data range: min={normalized_data.values.min():.6f}, max={normalized_data.values.max():.6f}")

    # Calculate distance matrix (1 - correlation)
    print(f"Calculating pearson distance matrix...")
    n_samples = len(df_enhanced)
    distance_matrix = np.zeros((n_samples, n_samples))
    available_indices = df_enhanced.index.tolist()

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
                    corr, _ = pearsonr(data1_valid, data2_valid)
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

    # Visualization
    print("Plotting heatmap and dendrogram using clustermap...")

    # Figure size
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
    normalized_data_named.index = available_indices

    # Plot clustered heatmap
    g = sns.clustermap(
        normalized_data_named,  # Use named version for visualization
        row_linkage=row_linkage,  # Precomputed hierarchical clustering for rows
        col_cluster=False,  # Disable column clustering (keep original order)
        cmap=cmap_custom,  # Custom color map for heatmap visualization
        vmin=0,  # Minimum value for color scale
        vmax=0.5,  # Maximum value for color scale
        figsize=(fig_size_width, fig_size_height),  # Dynamic figure dimensions
        dendrogram_ratio=0.2,  # Proportion of figure dedicated to dendrogram
        cbar_pos=(0.02, 0.8, 0.03, 0.15),  # Color bar position (left, bottom, width, height)
        yticklabels=False,  # Disable default y-axis labels (we'll add custom ones)
        xticklabels=False,  # Disable default x-axis labels (add custom ones)
        cbar_kws={"label": ""}  # Color bar settings (empty label)
    )

    # Calculate dynamic font sizes
    def calculate_optimal_fontsize(n_items, max_font=14, min_font=4):
        """Calculate optimal font size based on number of items"""
        import math
        if n_items <= 10:
            return max_font
        elif n_items <= 50:
            return max_font - (n_items - 10) * 0.1
        else:
            # For large datasets, use logarithmic scaling
            return max(min_font, max_font - math.log(n_items) * 2)

    y_fontsize = calculate_optimal_fontsize(len(available_indices))
    x_fontsize = calculate_optimal_fontsize(len(numeric_cols))
    print(f"Optimal font sizes - y: {y_fontsize:.1f}, x: {x_fontsize:.1f}")

    # Set y-axis labels at the center of each row with correct names
    if len(available_indices) > 0:
        # Get the reordered indices from clustering
        reordered_indices = g.dendrogram_row.reordered_ind
        reordered_names = [available_indices[i] for i in reordered_indices]

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
        g.ax_heatmap.set_xticklabels(numeric_cols, fontsize=x_fontsize, rotation=90, ha='center')

    # Get clustering order and draw enclosing rectangles for FINAL target rows
    reordered_indices = g.dendrogram_row.reordered_ind
    reordered_names = [available_indices[i] for i in reordered_indices]

    # Draw rectangles around FINAL target rows (after deduplication)
    for i, row_name in enumerate(reordered_names):
        if row_name in benchmark_contigs:
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

    # Save heatmap to output directory
    heatmap_filename = os.path.join(output_dir, f"{output_prefix}_enhanced_pearson_clustering_heatmap.svg")
    plt.savefig(heatmap_filename, bbox_inches="tight", dpi=300)
    print(f"Enhanced heatmap saved as: {heatmap_filename}")

    return normalized_data, distance_matrix, reordered_names


def main():
    # Parse arguments
    args = parse_arguments()

    # Create output directory
    output_dir = ensure_directory(args.output_dir)
    print(f"All output files will be saved to: {output_dir}")

    # Load and enhance features with length-based column merging
    df_enhanced, benchmark_contigs = load_and_enhance_features(
        args.input, args.features, args.benchmark_features, args.min_length, args.max_length
    )

    # Get numeric columns for correlation analysis
    numeric_cols = [col for col in df_enhanced.select_dtypes(include=[np.number]).columns.tolist()]

    print("=" * 60)
    print(f"ENHANCED CORRELATION ANALYSIS WITH LENGTH-BASED COLUMN MERGING")
    print("=" * 60)
    print(f"Normalized features: {args.input}")
    print(f"Non-benchmark features: {args.features}")
    print(f"Benchmark features: {args.benchmark_features}")
    print(f"Length range: {args.min_length}-{args.max_length} nt")
    print(f"Length columns: {args.max_length - args.min_length + 1}")
    print(f"Total dimensions: {len(numeric_cols)}")
    print(f"Total contigs: {len(df_enhanced)}")
    print(f"Benchmark contigs: {len(benchmark_contigs)}")
    print(f"Correlation method: PEARSON (fixed)")
    print(f"Output prefix: {args.output}")
    print(f"Output directory: {output_dir}")
    print("=" * 60)

    # Validate benchmark availability
    valid_benchmarks = [benchmark for benchmark in benchmark_contigs if benchmark in df_enhanced.index]
    print(f"Valid benchmark contigs: {len(valid_benchmarks)}/{len(benchmark_contigs)}")

    if not valid_benchmarks:
        print("Error: No valid benchmark contigs found in enhanced feature matrix")
        return

    # Calculate correlations for all benchmarks using PEARSON
    print(f"Calculating PEARSON correlations...")
    all_correlations = Parallel(n_jobs=args.threads)(
        delayed(calculate_all_correlations)(target_row, df_enhanced, numeric_cols)
        for target_row in valid_benchmarks
    )

    # Combine results from all benchmarks
    combined_scores = defaultdict(lambda: {'correlations': [], 'p_values': []})

    for i, (correlations, p_values) in enumerate(all_correlations):
        target_row = valid_benchmarks[i]
        for contig, corr in correlations.items():
            combined_scores[contig]['correlations'].append(corr)
            combined_scores[contig]['p_values'].append(p_values[contig])

    # Calculate mean correlation and combined p-value for each contig
    results = []
    for contig, scores in combined_scores.items():
        if scores['correlations'] and contig not in valid_benchmarks:  # Exclude benchmarks from results
            mean_corr = np.mean(scores['correlations'])

            # Combine p-values using Fisher's method
            clean_p_values = [max(min(p, 0.999), 0.001) for p in scores['p_values']]
            try:
                _, combined_p = combine_pvalues(clean_p_values, method='fisher')
            except:
                combined_p = min(1.0, len(clean_p_values) * min(clean_p_values))

            results.append({
                'Contig_Name': contig,
                'Mean_Pearson_Correlation': mean_corr,
                'Combined_P_Value': combined_p
            })

    # Sort by mean correlation in descending order
    results.sort(key=lambda x: x['Mean_Pearson_Correlation'], reverse=True)
    print(f"Total contigs analyzed: {len(results)}")
    print("No threshold filtering applied - saving all results")

    # Save correlation results to output directory
    results_df = pd.DataFrame(results)
    results_filename = os.path.join(output_dir, f"{args.output}_enhanced_pearson_correlation_results.csv")
    results_df.to_csv(results_filename, index=False)
    print(f"Correlation results saved as: {results_filename}")

    # Create enhanced heatmap
    normalized_data, distance_matrix, cluster_order = create_enhanced_heatmap(
        df_enhanced, valid_benchmarks, args.output, output_dir, args.min_length, args.max_length
    )

    # Save enhanced feature vectors to output directory
    enhanced_vectors_filename = os.path.join(output_dir, f"{args.output}_enhanced_feature_vectors.csv")
    normalized_data.to_csv(enhanced_vectors_filename)
    print(f"Enhanced feature vectors saved as: {enhanced_vectors_filename}")

    # Save distance matrix to output directory
    distance_df = pd.DataFrame(distance_matrix, index=df_enhanced.index, columns=df_enhanced.index)
    distance_filename = os.path.join(output_dir, f"{args.output}_enhanced_pearson_distance_matrix.csv")
    distance_df.to_csv(distance_filename)
    print(f"Distance matrix saved as: {distance_filename}")

    # Save clustering order to output directory
    cluster_df = pd.DataFrame({'Contig_Name': cluster_order, 'Cluster_Order': range(len(cluster_order))})
    cluster_filename = os.path.join(output_dir, f"{args.output}_enhanced_clustering_order.csv")
    cluster_df.to_csv(cluster_filename, index=False)
    print(f"Clustering order saved as: {cluster_filename}")

    # Save analysis summary to output directory
    summary_filename = os.path.join(output_dir, f"{args.output}_enhanced_analysis_summary.txt")
    with open(summary_filename, 'w') as f:
        f.write("Enhanced Correlation Analysis with Length-Based Column Merging\n")
        f.write("=" * 60 + "\n")
        f.write(f"Length range: {args.min_length}-{args.max_length} nt\n")
        f.write(f"Length columns: {args.max_length - args.min_length + 1}\n")
        f.write(f"Total dimensions: {len(numeric_cols)}\n")
        f.write(f"Correlation method: PEARSON\n")
        f.write(f"Total contigs: {len(df_enhanced)}\n")
        f.write(f"Benchmark contigs: {len(valid_benchmarks)}\n")
        f.write(f"All contigs analyzed: {len(results)}\n")
        f.write(f"Enhanced features calculation:\n")
        f.write(f"  AC_normalized = log10(Avg_Coverage + 1) / 10\n")
        f.write(f"  GC_content = GC_Content_percent / 100\n")
        f.write(f"  21-22-nt = sRNA_21_22nt_Percent / 100\n")
        f.write(f"Input files:\n")
        f.write(f"  All contigs feature vectors: {args.input}\n")
        f.write(f"  Non-benchmark genomic features: {args.features}\n")
        f.write(f"  Benchmark genomic features: {args.benchmark_features}\n")

    print(f"Analysis summary saved as: {summary_filename}")

    # Print top results
    print("\n" + "=" * 60)
    print("TOP 10 CORRELATED CONTIGS (NON-BENCHMARK)")
    print("=" * 60)
    for i, result in enumerate(results[:10], 1):
        print(f"{i:2d}. {result['Contig_Name']:20} "
              f"Mean Pearson: {result['Mean_Pearson_Correlation']:.4f} "
              f"P-value: {result['Combined_P_Value']:.2e}")

    print("\n" + "=" * 60)
    print("ENHANCED ANALYSIS COMPLETED SUCCESSFULLY")
    print("=" * 60)
    print(f"Results saved to directory: {output_dir}")
    print(f"All files use prefix: {args.output}")
    print(f"Length columns: {args.max_length - args.min_length + 1}")
    print(f"Total dimensions: {len(numeric_cols)}")
    print(f"Enhanced features: AC_normalized, GC_content, 21-22-nt")
    print(f"All contigs analyzed: {len(results)}")
    print(f"Correlation method: PEARSON")

    # List all output files
    print("\nOutput files:")
    for filename in os.listdir(output_dir):
        if filename.startswith(args.output):
            print(f"  - {filename}")


if __name__ == "__main__":
    main()