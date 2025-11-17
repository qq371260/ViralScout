#!/usr/bin/env python3
"""
Auto-select optimal benchmarks and vector mode based on correlation analysis
Updated with improved strand processing using correlation-based selection
and enhanced error handling for edge cases
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from itertools import combinations
import os
import argparse
from collections import defaultdict
import warnings
import re

warnings.filterwarnings('ignore')


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Auto-select optimal benchmarks and vector mode')
    parser.add_argument('--benchmarks', type=str, required=True,
                        help='Benchmark contigs file (txt or csv)')
    parser.add_argument('--vector_dir', type=str, default='.',
                        help='Directory containing vector files (default: current directory)')
    parser.add_argument('--output', type=str, default='optimal_benchmarks',
                        help='Output prefix (default: optimal_benchmarks)')
    parser.add_argument('--min_correlation', type=float, default=0.95,
                        help='Minimum correlation for benchmark selection')
    parser.add_argument('--min_correlation_mode', type=float, default=0.9,
                        help='Minimum correlation for mode selection')
    parser.add_argument('--target_correlation', type=float, default=0.95,
                        help='Target correlation for mode selection')
    parser.add_argument('--p_value', type=float, default=0.05,
                        help='P-value threshold')
    return parser.parse_args()


def extract_first_word(text):
    """Extract the first word from text, ignoring whitespace and common separators"""
    match = re.search(r'[^\s\t,;]+', text.strip())
    return match.group() if match else ""


def load_benchmarks(benchmarks_file):
    """Load benchmark contigs from file"""
    print(f"Loading benchmark contigs from: {benchmarks_file}")

    if benchmarks_file.endswith('.csv'):
        df = pd.read_csv(benchmarks_file)
        benchmark_names = df.iloc[:, 0].tolist()
    else:
        with open(benchmarks_file, 'r') as f:
            benchmark_names = [line.strip() for line in f if line.strip()]

    # Extract first word
    benchmark_names = [extract_first_word(name) for name in benchmark_names]
    print(f"Loaded {len(benchmark_names)} benchmark contigs")
    return benchmark_names


def load_vector_file(file_path, mode):
    """Load vector file and return DataFrame"""
    if not os.path.exists(file_path):
        return None

    print(f"Loading {mode} vectors from: {file_path}")
    df = pd.read_csv(file_path)
    df = df.set_index(df.columns[0])
    return df


def calculate_correlation(method, x, y):
    """Calculate correlation using specified method"""
    if method == 'spearman':
        return spearmanr(x, y)
    else:  # pearson
        return pearsonr(x, y)


def calculate_pairwise_correlations(df, contigs, min_correlation=0.95, p_value=0.05, method='pearson'):
    """Calculate pairwise correlations between contigs"""
    correlations = []
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

    for contig1, contig2 in combinations(contigs, 2):
        if contig1 not in df.index or contig2 not in df.index:
            continue

        try:
            vec1 = df.loc[contig1, numeric_cols].values.astype(float)
            vec2 = df.loc[contig2, numeric_cols].values.astype(float)

            # Remove NaN values
            mask = ~(np.isnan(vec1) | np.isnan(vec2))
            vec1_valid = vec1[mask]
            vec2_valid = vec2[mask]

            if len(vec1_valid) < 3:
                continue

            # Use selected correlation method
            corr, p_val = calculate_correlation(method, vec1_valid, vec2_valid)

            if not np.isnan(corr) and p_val < p_value:
                correlations.append({
                    'contig1': contig1,
                    'contig2': contig2,
                    'correlation': corr,
                    'p_value': p_val,
                    'method': method
                })
        except Exception as e:
            continue

    return correlations


def process_strand_vectors_correlation(df, selected_benchmarks, corr_method='pearson'):
    """Process strand vectors using correlation-based selection to unify strand polarity"""
    processed_df = df.copy()
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

    # Extract base names (remove _minus suffix if present)
    base_names = {}
    for contig in processed_df.index:
        if contig.endswith('_minus'):
            base_name = contig[:-6]
        else:
            base_name = contig
        if base_name not in base_names:
            base_names[base_name] = []
        base_names[base_name].append(contig)

    # Get valid target rows that exist in the dataframe
    valid_target_rows = [benchmark for benchmark in selected_benchmarks
                         if benchmark in processed_df.index or f"{benchmark}_minus" in processed_df.index]

    if not valid_target_rows:
        print("No valid target rows found for strand processing")
        return processed_df

    # Step 1: Select reference benchmark (first one, prefer positive strand)
    selected_strands = []

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

        if positive_strand and positive_strand in processed_df.index:
            selected_strands.append(positive_strand)
            print(f"  First benchmark: selected positive strand '{positive_strand}'")
        elif minus_strand and minus_strand in processed_df.index:
            selected_strands.append(minus_strand)
            print(f"  First benchmark: using minus strand '{minus_strand}' (no positive available)")
        else:
            # Fallback: use original if available
            if first_benchmark in processed_df.index:
                selected_strands.append(first_benchmark)
                print(f"  First benchmark: using original '{first_benchmark}' (fallback)")
    else:
        # Fallback: use original if available
        if first_benchmark in processed_df.index:
            selected_strands.append(first_benchmark)
            print(f"  First benchmark: using original '{first_benchmark}' (fallback)")

    # Step 2: For other benchmarks, select strand with highest correlation to first benchmark
    for i, benchmark in enumerate(valid_target_rows[1:], 1):
        benchmark_base = benchmark
        if benchmark_base.endswith('_minus'):
            benchmark_base = benchmark_base[:-6]

        if benchmark_base in base_names:
            available_strands = base_names[benchmark_base]
            strand_scores = {}

            for strand in available_strands:
                if strand in processed_df.index and selected_strands:  # Ensure we have reference
                    try:
                        # Calculate correlation with first benchmark
                        first_bench_vector = processed_df.loc[selected_strands[0], numeric_cols].values.astype(float)
                        strand_vector = processed_df.loc[strand, numeric_cols].values.astype(float)

                        valid_mask = ~(np.isnan(first_bench_vector) | np.isnan(strand_vector))
                        first_valid = first_bench_vector[valid_mask]
                        strand_valid = strand_vector[valid_mask]

                        if len(first_valid) > 3:
                            corr, _ = calculate_correlation(corr_method, first_valid, strand_valid)
                            if not np.isnan(corr):
                                strand_scores[strand] = corr
                    except Exception as e:
                        print(f"  Error calculating correlation for {strand}: {e}")

            if strand_scores:
                # Select strand with highest correlation to first benchmark
                best_strand = max(strand_scores.items(), key=lambda x: x[1])
                selected_strands.append(best_strand[0])
                print(f"  Benchmark {i}: selected '{best_strand[0]}' (corr: {best_strand[1]:.4f})")
            else:
                # Fallback: use positive strand if available
                positive_strand = None
                for strand in available_strands:
                    if not strand.endswith('_minus') and strand in processed_df.index:
                        positive_strand = strand
                        break
                if positive_strand:
                    selected_strands.append(positive_strand)
                    print(f"  Benchmark {i}: using positive strand '{positive_strand}' (fallback)")
                elif available_strands and available_strands[0] in processed_df.index:
                    selected_strands.append(available_strands[0])
                    print(f"  Benchmark {i}: using '{available_strands[0]}' (fallback)")
        else:
            # Fallback: use original if available
            if benchmark in processed_df.index:
                selected_strands.append(benchmark)
                print(f"  Benchmark {i}: using original '{benchmark}' (fallback)")

    print(f"  Selected {len(selected_strands)} benchmark strands")

    # Create new dataframe with only selected strands, renaming minus strands to base names
    new_index_map = {}
    strands_to_keep = set()

    for strand in selected_strands:
        if strand.endswith('_minus'):
            base_name = strand[:-6]
            new_index_map[strand] = base_name
        else:
            new_index_map[strand] = strand
        strands_to_keep.add(strand)

    # Filter dataframe to only include selected strands
    filtered_df = processed_df.loc[list(strands_to_keep)].copy()

    # Rename minus strands to base names
    filtered_df = filtered_df.rename(index=new_index_map)

    # Also include non-benchmark contigs that don't have strand variants
    non_strand_contigs = [contig for contig in processed_df.index
                          if not any(contig.endswith(suffix) for suffix in ['_minus'])]
    for contig in non_strand_contigs:
        if contig not in filtered_df.index:
            filtered_df = pd.concat([filtered_df, processed_df.loc[[contig]]])

    return filtered_df


def select_optimal_benchmarks(size_df, benchmarks, min_correlation=0.95, p_value=0.05):
    """Select optimal benchmarks based on pairwise correlations"""
    print("\n" + "=" * 60)
    print("Step 1: Selecting optimal benchmarks using size vectors")
    print("=" * 60)

    # Calculate all pairwise correlations using both Pearson and Spearman
    all_correlations_pearson = calculate_pairwise_correlations(
        size_df, benchmarks, min_correlation, p_value, method='pearson'
    )

    all_correlations_spearman = calculate_pairwise_correlations(
        size_df, benchmarks, min_correlation, p_value, method='spearman'
    )

    # Choose the method that gives more benchmarks
    if len(all_correlations_pearson) >= len(all_correlations_spearman):
        all_correlations = all_correlations_pearson
        selected_method = 'pearson'
    else:
        all_correlations = all_correlations_spearman
        selected_method = 'spearman'
    print(
        f"Selected {selected_method.upper()} correlation for benchmark selection (found {len(all_correlations)} valid pairs)")

    # Group correlations by contig
    contig_stats = defaultdict(lambda: {'count': 0, 'sum_corr': 0, 'pairs': []})

    for corr_data in all_correlations:
        contig1, contig2 = corr_data['contig1'], corr_data['contig2']
        correlation = corr_data['correlation']

        contig_stats[contig1]['count'] += 1
        contig_stats[contig1]['sum_corr'] += correlation
        contig_stats[contig1]['pairs'].append((contig2, correlation))

        contig_stats[contig2]['count'] += 1
        contig_stats[contig2]['sum_corr'] += correlation
        contig_stats[contig2]['pairs'].append((contig1, correlation))

    # Find contig with maximum number of high-correlation pairs
    if contig_stats:  # Ensure contig_stats is not empty
        max_count = max(stats['count'] for stats in contig_stats.values())
    else:
        max_count = 0

    # Get all contigs with maximum count
    candidates = [contig for contig, stats in contig_stats.items()
                  if stats['count'] == max_count]

    if len(candidates) == 0:
        # No candidates found (all correlations below threshold), select first benchmark
        optimal_contig = benchmarks[0]
        print(f"No high-correlation pairs found. Selected first benchmark: {optimal_contig}")
    elif len(candidates) == 1:
        optimal_contig = candidates[0]
        print(f"Selected optimal benchmark: {optimal_contig} with {max_count} high-correlation pairs")
    else:
        # If multiple candidates, select the one with highest average correlation
        best_avg_corr = -1
        optimal_contig = None

        for contig in candidates:
            if contig_stats[contig]['count'] > 0:
                avg_corr = contig_stats[contig]['sum_corr'] / contig_stats[contig]['count']
            else:
                avg_corr = 0  # If no correlation pairs, set average correlation to 0

            if avg_corr > best_avg_corr:
                best_avg_corr = avg_corr
                optimal_contig = contig

        print(
            f"Selected optimal benchmark: {optimal_contig} with {max_count} high-correlation pairs (avg correlation: {best_avg_corr:.4f})")

    # Get all contigs that have high correlation with the optimal contig
    optimal_pairs = contig_stats[optimal_contig]['pairs']
    selected_benchmarks = set([optimal_contig])

    for other_contig, correlation in optimal_pairs:
        if correlation >= min_correlation:
            selected_benchmarks.add(other_contig)

    # Convert to sorted list
    selected_benchmarks = sorted(list(selected_benchmarks))

    print(f"\nSelected {len(selected_benchmarks)} benchmarks:")
    for i, contig in enumerate(selected_benchmarks, 1):
        print(f"  {i}. {contig}")

    # Print rejected benchmarks
    rejected = set(benchmarks) - set(selected_benchmarks)
    if rejected:
        print(f"\nRejected {len(rejected)} benchmarks:")
        for i, contig in enumerate(sorted(list(rejected)), 1):
            print(f"  {i}. {contig}")

    return selected_benchmarks


def evaluate_vector_mode(df, selected_benchmarks, min_correlation=0.9, target_correlation=0.95, method='pearson'):
    """Evaluate vector mode based on benchmark correlations"""
    print(f"  Evaluating {len(selected_benchmarks)} benchmarks using {method.upper()} correlation...")

    # Calculate pairwise correlations
    correlations = calculate_pairwise_correlations(
        df, selected_benchmarks,
        min_correlation=min_correlation,
        method=method
    )

    if not correlations:
        print(f"  No valid correlations found for this mode")
        return None, None, False

    # Check if all pairs meet minimum correlation
    all_pairs_valid = True
    benchmark_pairs = list(combinations(selected_benchmarks, 2))
    min_corr = 1.0

    for contig1, contig2 in benchmark_pairs:
        pair_found = False
        for corr_data in correlations:
            if (corr_data['contig1'] == contig1 and corr_data['contig2'] == contig2) or \
                    (corr_data['contig1'] == contig2 and corr_data['contig2'] == contig1):
                pair_found = True
                corr_value = corr_data['correlation']
                if corr_value < min_corr:
                    min_corr = corr_value
                if corr_value < min_correlation:
                    all_pairs_valid = False
                break
        if not pair_found:
            all_pairs_valid = False

    if not all_pairs_valid:
        print(f"  Not all benchmark pairs meet minimum correlation ({min_correlation})")
        print(f"  Minimum correlation found: {min_corr:.4f}")
        return None, None, False

    # Calculate average correlation
    avg_correlation = np.mean([corr['correlation'] for corr in correlations])

    # Calculate distance from target correlation
    distance_from_target = abs(avg_correlation - target_correlation)

    print(f"  All pairs valid: ✓")
    print(f"  Minimum correlation found: {min_corr:.4f}")
    print(f"  Average correlation: {avg_correlation:.4f}")
    print(f"  Distance from target ({target_correlation}): {distance_from_target:.4f}")

    return avg_correlation, distance_from_target, True


def select_optimal_mode(vector_files, selected_benchmarks, min_correlation=0.9, target_correlation=0.95):
    """Select optimal vector mode and correlation method based on benchmark correlations"""
    print("\n" + "=" * 60)
    print("Step 2: Selecting optimal vector mode and correlation method")
    print("=" * 60)

    # Handle single benchmark case
    if len(selected_benchmarks) < 2:
        print(f"Warning: Only {len(selected_benchmarks)} benchmark(s) selected.")
        print("Cannot perform correlation-based mode selection with less than 2 benchmarks.")
        print("Using default: size mode with Pearson correlation")

        # Verify size mode file exists
        if os.path.exists(vector_files['size']):
            return 'size', 'pearson'
        else:
            # If size mode doesn't exist, try other modes
            for mode, file_path in vector_files.items():
                if os.path.exists(file_path):
                    print(f"Fallback: Using {mode} mode with Pearson correlation")
                    return mode, 'pearson'
            print("Error: No vector files found!")
            return None, None

    all_results = []

    # Evaluate all 8 combinations: 4 modes × 2 correlation methods
    for mode, file_path in vector_files.items():
        print(f"\nEvaluating {mode} mode...")

        df = load_vector_file(file_path, mode)
        if df is None:
            print(f"  Vector file not found: {file_path}")
            continue

        # Process strand vectors for strand-specific modes using correlation-based method
        if mode in ['sizeXstr', 'sizeX5ntXstr']:
            print(f"  Processing strand information using correlation-based selection...")
            # Use Pearson for strand selection since we're comparing with first benchmark
            df = process_strand_vectors_correlation(df, selected_benchmarks, corr_method='spearman')

        # Try both correlation methods
        for method in ['pearson', 'spearman']:
            print(f"  Testing {method.upper()} correlation...")

            # Evaluate the mode with this correlation method
            avg_corr, distance, valid = evaluate_vector_mode(
                df, selected_benchmarks, min_correlation, target_correlation, method
            )

            if valid:
                all_results.append({
                    'mode': mode,
                    'method': method,
                    'avg_correlation': avg_corr,
                    'distance_from_target': distance,
                    'valid': True
                })

    # Select best combination
    if not all_results:
        print(f"\nNo valid vector mode and correlation method combinations found!")
        print("Falling back to default mode...")

        # Try size mode first
        if os.path.exists(vector_files['size']):
            best_mode = 'size'
            best_method = 'pearson'
            print(f"\nFALLBACK: Using size mode with Pearson correlation")
        else:
            # If size mode doesn't exist, try other available modes
            for mode, file_path in vector_files.items():
                if os.path.exists(file_path):
                    best_mode = mode
                    best_method = 'pearson'
                    print(f"\nFALLBACK: Using {mode} mode with Pearson correlation")
                    break
            else:
                print(f"\nERROR: No vector files found for fallback!")
                return None, None
    else:
        # Select combination with distance closest to 0 (closest to target correlation)
        best_combination = min(all_results, key=lambda x: x['distance_from_target'])
        best_mode = best_combination['mode']
        best_method = best_combination['method']
        best_distance = best_combination['distance_from_target']

    print(f"\n" + "=" * 40)
    print("OPTIMAL COMBINATION SELECTED:")
    print("=" * 40)
    print(f"Vector Mode: {best_mode}")
    print(f"Correlation Method: {best_method.upper()}")
    if best_mode != 'size' or best_method != 'pearson':  # Only show correlation if not fallback
        print(f"Average correlation: {best_combination.get('avg_correlation', 'N/A'):.4f}")
        print(f"Distance from target ({target_correlation}): {best_distance:.4f}")

    # Print all combinations rankings (only if we have results)
    if all_results:
        print(f"\nAll combinations rankings (by distance from target):")
        sorted_combinations = sorted(all_results, key=lambda x: x['distance_from_target'])
        for i, combo in enumerate(sorted_combinations, 1):
            print(f"  {i}. {combo['mode']} + {combo['method'].upper()}: "
                  f"avg_corr={combo['avg_correlation']:.4f}, distance={combo['distance_from_target']:.4f}")

    return best_mode, best_method


def main():
    args = parse_arguments()

    # Load benchmarks
    benchmarks = load_benchmarks(args.benchmarks)

    # Define vector files
    vector_files = {
        'size': os.path.join(args.vector_dir, 'size_vector.csv'),
        'size_P_5nt': os.path.join(args.vector_dir, 'size_P_5nt_vector.csv'),
        'sizeXstr': os.path.join(args.vector_dir, 'sizeXstr_vector.csv'),
        'sizeX5nt': os.path.join(args.vector_dir, 'sizeX5nt_vector.csv'),
        'sizeX5ntXstr': os.path.join(args.vector_dir, 'sizeX5ntXstr_vector.csv')
    }

    # Step 1: Load size vectors and select optimal benchmarks
    size_df = load_vector_file(vector_files['size'], 'size')
    if size_df is None:
        print(f"Error: Size vector file not found: {vector_files['size']}")
        return

    selected_benchmarks = select_optimal_benchmarks(
        size_df, benchmarks,
        min_correlation=args.min_correlation,
        p_value=args.p_value
    )

    # Save selected benchmarks
    output_file = f"{args.output}_benchmarks.txt"
    with open(output_file, 'w') as f:
        for contig in selected_benchmarks:
            f.write(f"{contig}\n")
    print(f"\nSelected benchmarks saved to: {output_file}")

    # Step 2: Select optimal vector mode and correlation method
    best_mode, best_method = select_optimal_mode(
        vector_files, selected_benchmarks,
        min_correlation=args.min_correlation_mode,
        target_correlation=args.target_correlation
    )

    # Save mode selection results
    if best_mode:
        mode_output_file = f"{args.output}_optimal_mode.txt"
        with open(mode_output_file, 'w') as f:
            f.write(f"OPTIMAL_VECTOR_MODE: {best_mode}\n")
            f.write(f"OPTIMAL_CORRELATION_METHOD: {best_method}\n")
            f.write(f"SELECTED_BENCHMARKS: {len(selected_benchmarks)}\n")
            f.write(f"BENCHMARK_FILE: {output_file}\n")
        print(f"\nOptimal combination information saved to: {mode_output_file}")

    print("\n" + "=" * 60)
    print("Analysis completed successfully!")
    print("=" * 60)


if __name__ == '__main__':
    main()