#!/usr/bin/env python3
"""
RNA_or_sRNA-seq Differential Expression Analysis
Generate TPM and Fold Change data with visualization
Add benchmark contigs distance calculation function
Support starting analysis from existing results files
Simplified version: directly use read counts and contig length to calculate TPM
Retain all contigs for analysis
"""

import argparse
import numpy as np
import pandas as pd
import subprocess
import os
import sys
import matplotlib.pyplot as plt
from scipy.spatial import distance
from matplotlib import font_manager


class SeqAnalyzer:
    def __init__(self, positive_bam=None, negative_bam=None, contig_length_file=None,
                 benchmark_contigs_file=None, existing_results_file=None, nearest_percent=0.1):
        self.positive_bam = positive_bam
        self.negative_bam = negative_bam
        self.contig_lengths = self._load_contig_lengths(contig_length_file) if contig_length_file else None
        self.benchmark_contigs = self._load_benchmark_contigs(
            benchmark_contigs_file) if benchmark_contigs_file else None
        self.existing_results_file = existing_results_file
        self.nearest_percent = nearest_percent

        self.contig_names = None
        self.read_count_data = None
        self.tpm_data = None
        self.fold_changes = None
        self.zero_flags = None
        self.benchmark_indices = None
        self.nearest_indices = None
        self.distances = None

    def _load_contig_lengths(self, length_file):
        """Load contig length file - for TPM calculation"""
        if not length_file or not os.path.isfile(length_file):
            raise ValueError("Contig length file is required for TPM calculation")

        lengths = {}
        with open(length_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    lengths[parts[0]] = int(parts[1])
        print(f"Loaded length information for {len(lengths)} contigs")

        # Statistics on length distribution
        if lengths:
            lengths_array = np.array(list(lengths.values()))
            print(
                f"Contig length statistics: min={np.min(lengths_array)}nt, max={np.max(lengths_array)}nt, mean={np.mean(lengths_array):.1f}nt")

        return lengths

    def _load_benchmark_contigs(self, benchmark_file):
        """Load benchmark contigs file"""
        if not benchmark_file or not os.path.isfile(benchmark_file):
            print("Warning: No benchmark contigs file provided")
            return set()

        benchmark_contigs = set()
        with open(benchmark_file, 'r') as f:
            for line in f:
                contig = line.strip()
                if contig:
                    benchmark_contigs.add(contig)

        print(f"Loaded {len(benchmark_contigs)} benchmark contigs")
        return benchmark_contigs

    def load_existing_results(self, results_file):
        """Load data from existing results file"""
        print(f"Loading data from existing results file: {results_file}")

        if not os.path.isfile(results_file):
            raise FileNotFoundError(f"Results file does not exist: {results_file}")

        df = pd.read_csv(results_file)

        # Extract necessary data
        self.contig_names = df['contig'].values
        self.tpm_data = np.column_stack([df['negative_tpm'].values, df['positive_tpm'].values])
        self.fold_changes = df['log2_fold_change'].values

        # Optional fields
        if 'negative_read_count' in df.columns and 'positive_read_count' in df.columns:
            self.read_count_data = np.column_stack([df['negative_read_count'].values, df['positive_read_count'].values])

        # Mark zero values (if available)
        if 'has_zero_coverage' in df.columns:
            self.zero_flags = df['has_zero_coverage'].values
        else:
            # If no zero value markers, determine based on TPM values
            self.zero_flags = (self.tpm_data[:, 0] == 0) | (self.tpm_data[:, 1] == 0)

        print(f"Successfully loaded data for {len(self.contig_names)} contigs")
        print(f"TPM range: {np.min(self.tpm_data):.2f} - {np.max(self.tpm_data):.2f}")
        print(f"Fold Change range: {np.min(self.fold_changes):.2f} - {np.max(self.fold_changes):.2f}")

        return df

    def calculate_total_mapped_reads(self, bam_file):
        """Calculate total mapped reads in BAM file"""
        cmd = f"samtools flagstat {bam_file} | grep 'mapped (' | head -1 | cut -d' ' -f1"
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
            total_reads = int(result.stdout.strip())
            print(f"  {bam_file}: {total_reads:,} mapped reads")
            return total_reads
        except:
            print(f"Warning: Could not calculate total mapped reads for {bam_file}, using default value 1e6")
            return 1000000

    def calculate_read_counts_and_tpm(self, bam_file, total_reads=None):
        """
        Use samtools idxstats to get read counts, then calculate TPM
        TPM calculation steps:
          1. RPK = read_count / (contig_length / 1000)
          2. RPK_sum = sum of RPK for all contigs
          3. TPM = (RPK / RPK_sum) * 1e6
        """
        if total_reads is None:
            total_reads = self.calculate_total_mapped_reads(bam_file)

        print(f"Processing seq BAM file: {bam_file}")

        # Use samtools idxstats to get read count for each contig
        print("  Calculating read counts using samtools idxstats...")
        cmd_idxstats = f"samtools idxstats {bam_file}"

        read_count_dict = {}
        try:
            result = subprocess.run(cmd_idxstats, shell=True, capture_output=True, text=True, check=True)
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        contig = parts[0]
                        if contig != '*':  # Skip unmapped reads
                            mapped_reads = int(parts[2])
                            read_count_dict[contig] = mapped_reads  # Record all contigs, including 0
        except subprocess.CalledProcessError as e:
            print(f"Error using samtools idxstats: {e}")
            raise

        # Calculate TPM
        tpm_dict = {}
        total_rpk = 0
        rpk_dict = {}

        for contig, read_count in read_count_dict.items():
            length = self.contig_lengths.get(contig, 0)
            if length > 0:
                # RPK = reads per kilobase
                rpk = read_count / (length / 1000)
                rpk_dict[contig] = rpk
                total_rpk += rpk

        # TPM = (RPK / total_RPK) * 1e6
        scaling_factor = 1e6 / total_rpk if total_rpk > 0 else 0

        for contig, rpk in rpk_dict.items():
            tpm_dict[contig] = rpk * scaling_factor

        print(f"  Detected {len(read_count_dict)} contigs")
        print(f"  Total RPK: {total_rpk:.2f}, calculated TPM for {len(tpm_dict)} contigs")

        return tpm_dict, read_count_dict

    def process_seq_samples(self):
        """Process seq sample data, retain all contigs"""
        print("Processing seq positive/negative sample data...")

        # Calculate total mapped reads
        print("Calculating total mapped reads:")
        positive_total_reads = self.calculate_total_mapped_reads(self.positive_bam)
        negative_total_reads = self.calculate_total_mapped_reads(self.negative_bam)

        # Calculate TPM and read count
        positive_tpm, positive_reads = self.calculate_read_counts_and_tpm(
            self.positive_bam, positive_total_reads
        )
        negative_tpm, negative_reads = self.calculate_read_counts_and_tpm(
            self.negative_bam, negative_total_reads
        )

        all_contigs = sorted(
            set(positive_tpm.keys()) |
            set(negative_tpm.keys())
        )

        print(f"Total contigs found: {len(all_contigs)}")
        print(f"  Only in positive sample: {len(set(positive_tpm.keys()) - set(negative_tpm.keys()))}")
        print(f"  Only in negative sample: {len(set(negative_tpm.keys()) - set(positive_tpm.keys()))}")
        print(f"  In both samples: {len(set(positive_tpm.keys()) & set(negative_tpm.keys()))}")

        # Organize data
        self.contig_names = np.array(all_contigs)
        tpm_data = []
        read_count_data = []
        zero_flags = []

        for contig in all_contigs:
            pos_tpm = positive_tpm.get(contig, 0)
            neg_tpm = negative_tpm.get(contig, 0)
            tpm_data.append([neg_tpm, pos_tpm])  # Negative first, positive second

            pos_reads = positive_reads.get(contig, 0)
            neg_reads = negative_reads.get(contig, 0)
            read_count_data.append([neg_reads, pos_reads])

            # Mark zero values
            has_zero = pos_tpm == 0 or neg_tpm == 0
            zero_flags.append(has_zero)

        self.tpm_data = np.array(tpm_data)
        self.read_count_data = np.array(read_count_data)
        self.zero_flags = np.array(zero_flags)

        # Calculate fold change: positive/negative
        pseudo_count = 0.1
        fc_data = []

        for i in range(len(all_contigs)):
            # FC = positive TPM / negative TPM
            fc = (self.tpm_data[i, 1] + pseudo_count) / (self.tpm_data[i, 0] + pseudo_count)
            log2_fc = np.log2(fc)
            fc_data.append(log2_fc)

        self.fold_changes = np.array(fc_data)

        # Statistics
        n_zero = sum(self.zero_flags)
        print(f"Successfully processed {len(all_contigs)} contigs")
        print(f"Contigs with zero values: {n_zero} ({n_zero / len(all_contigs) * 100:.1f}%)")
        print(f"TPM range: {np.min(self.tpm_data):.2f} - {np.max(self.tpm_data):.2f}")
        print(f"Fold Change range: {np.min(self.fold_changes):.2f} - {np.max(self.fold_changes):.2f}")

        return True

    def calculate_benchmark_distances(self):
        """Calculate distances of all contigs to benchmark contigs"""
        if not self.benchmark_contigs:
            print("Warning: No benchmark contigs provided, skipping distance calculation")
            return None

        # Find indices of benchmark contigs in data
        self.benchmark_indices = []
        benchmark_found = []

        for i, contig in enumerate(self.contig_names):
            if contig in self.benchmark_contigs:
                self.benchmark_indices.append(i)
                benchmark_found.append(contig)

        print(f"Found {len(self.benchmark_indices)}/{len(self.benchmark_contigs)} benchmark contigs in data")

        if len(self.benchmark_indices) == 0:
            print("Error: No benchmark contigs found in data")
            return None

        # Calculate features: log(total TPM) and log FC
        total_tpm = np.sum(self.tpm_data, axis=1)
        log_total_tpm = np.log10(total_tpm + 1)  # Add 1 to avoid log(0)
        log_fc = self.fold_changes

        # Create feature matrix: [log(total TPM), log FC]
        features = np.column_stack([log_total_tpm, log_fc])

        # Calculate feature center of benchmark contigs
        benchmark_features = features[self.benchmark_indices]
        benchmark_center = np.mean(benchmark_features, axis=0)

        print(f"Benchmark contigs center: log_total_TPM={benchmark_center[0]:.3f}, log_FC={benchmark_center[1]:.3f}")

        # Calculate Euclidean distance of all contigs to benchmark center
        self.distances = np.array([distance.euclidean(features[i], benchmark_center)
                                   for i in range(len(features))])

        print(f"Distance range: {np.min(self.distances):.4f} - {np.max(self.distances):.4f}")

        return self.distances

    def select_nearest_contigs(self):
        """Select top_percent% contigs with smallest distances (excluding benchmark contigs)"""
        if self.distances is None:
            print("Error: Please run calculate_benchmark_distances() first")
            return None

        # Exclude benchmark contigs
        non_benchmark_mask = np.ones(len(self.distances), dtype=bool)
        non_benchmark_mask[self.benchmark_indices] = False
        non_benchmark_distances = self.distances[non_benchmark_mask]
        non_benchmark_indices = np.where(non_benchmark_mask)[0]

        # Calculate number of contigs to select
        n_select = max(1, int(len(non_benchmark_distances) * self.nearest_percent / 100))

        # Select contigs with smallest distances
        nearest_indices_in_subset = np.argsort(non_benchmark_distances)[:n_select]
        self.nearest_indices = non_benchmark_indices[nearest_indices_in_subset]

        print(f"Selected {n_select} nearest contigs from {len(non_benchmark_distances)} non-benchmark contigs ({self.nearest_percent}%)")
        print(
            f"Nearest contigs distance range: {np.min(self.distances[self.nearest_indices]):.4f} - {np.max(self.distances[self.nearest_indices]):.4f}")

        return self.nearest_indices

    def create_enhanced_visualization(self, output_prefix):
        """Create enhanced visualization, retaining only feature space and distance visualization"""
        font_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "times.ttf")
        if os.path.exists(font_path):
            font_manager.fontManager.addfont(font_path)
            # Set global font to Times New Roman
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['font.serif'] = ['Times New Roman']
            plt.rcParams['mathtext.fontset'] = 'stix'
        plt.rcParams['font.size'] = 14

        # Prepare data
        total_tpm = np.sum(self.tpm_data, axis=1)
        log_total_tpm = np.log10(total_tpm + 1)
        log2_fc = self.fold_changes

        # Create figure - retain only one subplot
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))

        if self.distances is not None:
            # Regular contigs
            scatter = ax.scatter(
                log_total_tpm, log2_fc,
                c=self.distances, cmap='viridis', alpha=0.6, s=30
            )

            # Mark benchmark contigs
            if self.benchmark_indices is not None:
                ax.scatter(
                    log_total_tpm[self.benchmark_indices], log2_fc[self.benchmark_indices],
                    c='red', alpha=0.6, s=80, marker='*', label='Benchmark contigs'
                )

            # Mark nearest contigs
            if self.nearest_indices is not None:
                ax.scatter(
                    log_total_tpm[self.nearest_indices], log2_fc[self.nearest_indices],
                    c='orange', alpha=0.8, s=30, marker='^', label=f'Nearest {self.nearest_percent}% contigs'
                )

            # Mark benchmark center
            benchmark_center_log_tpm = np.mean(log_total_tpm[self.benchmark_indices])
            benchmark_center_log_fc = np.mean(log2_fc[self.benchmark_indices])
            ax.scatter(
                [benchmark_center_log_tpm], [benchmark_center_log_fc],
                c='black', alpha=0.6, s=80, marker='X', label='Benchmark center'
            )

            # Set color bar and labels
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('Distance to benchmark center', fontsize=16)

            # Set axis labels and title
            ax.set_xlabel('Log10(TPM + 1)', fontsize=16)
            ax.set_ylabel('Log2(FC)', fontsize=16)
            ax.set_title('', fontsize=20)

            # Set legend
            legend = ax.legend()
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(f"{output_prefix}_feature_space.svg", bbox_inches='tight')
        plt.show()

    def run_complete_seq_analysis(self, output_prefix):
        """Run complete seq analysis"""
        # Process data (if BAM files provided) or load existing results
        if self.existing_results_file:
            df = self.load_existing_results(self.existing_results_file)
        else:
            if not self.positive_bam or not self.negative_bam:
                raise ValueError("Must provide BAM files or existing results file")
            if not self.contig_lengths:
                raise ValueError("Must provide contig length file for TPM calculation")

            self.process_seq_samples()

            # Save processed data
            data_dict = {
                'contig': self.contig_names,
                'length': [self.contig_lengths.get(c, 0) for c in self.contig_names],
                'negative_read_count': self.read_count_data[:, 0],
                'positive_read_count': self.read_count_data[:, 1],
                'negative_tpm': self.tpm_data[:, 0],
                'positive_tpm': self.tpm_data[:, 1],
                'log2_fold_change': self.fold_changes,
                'has_zero_coverage': self.zero_flags
            }

            df = pd.DataFrame(data_dict)

            # Add differential expression status
            df['DE_status'] = 'Not DE'
            df.loc[df['log2_fold_change'] > 1, 'DE_status'] = 'Up-regulated'
            df.loc[df['log2_fold_change'] < -1, 'DE_status'] = 'Down-regulated'

            # Sorted based on |fold change|
            df['abs_log2_fc'] = np.abs(df['log2_fold_change'])
            df = df.sort_values('abs_log2_fc', ascending=False)
            df = df.drop('abs_log2_fc', axis=1)

            df.to_csv(f"{output_prefix}_seq_data.csv", index=False)

        # If benchmark contigs available, calculate distances and select nearest contigs
        if self.benchmark_contigs:
            print("\n=== Benchmark Contigs Distance Analysis ===")
            distances = self.calculate_benchmark_distances()
            if distances is not None:
                # Select nearest contigs
                nearest_indices = self.select_nearest_contigs()

                # Add distance information — map by contig name to preserve correct alignment
                # Create a Series mapping contig -> distance
                dist_series = pd.Series(self.distances, index=self.contig_names)

                # Create a Series mapping contig -> distance
                dist_series = pd.Series(self.distances, index=self.contig_names)

                # Map distances to df by contig name (regardless of df sorting)
                df['distance_to_benchmark'] = df['contig'].map(dist_series)

                # Initialize flags
                df['is_benchmark'] = False
                df['is_nearest'] = False

                # Mark benchmark contigs by contig name
                if self.benchmark_indices is not None and len(self.benchmark_indices) > 0:
                    benchmark_contigs_in_data = list(self.contig_names[self.benchmark_indices])
                    df.loc[df['contig'].isin(benchmark_contigs_in_data), 'is_benchmark'] = True

                # Mark nearest contigs by contig name
                if self.nearest_indices is not None and len(self.nearest_indices) > 0:
                    nearest_contigs_in_data = list(self.contig_names[self.nearest_indices])
                    df.loc[df['contig'].isin(nearest_contigs_in_data), 'is_nearest'] = True

                # --- Reorder df to match the original contig order used in the figure ---
                df['order_index'] = df['contig'].map({c: i for i, c in enumerate(self.contig_names)})
                df = df.sort_values('order_index').drop(columns=['order_index'])

                # Save results with distance info (in same order as figure)
                df.to_csv(f"{output_prefix}_with_distances.csv", index=False)

                # Save nearest contigs file (keep same order as figure, not sorted by distance)
                if self.nearest_indices is not None and len(self.nearest_indices) > 0:
                    nearest_contigs = list(self.contig_names[self.nearest_indices])
                    nearest_df = df[df['contig'].isin(nearest_contigs)].copy()
                    nearest_df.to_csv(f"{output_prefix}_nearest_contigs.csv", index=False)

                    # Save list of nearest contigs in same figure order
                    with open(f"{output_prefix}_nearest_contigs.txt", 'w') as f:
                        for contig in nearest_contigs:
                            f.write(f"{contig}\n")

                # Save benchmark contigs file (keep same order as figure)
                if self.benchmark_indices is not None and len(self.benchmark_indices) > 0:
                    benchmark_contigs = list(self.contig_names[self.benchmark_indices])
                    benchmark_df = df[df['contig'].isin(benchmark_contigs)].copy()
                    benchmark_df.to_csv(f"{output_prefix}_benchmark_contigs.csv", index=False)

        # Create enhanced visualization
        self.create_enhanced_visualization(output_prefix)

        print(f"\nseq Analysis Complete!")
        print(f"Total contigs: {len(df)}")

        if 'DE_status' in df.columns:
            up_regulated = df[df['DE_status'] == 'Up-regulated']
            down_regulated = df[df['DE_status'] == 'Down-regulated']
            print(f"Up-regulated contigs (FC > 2): {len(up_regulated)}")
            print(f"Down-regulated contigs (FC < 0.5): {len(down_regulated)}")
            print(f"Non-differentially expressed contigs: {len(df) - len(up_regulated) - len(down_regulated)}")

        if self.benchmark_indices is not None:
            print(f"Benchmark contigs: {len(self.benchmark_indices)}")
        if self.nearest_indices is not None:
            print(f"Nearest contigs ({self.nearest_percent}%): {len(self.nearest_indices)}")

        return df


def main():
    parser = argparse.ArgumentParser(description='seq Differential Expression Analysis')
    parser.add_argument('--positive_bam', help='Positive sample BAM file')
    parser.add_argument('--negative_bam', help='Negative sample BAM file')
    parser.add_argument('--contig_lengths', required=True, help='Contig length file (required)')
    parser.add_argument('--benchmark_contigs', required=True, help='Benchmark contigs file')
    parser.add_argument('--existing_results', help='Existing results file (if provided, skip BAM processing)')
    parser.add_argument('--output', required=True, help='Output file prefix')
    parser.add_argument('--nearest_percent', type=float, default=0.1,
                       help='Percentage of nearest contigs to select (default: 0.1)')

    args = parser.parse_args()

    # Check parameters
    if not args.existing_results and (not args.positive_bam or not args.negative_bam):
        parser.error("Must provide --existing_results or both --positive_bam and --negative_bam")

    # Create seq analyzer
    analyzer = SeqAnalyzer(
        positive_bam=args.positive_bam,
        negative_bam=args.negative_bam,
        contig_length_file=args.contig_lengths,
        benchmark_contigs_file=args.benchmark_contigs,
        existing_results_file=args.existing_results,
        nearest_percent=args.nearest_percent
    )

    # Run analysis
    try:
        processed_df = analyzer.run_complete_seq_analysis(args.output)

        # Output summary
        print(f"\nTop 10 most significantly differentially expressed contigs:")
        for i, row in processed_df.head(10).iterrows():
            length_info = f" [{row['length']}nt]" if 'length' in row and row['length'] > 0 else ""
            zero_marker = " [ZERO]" if 'has_zero_coverage' in row and row['has_zero_coverage'] else ""
            benchmark_marker = " [BENCHMARK]" if 'is_benchmark' in row and row['is_benchmark'] else ""
            nearest_marker = " [NEAREST]" if 'is_nearest' in row and row['is_nearest'] else ""

            print(f"  {row['contig']}: "
                  f"FC={row['log2_fold_change']:.2f}, "
                  f"TPM={row['positive_tpm']:.1f}/{row['negative_tpm']:.1f}, "
                  f"{row.get('DE_status', 'N/A')}{length_info}{zero_marker}{benchmark_marker}{nearest_marker}")

    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()