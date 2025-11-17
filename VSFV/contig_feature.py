#!/usr/bin/env python3
"""
Extract contig metrics (size, coverage, GC content) for filtered contigs
and merge with existing filtered_contigs.csv
"""

import pandas as pd
import pysam
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
import argparse
import os
import sys
import glob
from collections import defaultdict
import logging
import re

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Extract contig metrics (size, coverage, GC content) and merge with filtered contigs',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-f', '--fasta', required=True,
                        help='Reference genome FASTA file')
    parser.add_argument('-b', '--bam', required=True,
                        help='Sorted BAM file with mapped reads (must be indexed)')
    parser.add_argument('-c', '--contigs', required=True,
                        help='Input filtered_contigs.csv file from correlation analysis')
    parser.add_argument('-o', '--output', required=True,
                        help='Output enhanced CSV file with contig metrics')
    parser.add_argument('--min-coverage', type=float, default=0.1,
                        help='Minimum coverage threshold for reporting')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of threads for BAM processing')
    parser.add_argument('--coverage-bins', type=int, default=1000,
                        help='Number of bins for coverage calculation in large contigs')
    parser.add_argument('--benchmarks-file', type=str, default=None,
                        help='Custom benchmarks file name (default: auto-detect any benchmark*.txt file)')
    parser.add_argument('--min-length', type=int, default=18,
                        help='Minimum sRNA size')
    parser.add_argument('--max-length', type=int, default=30,
                        help='Maximum sRNA size')
    return parser.parse_args()


def detect_correlation_method(df):
    """Detect whether the data uses Spearman or Pearson correlation"""
    if 'Mean_Spearman_Correlation' in df.columns:
        return 'spearman'
    elif 'Mean_Pearson_Correlation' in df.columns:
        return 'pearson'
    else:
        # Try to detect from other column names or filename patterns
        for col in df.columns:
            if 'spearman' in col.lower():
                return 'spearman'
            elif 'pearson' in col.lower():
                return 'pearson'
        # Default to spearman if cannot detect
        return 'spearman'


def extract_first_word(name):
    """Extract only the first continuous word from the name (including - and _) and remove _minus suffix"""
    first_word = re.split(r'[\s\t,;]', str(name).strip())[0]

    # Remove _minus suffix if present
    if first_word.endswith('_minus'):
        first_word = first_word[:-6]  # Remove the last 6 characters '_minus'

    return first_word


def extract_contig_sequences(fasta_file, contig_names, output_fasta):
    """
    Extract sequences for filtered contigs to a separate FASTA file
    Matches contig names by the first word (including - and _ connected)
    """
    logger.info(f"Extracting contig sequences from {fasta_file} to {output_fasta}")

    # Create a mapping from first word to full header and sequence
    contig_name_set = set(contig_names)
    extracted_count = 0

    with open(output_fasta, 'w') as out_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Extract first word from the FASTA header
            first_word = extract_first_word(record.id)

            if first_word in contig_name_set:
                # Write the complete record with original header
                SeqIO.write(record, out_file, "fasta")
                extracted_count += 1
                logger.debug(f"Extracted sequence for: {record.id} (matched: {first_word})")

    logger.info(f"Extracted {extracted_count} contig sequences")

    # Check for missing contigs
    found_first_words = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        found_first_words.add(extract_first_word(record.id))

    missing_contigs = contig_name_set - found_first_words
    if missing_contigs:
        logger.warning(f"Missing {len(missing_contigs)} contigs in FASTA file")
        for contig in list(missing_contigs)[:5]:
            logger.warning(f"  - {contig}")
        if len(missing_contigs) > 5:
            logger.warning(f"  ... and {len(missing_contigs) - 5} more")

    return extracted_count


def load_filtered_contigs(contigs_file):
    """Load filtered contigs from CSV file"""
    logger.info(f"Loading filtered contigs from: {contigs_file}")
    df = pd.read_csv(contigs_file)

    # Detect correlation method
    correlation_method = detect_correlation_method(df)
    logger.info(f"Detected correlation method: {correlation_method}")

    # Set required columns based on detected method
    if correlation_method == 'pearson':
        required_columns = ['Contig_Name', 'Mean_Pearson_Correlation', 'P_Value']
        correlation_column = 'Mean_Pearson_Correlation'
    else:
        required_columns = ['Contig_Name', 'Mean_Spearman_Correlation', 'P_Value']
        correlation_column = 'Mean_Spearman_Correlation'

    for col in required_columns:
        if col not in df.columns:
            logger.error(f"Required column '{col}' not found in {contigs_file}")
            sys.exit(1)

    # Create a mapping from processed name to original name
    df['Processed_Name'] = df['Contig_Name'].apply(extract_first_word)

    logger.info(f"Loaded {len(df)} contigs")
    return df, correlation_method, correlation_column


def load_benchmarks(benchmarks_file):
    """Load benchmark contigs from text file (one contig name per line)"""
    logger.info(f"Loading benchmark contigs from: {benchmarks_file}")

    if not os.path.exists(benchmarks_file):
        logger.warning(f"Benchmarks file not found: {benchmarks_file}")
        return []

    with open(benchmarks_file, 'r') as f:
        benchmarks = [extract_first_word(line.strip()) for line in f if line.strip()]

    logger.info(f"Loaded {len(benchmarks)} benchmark contigs")
    return benchmarks


def get_contig_sizes_and_gc(fasta_file, contig_df):
    """Extract contig sizes and GC content from FASTA file"""
    logger.info(f"Reading contig sizes and GC content from: {fasta_file}")

    contig_metrics = {}
    found_contigs = set()

    # Create mapping from processed name to original name
    processed_to_original = dict(zip(contig_df['Processed_Name'], contig_df['Contig_Name']))
    processed_names = set(contig_df['Processed_Name'])

    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract first word for matching
        first_word = extract_first_word(record.id)

        if first_word in processed_names:
            size = len(record.seq)
            gc_content = GC(record.seq)  # Convert to percentage

            # Store with original contig name
            original_name = processed_to_original[first_word]
            contig_metrics[original_name] = {
                'Size_bp': size,
                'GC_Content_percent': gc_content
            }
            found_contigs.add(original_name)

    # Check for missing contigs
    missing_contigs = set(contig_df['Contig_Name']) - found_contigs
    if missing_contigs:
        logger.warning(f"Missing {len(missing_contigs)} contigs in FASTA file")
        for contig in list(missing_contigs)[:5]:  # Show first 5 missing
            logger.warning(f"  - {contig}")
        if len(missing_contigs) > 5:
            logger.warning(f"  ... and {len(missing_contigs) - 5} more")

    logger.info(f"Found metrics for {len(found_contigs)} contigs in FASTA")
    return contig_metrics


def calculate_sRNA_length_distribution(bam_file, processed_names, min_length=18, max_length=30, threads=1):
    """
    Calculate sRNA length distribution for each contig
    Returns dictionary with total reads, 21nt percentage, and 22nt percentage
    """
    logger.info(f"Calculating sRNA length distribution from BAM file: {bam_file}")
    logger.info(f"Length range: {min_length}-{max_length} nt")

    # Check if BAM file is indexed
    if not os.path.exists(bam_file + '.bai'):
        logger.error(f"BAM index not found: {bam_file}.bai")
        logger.error("Please index the BAM file with: samtools index your_file.bam")
        sys.exit(1)

    sRNA_data = {}
    bam = pysam.AlignmentFile(bam_file, "rb", threads=threads)

    # Create mapping from first word to full contig names in BAM
    bam_contig_map = {}
    for seq in bam.header['SQ']:
        full_name = seq['SN']
        first_word = extract_first_word(full_name)
        bam_contig_map[first_word] = full_name

    processed = 0
    total_contigs = len(processed_names)

    for contig_first_word in processed_names:
        processed += 1
        if processed % 100 == 0:
            logger.info(f"Processed {processed}/{total_contigs} contigs for sRNA length distribution")

        try:
            # Get full contig name from BAM
            if contig_first_word not in bam_contig_map:
                logger.warning(f"Contig {contig_first_word} not found in BAM header, skipping sRNA analysis")
                sRNA_data[contig_first_word] = {
                    'total_sRNA_reads': 0,
                    'sRNA_21nt_percent': 0.0,
                    'sRNA_22nt_percent': 0.0,
                    'sRNA_21_22nt_percent': 0.0
                }
                continue

            full_contig_name = bam_contig_map[contig_first_word]

            # Initialize counters
            total_reads = 0
            count_21nt = 0
            count_22nt = 0

            # Iterate through reads mapped to this contig
            for read in bam.fetch(full_contig_name):
                read_length = read.query_length

                # Only count reads within the specified length range
                if min_length <= read_length <= max_length:
                    total_reads += 1
                    if read_length == 21:
                        count_21nt += 1
                    elif read_length == 22:
                        count_22nt += 1

            # Calculate percentages
            sRNA_21nt_percent = (count_21nt / total_reads * 100) if total_reads > 0 else 0.0
            sRNA_22nt_percent = (count_22nt / total_reads * 100) if total_reads > 0 else 0.0
            sRNA_21_22nt_percent = sRNA_21nt_percent + sRNA_22nt_percent

            sRNA_data[contig_first_word] = {
                'total_sRNA_reads': total_reads,
                'sRNA_21nt_percent': sRNA_21nt_percent,
                'sRNA_22nt_percent': sRNA_22nt_percent,
                'sRNA_21_22nt_percent': sRNA_21_22nt_percent
            }

        except ValueError as e:
            logger.warning(f"Error calculating sRNA distribution for {contig_first_word}: {e}")
            sRNA_data[contig_first_word] = {
                'total_sRNA_reads': 0,
                'sRNA_21nt_percent': 0.0,
                'sRNA_22nt_percent': 0.0,
                'sRNA_21_22nt_percent': 0.0
            }
        except Exception as e:
            logger.warning(f"Unexpected error for {contig_first_word}: {e}")
            sRNA_data[contig_first_word] = {
                'total_sRNA_reads': 0,
                'sRNA_21nt_percent': 0.0,
                'sRNA_22nt_percent': 0.0,
                'sRNA_21_22nt_percent': 0.0
            }

    bam.close()

    # Log summary statistics
    total_sRNA_reads = sum(data['total_sRNA_reads'] for data in sRNA_data.values())
    avg_21nt_percent = np.mean(
        [data['sRNA_21nt_percent'] for data in sRNA_data.values() if data['total_sRNA_reads'] > 0])
    avg_22nt_percent = np.mean(
        [data['sRNA_22nt_percent'] for data in sRNA_data.values() if data['total_sRNA_reads'] > 0])

    logger.info(f"sRNA length distribution analysis completed for {len(sRNA_data)} contigs")
    logger.info(f"Total sRNA reads: {total_sRNA_reads}")
    logger.info(f"Average 21nt percentage: {avg_21nt_percent:.2f}%")
    logger.info(f"Average 22nt percentage: {avg_22nt_percent:.2f}%")

    return sRNA_data


def calculate_coverage(bam_file, processed_names, min_coverage=0.1, bins=1000, threads=1):
    """Calculate average coverage for contigs from BAM file"""
    logger.info(f"Calculating coverage from BAM file: {bam_file}")

    # Check if BAM file is indexed
    if not os.path.exists(bam_file + '.bai'):
        logger.error(f"BAM index not found: {bam_file}.bai")
        logger.error("Please index the BAM file with: samtools index your_file.bam")
        sys.exit(1)

    coverage_data = {}
    bam = pysam.AlignmentFile(bam_file, "rb", threads=threads)

    # Create mapping from first word to full contig names in BAM
    bam_contig_map = {}
    for seq in bam.header['SQ']:
        full_name = seq['SN']
        first_word = extract_first_word(full_name)
        bam_contig_map[first_word] = full_name

    processed = 0
    total_contigs = len(processed_names)

    for contig_first_word in processed_names:
        processed += 1
        if processed % 100 == 0:
            logger.info(f"Processed {processed}/{total_contigs} contigs for coverage")

        try:
            # Get full contig name from BAM
            if contig_first_word not in bam_contig_map:
                logger.warning(f"Contig {contig_first_word} not found in BAM header, skipping coverage calculation")
                coverage_data[contig_first_word] = 0.0
                continue

            full_contig_name = bam_contig_map[contig_first_word]

            # Get contig length from BAM header
            contig_length = None
            for seq in bam.header['SQ']:
                if seq['SN'] == full_contig_name:
                    contig_length = seq['LN']
                    break

            if contig_length is None:
                logger.warning(f"Contig {full_contig_name} not found in BAM header, skipping coverage calculation")
                coverage_data[contig_first_word] = 0.0
                continue

            # For large contigs, use binning to save memory
            if contig_length > 100000:  # Large contig, use binning
                bin_size = max(1, contig_length // bins)
                coverage_values = []

                for start in range(0, contig_length, bin_size):
                    end = min(start + bin_size, contig_length)
                    # Use samtools depth-like approach with pysam
                    depth_count = 0
                    bases_covered = 0

                    for pileupcolumn in bam.pileup(full_contig_name, start, end, stepper='nofilter'):
                        if start <= pileupcolumn.pos < end:
                            depth_count += pileupcolumn.n
                            bases_covered += 1

                    if bases_covered > 0:
                        coverage_values.append(depth_count / bases_covered)
                    else:
                        coverage_values.append(0)

                avg_coverage = np.mean(coverage_values) if coverage_values else 0
            else:
                # For small contigs, calculate coverage directly
                coverage_values = []
                for pileupcolumn in bam.pileup(full_contig_name, stepper='nofilter'):
                    coverage_values.append(pileupcolumn.n)

                avg_coverage = np.mean(coverage_values) if coverage_values else 0

            # Apply minimum coverage threshold
            if avg_coverage < min_coverage:
                avg_coverage = 0.0

            coverage_data[contig_first_word] = avg_coverage

        except ValueError as e:
            logger.warning(f"Error calculating coverage for {contig_first_word}: {e}")
            coverage_data[contig_first_word] = 0.0
        except Exception as e:
            logger.warning(f"Unexpected error for {contig_first_word}: {e}")
            coverage_data[contig_first_word] = 0.0

    bam.close()
    logger.info(f"Coverage calculation completed for {len(coverage_data)} contigs")
    return coverage_data


def merge_metrics(filtered_df, size_gc_metrics, coverage_data, sRNA_data, correlation_column):
    """Merge all metrics into the filtered contigs DataFrame"""
    logger.info("Merging all contig metrics...")

    # Create new columns
    sizes = []
    gc_contents = []
    coverages = []
    total_sRNA_reads = []
    sRNA_21nt_percents = []
    sRNA_22nt_percents = []
    sRNA_21_22nt_percents = []

    for original_name, processed_name in zip(filtered_df['Contig_Name'], filtered_df['Processed_Name']):
        if original_name in size_gc_metrics:
            sizes.append(size_gc_metrics[original_name]['Size_bp'])
            gc_contents.append(size_gc_metrics[original_name]['GC_Content_percent'])
        else:
            sizes.append(0)
            gc_contents.append(0.0)

        if processed_name in coverage_data:
            coverages.append(coverage_data[processed_name])
        else:
            coverages.append(0.0)

        if processed_name in sRNA_data:
            total_sRNA_reads.append(sRNA_data[processed_name]['total_sRNA_reads'])
            sRNA_21nt_percents.append(sRNA_data[processed_name]['sRNA_21nt_percent'])
            sRNA_22nt_percents.append(sRNA_data[processed_name]['sRNA_22nt_percent'])
            sRNA_21_22nt_percents.append(sRNA_data[processed_name]['sRNA_21_22nt_percent'])
        else:
            total_sRNA_reads.append(0)
            sRNA_21nt_percents.append(0.0)
            sRNA_22nt_percents.append(0.0)
            sRNA_21_22nt_percents.append(0.0)

    # Add new columns to DataFrame
    enhanced_df = filtered_df.copy()
    enhanced_df['Size_bp'] = sizes
    enhanced_df['GC_Content_percent'] = gc_contents
    enhanced_df['Avg_Coverage'] = coverages
    enhanced_df['Total_sRNA_Reads'] = total_sRNA_reads
    enhanced_df['sRNA_21nt_Percent'] = sRNA_21nt_percents
    enhanced_df['sRNA_22nt_Percent'] = sRNA_22nt_percents
    enhanced_df['sRNA_21_22nt_Percent'] = sRNA_21_22nt_percents

    # 删除临时列
    if 'Processed_Name' in enhanced_df.columns:
        enhanced_df = enhanced_df.drop('Processed_Name', axis=1)

    # Reorder columns for better readability
    column_order = [
        'Contig_Name', 'Size_bp', 'Avg_Coverage', 'GC_Content_percent',
        'Total_sRNA_Reads', 'sRNA_21nt_Percent', 'sRNA_22nt_Percent', 'sRNA_21_22nt_Percent',
        correlation_column, 'P_Value'  # Use dynamic correlation column name
    ]

    # Keep any additional columns that might be in the original file
    extra_columns = [col for col in enhanced_df.columns if col not in column_order]
    final_column_order = column_order + extra_columns

    return enhanced_df[final_column_order]


def generate_benchmarks_features(fasta_file, bam_file, benchmarks, min_coverage=0.1, bins=1000, threads=1,
                                 min_length=18, max_length=30):
    """Generate feature metrics for benchmark contigs"""
    if not benchmarks:
        logger.warning("No benchmark contigs provided")
        return None

    logger.info("Generating feature metrics for benchmark contigs...")

    # Get size and GC content
    # 为benchmarks创建一个临时的DataFrame结构
    benchmarks_df_temp = pd.DataFrame({'Contig_Name': benchmarks, 'Processed_Name': benchmarks})
    size_gc_metrics = get_contig_sizes_and_gc(fasta_file, benchmarks_df_temp)

    # Calculate coverage
    coverage_data = calculate_coverage(bam_file, benchmarks, min_coverage, bins, threads)

    # Calculate sRNA length distribution
    sRNA_data = calculate_sRNA_length_distribution(bam_file, benchmarks, min_length, max_length, threads)

    # Create DataFrame
    benchmarks_data = []
    for contig in benchmarks:
        if contig in size_gc_metrics and contig in coverage_data and contig in sRNA_data:
            benchmarks_data.append({
                'Contig_Name': contig,
                'Size_bp': size_gc_metrics[contig]['Size_bp'],
                'Avg_Coverage': coverage_data[contig],
                'GC_Content_percent': size_gc_metrics[contig]['GC_Content_percent'],
                'Total_sRNA_Reads': sRNA_data[contig]['total_sRNA_reads'],
                'sRNA_21nt_Percent': sRNA_data[contig]['sRNA_21nt_percent'],
                'sRNA_22nt_Percent': sRNA_data[contig]['sRNA_22nt_percent'],
                'sRNA_21_22nt_Percent': sRNA_data[contig]['sRNA_21_22nt_percent'],
                'Type': 'Benchmark'
            })
        elif contig in size_gc_metrics:
            benchmarks_data.append({
                'Contig_Name': contig,
                'Size_bp': size_gc_metrics[contig]['Size_bp'],
                'Avg_Coverage': 0.0,
                'GC_Content_percent': size_gc_metrics[contig]['GC_Content_percent'],
                'Total_sRNA_Reads': sRNA_data.get(contig, {}).get('total_sRNA_reads', 0),
                'sRNA_21nt_Percent': sRNA_data.get(contig, {}).get('sRNA_21nt_percent', 0.0),
                'sRNA_22nt_Percent': sRNA_data.get(contig, {}).get('sRNA_22nt_percent', 0.0),
                'sRNA_21_22nt_Percent': sRNA_data.get(contig, {}).get('sRNA_21_22nt_percent', 0.0),
                'Type': 'Benchmark'
            })

    if benchmarks_data:
        benchmarks_df = pd.DataFrame(benchmarks_data)
        logger.info(f"Generated feature metrics for {len(benchmarks_df)} benchmark contigs")
        return benchmarks_df
    else:
        logger.warning("No benchmark contigs found with valid metrics")
        return None


def generate_summary(enhanced_df, correlation_column, benchmarks_df=None, benchmarks_file=None):
    """Generate summary statistics"""
    logger.info("Generating summary statistics...")

    summary = {
        'total_contigs': len(enhanced_df),
        'contigs_with_coverage': len(enhanced_df[enhanced_df['Avg_Coverage'] > 0]),
        'contigs_with_size': len(enhanced_df[enhanced_df['Size_bp'] > 0]),
        'contigs_with_sRNA': len(enhanced_df[enhanced_df['Total_sRNA_Reads'] > 0]),
        'avg_size': enhanced_df['Size_bp'].mean(),
        'median_size': enhanced_df['Size_bp'].median(),
        'avg_coverage': enhanced_df['Avg_Coverage'].mean(),
        'median_coverage': enhanced_df['Avg_Coverage'].median(),
        'avg_gc': enhanced_df['GC_Content_percent'].mean(),
        'median_gc': enhanced_df['GC_Content_percent'].median(),
        'avg_total_sRNA': enhanced_df['Total_sRNA_Reads'].mean(),
        'median_total_sRNA': enhanced_df['Total_sRNA_Reads'].median(),
        'avg_21nt_percent': enhanced_df['sRNA_21nt_Percent'].mean(),
        'median_21nt_percent': enhanced_df['sRNA_21nt_Percent'].median(),
        'avg_22nt_percent': enhanced_df['sRNA_22nt_Percent'].mean(),
        'median_22nt_percent': enhanced_df['sRNA_22nt_Percent'].median(),
        'avg_21_22nt_percent': enhanced_df['sRNA_21_22nt_Percent'].mean(),
        'median_21_22nt_percent': enhanced_df['sRNA_21_22nt_Percent'].median(),
        'avg_correlation': enhanced_df[correlation_column].mean(),  # Use dynamic column name
        'correlation_method': 'Pearson' if 'pearson' in correlation_column.lower() else 'Spearman'  # Record method
    }

    if benchmarks_df is not None:
        summary['benchmark_contigs'] = len(benchmarks_df)
        summary['benchmark_avg_size'] = benchmarks_df['Size_bp'].mean()
        summary['benchmark_avg_coverage'] = benchmarks_df['Avg_Coverage'].mean()
        summary['benchmark_avg_gc'] = benchmarks_df['GC_Content_percent'].mean()
        summary['benchmark_avg_total_sRNA'] = benchmarks_df['Total_sRNA_Reads'].mean()
        summary['benchmark_avg_21nt_percent'] = benchmarks_df['sRNA_21nt_Percent'].mean()
        summary['benchmark_avg_22nt_percent'] = benchmarks_df['sRNA_22nt_Percent'].mean()
        summary['benchmark_avg_21_22nt_percent'] = benchmarks_df['sRNA_21_22nt_Percent'].mean()
        summary['benchmark_file'] = benchmarks_file

    return summary


def main():
    """Main function"""
    args = parse_arguments()

    # Check input files exist
    for file_path in [args.fasta, args.bam, args.contigs]:
        if not os.path.exists(file_path):
            logger.error(f"Input file not found: {file_path}")
            sys.exit(1)

    logger.info("Starting contig feature analysis")
    logger.info(f"FASTA file: {args.fasta}")
    logger.info(f"BAM file: {args.bam}")
    logger.info(f"Contigs file: {args.contigs}")
    logger.info(f"Output file: {args.output}")
    logger.info(f"sRNA size range: {args.min_length}-{args.max_length} nt")

    # Step 1: Load filtered contigs
    filtered_df, correlation_method, correlation_column = load_filtered_contigs(args.contigs)

    # Use processed names for matching, but keep original names in output
    processed_names = filtered_df['Processed_Name'].tolist()
    original_names = filtered_df['Contig_Name'].tolist()

    # Step 2: Extract sequences to FASTA file
    output_fasta = args.output.replace('.csv', '_sequences.fasta')
    extracted_count = extract_contig_sequences(args.fasta, processed_names, output_fasta)
    logger.info(f"Extracted {extracted_count} contig sequences to {output_fasta}")

    # Step 3: Extract size and GC content from FASTA
    size_gc_metrics = get_contig_sizes_and_gc(args.fasta, filtered_df)

    # Step 4: Calculate sRNA size distribution
    sRNA_data = calculate_sRNA_length_distribution(
        args.bam, processed_names,
        min_length=args.min_length,
        max_length=args.max_length,
        threads=args.threads
    )

    # Step 5: Calculate coverage from BAM
    coverage_data = calculate_coverage(
        args.bam, processed_names,
        min_coverage=args.min_coverage,
        bins=args.coverage_bins,
        threads=args.threads
    )

    # Step 6: Merge all metrics (使用原始名称)
    enhanced_df = merge_metrics(filtered_df, size_gc_metrics, coverage_data, sRNA_data, correlation_column)

    # Step 7: Auto-detect and process benchmarks file
    benchmarks_df = None
    benchmarks_file_used = None

    # Check for benchmarks files
    if args.benchmarks_file:
        # Use user-specified benchmarks file - search in multiple locations
        possible_locations = [
            args.benchmarks_file,
            f"../{args.benchmarks_file}",
            f"./{args.benchmarks_file}"
        ]

        found_benchmark_file = None
        for location in possible_locations:
            if os.path.exists(location):
                found_benchmark_file = location
                break

        if found_benchmark_file:
            benchmarks_file_used = found_benchmark_file
            logger.info(f"Using specified benchmarks file: {benchmarks_file_used}")
            benchmarks = load_benchmarks(benchmarks_file_used)
            benchmarks_df = generate_benchmarks_features(
                args.fasta, args.bam, benchmarks,
                min_coverage=args.min_coverage,
                bins=args.coverage_bins,
                threads=args.threads,
                min_length=args.min_length,
                max_length=args.max_length
            )
        else:
            logger.warning(f"Specified benchmarks file not found in any location: {args.benchmarks_file}")
            logger.warning("Checked locations:")
            for location in possible_locations:
                logger.warning(f"  - {location}")

    # Save benchmarks features if available
    if benchmarks_df is not None:
        # Create output filename based on input benchmarks file
        if benchmarks_file_used:
            base_name = os.path.splitext(os.path.basename(benchmarks_file_used))[0]
            benchmarks_output = f"{base_name}_features.csv"
        else:
            benchmarks_output = "benchmark_features.csv"

        benchmarks_df.to_csv(benchmarks_output, index=False)
        logger.info(f"Benchmark features saved to: {benchmarks_output}")

    # Step 8: Generate summary
    summary = generate_summary(enhanced_df, correlation_column, benchmarks_df, benchmarks_file_used)

    # Step 9: Save results
    enhanced_df.to_csv(args.output, index=False)
    logger.info(f"Enhanced contig data saved to: {args.output}")

    # Save summary as well
    summary_file = args.output.replace('.csv', '_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("Contig Feature Analysis Summary\n")
        f.write("=" * 50 + "\n")
        f.write(f"Correlation method: {summary['correlation_method']}\n")
        f.write(f"Total contigs analyzed: {summary['total_contigs']}\n")
        f.write(f"Contigs with coverage data: {summary['contigs_with_coverage']}\n")
        f.write(f"Contigs with size data: {summary['contigs_with_size']}\n")
        f.write(f"Contigs with sRNA data: {summary['contigs_with_sRNA']}\n")
        f.write(f"Average contig size: {summary['avg_size']:.0f} bp\n")
        f.write(f"Median contig size: {summary['median_size']:.0f} bp\n")
        f.write(f"Average coverage: {summary['avg_coverage']:.2f}\n")
        f.write(f"Median coverage: {summary['median_coverage']:.2f}\n")
        f.write(f"Average GC content: {summary['avg_gc']:.2f}%\n")
        f.write(f"Median GC content: {summary['median_gc']:.2f}%\n")
        f.write(f"Average total sRNA reads: {summary['avg_total_sRNA']:.1f}\n")
        f.write(f"Median total sRNA reads: {summary['median_total_sRNA']:.1f}\n")
        f.write(f"Average 21nt percentage: {summary['avg_21nt_percent']:.2f}%\n")
        f.write(f"Median 21nt percentage: {summary['median_21nt_percent']:.2f}%\n")
        f.write(f"Average 22nt percentage: {summary['avg_22nt_percent']:.2f}%\n")
        f.write(f"Median 22nt percentage: {summary['median_22nt_percent']:.2f}%\n")
        f.write(f"Average 21+22nt percentage: {summary['avg_21_22nt_percent']:.2f}%\n")
        f.write(f"Median 21+22nt percentage: {summary['median_21_22nt_percent']:.2f}%\n")
        f.write(f"Average {summary['correlation_method']} correlation: {summary['avg_correlation']:.4f}\n")
        f.write(f"Contig sequences extracted: {extracted_count}\n")

        if benchmarks_df is not None:
            f.write("\nBenchmark Contigs Summary:\n")
            f.write(f"Benchmark file: {summary.get('benchmark_file', 'Unknown')}\n")
            f.write(f"Total benchmark contigs: {summary['benchmark_contigs']}\n")
            f.write(f"Average benchmark size: {summary['benchmark_avg_size']:.0f} bp\n")
            f.write(f"Average benchmark coverage: {summary['benchmark_avg_coverage']:.2f}\n")
            f.write(f"Average benchmark GC content: {summary['benchmark_avg_gc']:.2f}%\n")
            f.write(f"Average benchmark total sRNA: {summary['benchmark_avg_total_sRNA']:.1f}\n")
            f.write(f"Average benchmark 21nt percentage: {summary['benchmark_avg_21nt_percent']:.2f}%\n")
            f.write(f"Average benchmark 22nt percentage: {summary['benchmark_avg_22nt_percent']:.2f}%\n")
            f.write(f"Average benchmark 21+22nt percentage: {summary['benchmark_avg_21_22nt_percent']:.2f}%\n")

    logger.info(f"Summary saved to: {summary_file}")

    # Print quick summary to console
    print("\n" + "=" * 60)
    print("QUICK SUMMARY")
    print("=" * 60)
    print(f"Correlation method: {summary['correlation_method']}")
    print(f"Total contigs: {summary['total_contigs']}")
    print(f"Average size: {summary['avg_size']:.0f} bp")
    print(f"Average coverage: {summary['avg_coverage']:.2f}")
    print(f"Average GC content: {summary['avg_gc']:.2f}%")
    print(f"Average total sRNA reads: {summary['avg_total_sRNA']:.1f}")
    print(f"Average 21nt percentage: {summary['avg_21nt_percent']:.2f}%")
    print(f"Average 22nt percentage: {summary['avg_22nt_percent']:.2f}%")
    print(f"Average 21+22nt percentage: {summary['avg_21_22nt_percent']:.2f}%")
    print(f"Average correlation: {summary['avg_correlation']:.4f}")
    print(f"Contigs with zero coverage: {summary['total_contigs'] - summary['contigs_with_coverage']}")
    print(f"Contigs with zero sRNA: {summary['total_contigs'] - summary['contigs_with_sRNA']}")
    print(f"Contig sequences extracted: {extracted_count}")

    if benchmarks_df is not None:
        print(f"Benchmark contigs analyzed: {summary['benchmark_contigs']}")
        print(f"Benchmark file: {summary.get('benchmark_file', 'Unknown')}")

    print("=" * 60)

    # Show top 10 contigs by correlation
    print(f"\nTop 10 contigs by {summary['correlation_method']} correlation:")
    print("-" * 100)
    top_contigs = enhanced_df.nlargest(10, correlation_column)[
        ['Contig_Name', 'Size_bp', 'Avg_Coverage', 'GC_Content_percent',
         'Total_sRNA_Reads', 'sRNA_21nt_Percent', 'sRNA_22nt_Percent', correlation_column]
    ]
    for idx, row in top_contigs.iterrows():
        print(f"{row['Contig_Name']:20} | Size: {row['Size_bp']:8} bp | "
              f"Cov: {row['Avg_Coverage']:6.2f} | GC: {row['GC_Content_percent']:5.1f}% | "
              f"sRNA: {row['Total_sRNA_Reads']:6} | "
              f"21nt: {row['sRNA_21nt_Percent']:5.1f}% | 22nt: {row['sRNA_22nt_Percent']:5.1f}% | "
              f"Corr: {row[correlation_column]:.4f}")

    # Show benchmark contigs if available
    if benchmarks_df is not None:
        print("\nBenchmark contigs:")
        print("-" * 100)
        for idx, row in benchmarks_df.iterrows():
            print(f"{row['Contig_Name']:20} | Size: {row['Size_bp']:8} bp | "
                  f"Cov: {row['Avg_Coverage']:6.2f} | GC: {row['GC_Content_percent']:5.1f}% | "
                  f"sRNA: {row['Total_sRNA_Reads']:6} | "
                  f"21nt: {row['sRNA_21nt_Percent']:5.1f}% | 22nt: {row['sRNA_22nt_Percent']:5.1f}%")


if __name__ == "__main__":
    main()