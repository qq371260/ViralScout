#!/usr/bin/env python3
import sys
import gzip
import argparse
from collections import Counter


def parse_lengths(spec):
    """Parse length specification, support comma-separated list or range"""
    if not spec:
        return []

    lengths = set()

    # Split comma-separated parts
    for part in spec.split(','):
        part = part.strip()
        if '-' in part:
            # Handle range, e.g., "23-24"
            try:
                start, end = map(int, part.split('-'))
                lengths.update(range(start, end + 1))
            except ValueError:
                print(f"Warning: Invalid range '{part}'", file=sys.stderr)
        else:
            # Handle single number, e.g., "23"
            try:
                lengths.add(int(part))
            except ValueError:
                print(f"Warning: Invalid length '{part}'", file=sys.stderr)

    return sorted(lengths)


def calculate_ratio(fastq_file, indicator_lengths=None):
    """Calculate proportion of specified sRNA lengths"""
    if indicator_lengths is None:
        indicator_lengths = [23, 24]  # Default value

    counts = Counter()
    total = 0

    # Determine if file is gzipped
    is_gzipped = fastq_file.endswith('.gz')

    try:
        opener = gzip.open if is_gzipped else open
        mode = 'rt' if is_gzipped else 'r'

        with opener(fastq_file, mode) as fh:
            for i, line in enumerate(fh):
                if (i % 4) == 1:  # Sequence line
                    total += 1
                    counts[len(line.strip())] += 1

        # Calculate sum of specified lengths
        indicator_sum = sum(counts.get(length, 0) for length in indicator_lengths)
        ratio = indicator_sum / total if total > 0 else 0.0

        return ratio, total, indicator_sum

    except Exception as e:
        print(f"Error processing {fastq_file}: {e}", file=sys.stderr)
        return 0.0, 0, 0


def main():
    parser = argparse.ArgumentParser(
        description='Calculate proportion of indicator sRNA lengths in a FASTQ file'
    )

    parser.add_argument(
        'fastq_file',
        help='Input FASTQ file (can be gzipped)'
    )

    parser.add_argument(
        '--indicator_sRNA',
        default='23,24',
        help='Comma-separated lengths or range (e.g., "23,24" or "23-24"). Default: 23,24'
    )

    args = parser.parse_args()

    # Parse length parameter
    indicator_lengths = parse_lengths(args.indicator_sRNA)

    if not indicator_lengths:
        print(f"Error: No valid lengths in '{args.indicator_sRNA}'", file=sys.stderr)
        sys.exit(1)

    # Calculate proportion
    ratio, total, indicator_sum = calculate_ratio(args.fastq_file, indicator_lengths)

    # Output proportion (easy for bash parsing)
    print(f"{ratio:.6f}")

    # Write result file
    lengths_str = args.indicator_sRNA
    with open("sRNA_indicator_proportion.txt", "w") as f:
        f.write(f"Indicator lengths: {lengths_str}nt\n")
        f.write(f"Total reads: {total}\n")
        f.write(f"Indicator reads: {indicator_sum}\n")
        f.write(f"Proportion: {ratio:.6f}\n")
        f.write(f"23–24-nt proportion < 0.3 indicates high virus-plant sRNA pattern convergence\n")
        f.write(f"No reference yet for other virus-host systems")

if __name__ == "__main__":
    main()