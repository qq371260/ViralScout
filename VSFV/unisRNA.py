import pysam
import argparse
import os
import sys


def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Filter BAM file to remove multi-mapped reads")
    parser.add_argument('-i', '--input_bam', required=True,
                        help="Path to the input BAM file")
    parser.add_argument('-o', '--output_bam', default='unisRNA.csv',
                        help="Path to the output filtered BAM file")
    parser.add_argument('-c', '--cores', type=int, default=os.cpu_count(),
                        help=f"Number of CPU cores to use (default: {os.cpu_count()} - all available cores)")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Input validation
    if not os.path.exists(args.input_bam):
        print(f"Error: Input BAM file '{args.input_bam}' does not exist", file=sys.stderr)
        sys.exit(1)

    try:
        # Open the input BAM file with multi-threading support using context manager
        with pysam.AlignmentFile(args.input_bam, "rb", threads=args.cores) as samfile:
            # Open the output BAM file with the same header
            with pysam.AlignmentFile(args.output_bam, "wb", header=samfile.header, threads=args.cores) as output:

                # Create a set to track reads with secondary alignments
                marked_reads = set()

                print("First pass: Identifying multi-mapped reads...", file=sys.stderr)

                # First pass: Identify reads with secondary alignments
                # These represent multi-mapped reads that we want to filter out
                for read in samfile.fetch(until_eof=True):
                    if read.is_secondary:  # Check for secondary alignment flag
                        marked_reads.add(read.query_name)  # Mark read by its name

                print(f"Found {len(marked_reads)} multi-mapped reads to filter", file=sys.stderr)

                # Reset the file pointer to the beginning of the BAM file for second pass
                samfile.reset()

                print("Second pass: Writing filtered BAM...", file=sys.stderr)
                written_count = 0
                total_count = 0

                # Second pass: Write only the reads that are not in the 'marked_reads' set
                # This removes all alignments (both primary and secondary) of multi-mapped reads
                for read in samfile.fetch(until_eof=True):
                    total_count += 1
                    if read.query_name not in marked_reads:
                        output.write(read)
                        written_count += 1

                # Print detailed statistics
                filtered_count = total_count - written_count
                print(f"Processed {total_count} total reads", file=sys.stderr)
                print(f"Written {written_count} unique reads ({written_count / total_count * 100:.1f}%)",
                      file=sys.stderr)
                print(f"Filtered {filtered_count} reads ({filtered_count / total_count * 100:.1f}%)", file=sys.stderr)

        print(f"Filtered BAM file successfully saved as {args.output_bam}")

    except Exception as e:
        print(f"Error processing BAM files: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()