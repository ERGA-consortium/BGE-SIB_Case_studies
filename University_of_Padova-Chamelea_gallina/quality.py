#!/usr/bin/env python3

import argparse
from gzip import open as gzopen

def calc_q_stats(fastq_file, q_thresholds=(20, 30), max_reads=None):
    total_reads = 0
    counts = {q: 0 for q in q_thresholds}

    open_func = gzopen if fastq_file.endswith(".gz") else open

    with open_func(fastq_file, 'rt') as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()

            phred_scores = [ord(c) - 33 for c in qual.strip()]
            avg_q = sum(phred_scores) / len(phred_scores)

            for q in q_thresholds:
                if avg_q > q:
                    counts[q] += 1

            total_reads += 1
            if max_reads and total_reads >= max_reads:
                break

    print(f"\n📊 Total reads processed: {total_reads}")
    for q in sorted(q_thresholds):
        pct = 100 * counts[q] / total_reads if total_reads else 0
        print(f"Reads with avg Q > {q}: {counts[q]} ({pct:.2f}%)")

def main():
    parser = argparse.ArgumentParser(description="Calculate read-level average quality statistics from a FASTQ file.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ or FASTQ.gz file")
    parser.add_argument("--q-thresholds", nargs='+', type=int, default=[20, 30], help="Quality thresholds to evaluate (default: 20 30)")
    parser.add_argument("--max-reads", type=int, default=None, help="Limit to first N reads (optional)")

    args = parser.parse_args()
    calc_q_stats(args.input, q_thresholds=args.q_thresholds, max_reads=args.max_reads)

if __name__ == "__main__":
    main()
