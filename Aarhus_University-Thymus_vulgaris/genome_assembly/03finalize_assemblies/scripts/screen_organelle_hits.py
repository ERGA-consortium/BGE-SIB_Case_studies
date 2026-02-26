#!/usr/bin/env python3
"""Classify organelle-like contigs from minimap2 PAF alignments."""

from __future__ import annotations

import argparse
from collections import defaultdict


def read_fasta_lengths(path: str) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current_id = None
    current_len = 0
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)
    if current_id is not None:
        lengths[current_id] = current_len
    return lengths


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
            continue
        merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def organelle_label(
    target_name: str,
    plastid_complete_name: str,
    mitochondria_complete_name: str,
    plastid_fragment_prefix: str,
    mitochondria_fragment_prefix: str,
) -> str:
    if target_name == plastid_complete_name or target_name.startswith(plastid_fragment_prefix):
        return "plastid"
    if target_name == mitochondria_complete_name or target_name.startswith(mitochondria_fragment_prefix):
        return "mitochondria"
    return "OTHER"


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize organelle-like contigs from PAF")
    parser.add_argument("--paf", required=True)
    parser.add_argument("--query-fa", required=True)
    parser.add_argument("--out-remove", required=True)
    parser.add_argument("--out-bed", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--min-identity", required=True, type=float)
    parser.add_argument("--min-aln-bp", required=True, type=int)
    parser.add_argument("--min-mapq", required=True, type=int)
    parser.add_argument("--min-cov-frac", required=True, type=float)
    parser.add_argument("--min-cov-bp", required=True, type=int)
    parser.add_argument("--plastid-complete-name", default="chrC")
    parser.add_argument("--mitochondria-complete-name", default="chrM")
    parser.add_argument("--plastid-fragment-prefix", default="pt")
    parser.add_argument("--mitochondria-fragment-prefix", default="mt")
    args = parser.parse_args()

    query_lengths = read_fasta_lengths(args.query_fa)

    hit_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    hit_target_bp: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    filtered_hit_counts: dict[str, int] = defaultdict(int)

    with open(args.paf) as paf_handle:
        for raw_line in paf_handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue

            query = fields[0]
            if query not in query_lengths:
                continue

            q_start = int(fields[2])
            q_end = int(fields[3])
            target = fields[5]
            matches = int(fields[9])
            aln_len = int(fields[10])
            mapq = int(fields[11])

            if aln_len <= 0:
                continue
            identity = matches / aln_len

            if identity < args.min_identity:
                continue
            if aln_len < args.min_aln_bp:
                continue
            if mapq < args.min_mapq:
                continue
            if q_end <= q_start:
                continue

            hit_intervals[query].append((q_start, q_end))
            hit_target_bp[query][target] += aln_len
            filtered_hit_counts[query] += 1

    remove_contigs: list[str] = []

    with open(args.out_summary, "w") as summary_handle, open(args.out_bed, "w") as bed_handle:
        summary_handle.write(
            "contig\tcontig_bp\torganelle_aligned_bp\torganelle_aligned_frac\tfiltered_hits\ttop_target\ttop_target_bp\ttop_target_type\tremove\n"
        )

        for contig, contig_bp in query_lengths.items():
            intervals = merge_intervals(hit_intervals.get(contig, []))
            aligned_bp = sum(end - start for start, end in intervals)
            aligned_frac = aligned_bp / contig_bp if contig_bp else 0.0
            target_bp = hit_target_bp.get(contig, {})

            if target_bp:
                top_target = max(target_bp.items(), key=lambda item: item[1])[0]
                top_target_bp = target_bp[top_target]
                top_type = organelle_label(
                    top_target,
                    args.plastid_complete_name,
                    args.mitochondria_complete_name,
                    args.plastid_fragment_prefix,
                    args.mitochondria_fragment_prefix,
                )
            else:
                top_target = "."
                top_target_bp = 0
                top_type = "."

            remove = aligned_frac >= args.min_cov_frac and aligned_bp >= args.min_cov_bp
            if remove:
                remove_contigs.append(contig)

            summary_handle.write(
                f"{contig}\t{contig_bp}\t{aligned_bp}\t{aligned_frac:.6f}\t{filtered_hit_counts.get(contig, 0)}"
                f"\t{top_target}\t{top_target_bp}\t{top_type}\t{int(remove)}\n"
            )

            for start, end in intervals:
                bed_handle.write(f"{contig}\t{start}\t{end}\torganelle_like\t{top_type}\n")

    with open(args.out_remove, "w") as remove_handle:
        for contig in remove_contigs:
            remove_handle.write(contig + "\n")


if __name__ == "__main__":
    main()
