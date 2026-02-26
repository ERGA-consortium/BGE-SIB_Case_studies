#!/usr/bin/env python3
"""Rewrite AGP object names to finalized IDs and drop removed scaffolds."""

from __future__ import annotations

import argparse
import csv


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


def read_map(path: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            old_id = row.get("old_id")
            new_id = row.get("new_id")
            if old_id and new_id:
                mapping[old_id] = new_id
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser(description="Rewrite AGP object IDs to finalized names")
    parser.add_argument("--input-agp", required=True)
    parser.add_argument("--map-tsv", required=True)
    parser.add_argument("--keep-fa", required=True)
    parser.add_argument("--append-fa", default="")
    parser.add_argument("--output-agp", required=True)
    args = parser.parse_args()

    mapping = read_map(args.map_tsv)
    keep_lengths = read_fasta_lengths(args.keep_fa)
    keep_ids = set(keep_lengths.keys())
    emitted_objects: set[str] = set()

    with open(args.input_agp) as agp_in, open(args.output_agp, "w") as agp_out:
        for raw_line in agp_in:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                agp_out.write(line + "\n")
                continue

            fields = line.split("\t")
            if len(fields) < 9:
                continue

            old_object = fields[0]
            new_object = mapping.get(old_object, old_object)
            fields[0] = new_object

            if new_object not in keep_ids:
                continue

            emitted_objects.add(new_object)
            agp_out.write("\t".join(fields) + "\n")

        if args.append_fa:
            append_lengths = read_fasta_lengths(args.append_fa)
            for seq_id, seq_len in append_lengths.items():
                if seq_id not in keep_ids:
                    continue
                if seq_id in emitted_objects:
                    continue
                if seq_len <= 0:
                    continue
                agp_out.write(
                    f"{seq_id}\t1\t{seq_len}\t1\tW\t{seq_id}\t1\t{seq_len}\t+\n"
                )


if __name__ == "__main__":
    main()
