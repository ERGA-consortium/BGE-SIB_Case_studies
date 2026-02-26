#!/usr/bin/env python3
"""Filter contigs from FASTA according to a remove-list file."""

from __future__ import annotations

import argparse


def iter_fasta(path: str):
    record_id = None
    seq_chunks = []
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_id is not None:
                    yield record_id, "".join(seq_chunks)
                record_id = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if record_id is not None:
        yield record_id, "".join(seq_chunks)


def write_wrapped(handle, sequence: str, width: int = 80) -> None:
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def read_remove_ids(path: str) -> set[str]:
    remove_ids: set[str] = set()
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            remove_ids.add(line.split()[0])
    return remove_ids


def main() -> None:
    parser = argparse.ArgumentParser(description="Write filtered FASTA")
    parser.add_argument("--input-fa", required=True)
    parser.add_argument("--remove-list", required=True)
    parser.add_argument("--output-fa", required=True)
    parser.add_argument("--mode", required=True, choices=["remove_contigs", "keep_and_mask"])
    args = parser.parse_args()

    remove_ids = read_remove_ids(args.remove_list)
    should_remove = args.mode == "remove_contigs"

    with open(args.output_fa, "w") as out_handle:
        for record_id, sequence in iter_fasta(args.input_fa):
            if should_remove and record_id in remove_ids:
                continue
            out_handle.write(f">{record_id}\n")
            write_wrapped(out_handle, sequence)


if __name__ == "__main__":
    main()
