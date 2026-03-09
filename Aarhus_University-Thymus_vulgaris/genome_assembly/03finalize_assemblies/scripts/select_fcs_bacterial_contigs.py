#!/usr/bin/env python3
"""Select contigs for removal from an FCS-GX report."""

from __future__ import annotations

import argparse


def iter_fasta_lengths(path: str) -> dict[str, int]:
    lengths: dict[str, int] = {}
    record_id = None
    seq_len = 0
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_id is not None:
                    lengths[record_id] = seq_len
                record_id = line[1:].split()[0]
                seq_len = 0
            else:
                seq_len += len(line)
    if record_id is not None:
        lengths[record_id] = seq_len
    return lengths


def norm(s: str) -> str:
    return s.strip().lstrip("#").lower().replace("_", "-")


def find_idx(header: list[str], exact: set[str], contains: str | None = None) -> int:
    normalized = [norm(x) for x in header]
    for i, name in enumerate(normalized):
        if name in exact:
            return i
    if contains is not None:
        for i, name in enumerate(normalized):
            if contains in name:
                return i
    return -1


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract bacterial contaminants from FCS-GX report")
    parser.add_argument("--report", required=True)
    parser.add_argument("--input-fa", required=True)
    parser.add_argument("--out-remove", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--remove-action", action="append", default=["EXCLUDE"])
    args = parser.parse_args()

    actions = {x.upper() for x in args.remove_action}

    header = None
    data_lines: list[str] = []
    with open(args.report) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if header is None and any(norm(x) == "action" for x in fields):
                header = fields
                continue
            if header is not None:
                if line.startswith("#"):
                    continue
                data_lines.append(line)

    if header is None:
        raise ValueError(f"Could not find header with an 'action' column in report: {args.report}")

    seq_idx = find_idx(
        header,
        {"seq-id", "seqid", "seq_id", "accession", "contig", "sequence", "name"},
    )
    action_idx = find_idx(header, {"action"})

    if seq_idx < 0 or action_idx < 0:
        raise ValueError(
            "Could not resolve required columns in FCS report header. "
            f"Header columns: {header}"
        )

    remove_ids: set[str] = set()
    rows_scanned = 0
    rows_selected = 0
    for line in data_lines:
        fields = line.split("\t")
        max_idx = max(seq_idx, action_idx)
        if len(fields) <= max_idx:
            continue
        rows_scanned += 1
        seq_id = fields[seq_idx].strip()
        action = fields[action_idx].strip().upper()
        if action in actions:
            remove_ids.add(seq_id)
            rows_selected += 1

    with open(args.out_remove, "w") as out_handle:
        for seq_id in sorted(remove_ids):
            out_handle.write(seq_id + "\n")

    lengths = iter_fasta_lengths(args.input_fa)
    removed_bp = sum(lengths.get(seq_id, 0) for seq_id in remove_ids)
    missing = sum(1 for seq_id in remove_ids if seq_id not in lengths)

    with open(args.out_summary, "w") as out_handle:
        out_handle.write("metric\tvalue\n")
        out_handle.write(f"rows_scanned\t{rows_scanned}\n")
        out_handle.write(f"rows_selected\t{rows_selected}\n")
        out_handle.write(f"contigs_selected\t{len(remove_ids)}\n")
        out_handle.write(f"contigs_missing_in_fasta\t{missing}\n")
        out_handle.write(f"bp_selected\t{removed_bp}\n")
        out_handle.write(f"remove_actions\t{','.join(sorted(actions))}\n")


if __name__ == "__main__":
    main()
