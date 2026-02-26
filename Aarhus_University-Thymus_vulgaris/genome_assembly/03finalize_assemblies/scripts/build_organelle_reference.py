#!/usr/bin/env python3
"""Build organelle FASTA with consistent chromosome-style naming."""

from __future__ import annotations

import argparse


def iter_fasta(path: str):
    name = None
    seq_chunks = []
    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks)
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if name is not None:
        yield name, "".join(seq_chunks)


def write_wrapped(handle, sequence: str, width: int = 80) -> None:
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def ensure_unique(name: str, seen_ids: set[str]) -> str:
    if name not in seen_ids:
        seen_ids.add(name)
        return name
    suffix = 2
    candidate = f"{name}x{suffix}"
    while candidate in seen_ids:
        suffix += 1
        candidate = f"{name}x{suffix}"
    seen_ids.add(candidate)
    return candidate


def build_names(
    records: list[tuple[str, str]],
    complete_name: str,
    fragment_prefix: str,
    min_complete_bp: int,
) -> list[tuple[str, str]]:
    if len(records) == 1 and len(records[0][1]) >= min_complete_bp:
        return [(complete_name, records[0][1])]

    pad = max(2, len(str(max(1, len(records)))))
    out = []
    for idx, (_, seq) in enumerate(records, start=1):
        out.append((f"{fragment_prefix}{idx:0{pad}d}", seq))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Combine plastid and mitochondrial FASTAs")
    parser.add_argument("--plastid-fa", required=True)
    parser.add_argument("--mitochondria-fa", required=True)
    parser.add_argument("--output-fa", required=True)
    parser.add_argument("--plastid-complete-name", default="chrC")
    parser.add_argument("--mitochondria-complete-name", default="chrM")
    parser.add_argument("--plastid-fragment-prefix", default="pt")
    parser.add_argument("--mitochondria-fragment-prefix", default="mt")
    parser.add_argument("--plastid-min-complete-bp", type=int, default=0)
    parser.add_argument("--mitochondria-min-complete-bp", type=int, default=0)
    args = parser.parse_args()

    plastid_records = list(iter_fasta(args.plastid_fa))
    mito_records = list(iter_fasta(args.mitochondria_fa))

    named_plastid = build_names(
        records=plastid_records,
        complete_name=args.plastid_complete_name,
        fragment_prefix=args.plastid_fragment_prefix,
        min_complete_bp=args.plastid_min_complete_bp,
    )
    named_mito = build_names(
        records=mito_records,
        complete_name=args.mitochondria_complete_name,
        fragment_prefix=args.mitochondria_fragment_prefix,
        min_complete_bp=args.mitochondria_min_complete_bp,
    )

    seen_ids: set[str] = set()
    with open(args.output_fa, "w") as out_handle:
        for name, seq in named_plastid + named_mito:
            final_name = ensure_unique(name, seen_ids)
            out_handle.write(f">{final_name}\n")
            write_wrapped(out_handle, seq)


if __name__ == "__main__":
    main()
