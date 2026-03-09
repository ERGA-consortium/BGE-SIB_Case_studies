#!/usr/bin/env python3

import argparse
import csv


def fasta_records(path):
    name = None
    seq_chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
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


def write_fasta(records, path, width=80):
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for start in range(0, len(seq), width):
                handle.write(seq[start : start + width] + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fa", required=True)
    parser.add_argument("--output-lg-fa", required=True)
    parser.add_argument("--output-map", required=True)
    parser.add_argument("--output-unplaced", required=True)
    parser.add_argument("--lg-min-bp", type=int, required=True)
    parser.add_argument("--lg-prefix", required=True)
    parser.add_argument("--unplaced-prefix", required=True)
    args = parser.parse_args()

    all_records = list(fasta_records(args.input_fa))
    large = [(name, seq) for name, seq in all_records if len(seq) >= args.lg_min_bp]
    small = [(name, seq) for name, seq in all_records if len(seq) < args.lg_min_bp]

    if not large:
        raise ValueError(
            f"No reference scaffolds >= lg_min_bp ({args.lg_min_bp}). "
            "Lower lg_min_bp or provide a different reference assembly."
        )

    large.sort(key=lambda item: (-len(item[1]), item[0]))
    small.sort(key=lambda item: (-len(item[1]), item[0]))

    un_pad = max(3, len(str(max(1, len(small)))))

    map_rows = []
    lg_records = []
    order_index = 1

    for idx, (old_id, seq) in enumerate(large, start=1):
        new_id = f"{args.lg_prefix}{idx}"
        map_rows.append(
            {
                "old_id": old_id,
                "new_id": new_id,
                "length": len(seq),
                "group": new_id,
                "status": "reference_lg",
                "strand": ".",
                "order_index": order_index,
            }
        )
        lg_records.append((new_id, seq))
        order_index += 1

    unplaced_rows = []
    for idx, (old_id, seq) in enumerate(small, start=1):
        new_id = f"{args.unplaced_prefix}{idx:0{un_pad}d}"
        row = {
            "old_id": old_id,
            "new_id": new_id,
            "length": len(seq),
            "group": "UN",
            "status": "reference_unplaced",
            "strand": ".",
            "order_index": order_index,
        }
        map_rows.append(row)
        unplaced_rows.append(row)
        order_index += 1

    write_fasta(lg_records, args.output_lg_fa)

    with open(args.output_map, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "old_id",
                "new_id",
                "length",
                "group",
                "status",
                "strand",
                "order_index",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(map_rows)

    with open(args.output_unplaced, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["old_id", "new_id", "length", "status"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in unplaced_rows:
            writer.writerow(
                {
                    "old_id": row["old_id"],
                    "new_id": row["new_id"],
                    "length": row["length"],
                    "status": row["status"],
                }
            )


if __name__ == "__main__":
    main()
