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


def reverse_complement(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def write_fasta(records, path, width=80):
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for start in range(0, len(seq), width):
                handle.write(seq[start : start + width] + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fa", required=True)
    parser.add_argument("--map-tsv", required=True)
    parser.add_argument("--output-fa", required=True)
    args = parser.parse_args()

    mapping = {}
    with open(args.map_tsv) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            mapping[row["old_id"]] = {
                "new_id": row["new_id"],
                "strand": row.get("strand", "."),
                "order_index": int(row.get("order_index", 10**18)),
            }

    output_records = []
    used_names = set()
    missing_order = 10**18

    for old_id, seq in fasta_records(args.input_fa):
        if old_id in mapping:
            target = mapping[old_id]
            new_id = target["new_id"]
            strand = target["strand"]
            order_index = target["order_index"]
        else:
            new_id = old_id
            strand = "."
            order_index = missing_order
            missing_order += 1

        if strand == "-":
            seq = reverse_complement(seq)

        if new_id in used_names:
            raise ValueError(f"Duplicate output sequence name: {new_id}")
        used_names.add(new_id)
        output_records.append((order_index, new_id, seq))

    output_records.sort(key=lambda item: (item[0], item[1]))
    write_fasta([(name, seq) for _, name, seq in output_records], args.output_fa)


if __name__ == "__main__":
    main()
