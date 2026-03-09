#!/usr/bin/env python3

import argparse
import csv
import re
from collections import Counter, defaultdict


def merged_bp(intervals):
    if not intervals:
        return 0
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return sum(end - start for start, end in merged)


def fasta_lengths(path):
    lengths = {}
    order = []
    name = None
    length = 0
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                    order.append(name)
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
    if name is not None:
        lengths[name] = length
        order.append(name)
    return lengths, order


def read_paf(path):
    stats = defaultdict(
        lambda: defaultdict(
            lambda: {
                "bp": 0,
                "nmatch": 0,
                "alen": 0,
                "plus_bp": 0,
                "minus_bp": 0,
                "min_tstart": 10**18,
                "intervals": [],
            }
        )
    )

    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue

            qname = cols[0]
            qstart = int(cols[2])
            qend = int(cols[3])
            strand = cols[4]
            tname = cols[5]
            tstart = int(cols[7])
            nmatch = int(cols[9])
            alen = int(cols[10])

            bp = max(0, qend - qstart)
            entry = stats[qname][tname]
            entry["bp"] += bp
            entry["nmatch"] += nmatch
            entry["alen"] += alen
            if strand == "+":
                entry["plus_bp"] += bp
            else:
                entry["minus_bp"] += bp
            if tstart < entry["min_tstart"]:
                entry["min_tstart"] = tstart
            if bp > 0:
                entry["intervals"].append((qstart, qend))
    return stats


def group_sort_key(name):
    match = re.match(r"^([A-Za-z]+)(\d+)$", name)
    if match:
        return (match.group(1), int(match.group(2)))
    return (name, 10**18)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf", required=True)
    parser.add_argument("--query-fa", required=True)
    parser.add_argument("--out-map", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--min-query-bp", type=int, required=True)
    parser.add_argument("--unplaced-prefix", required=True)
    args = parser.parse_args()

    query_lengths, query_order = fasta_lengths(args.query_fa)
    paf_stats = read_paf(args.paf)

    assigned_by_group = defaultdict(list)
    unassigned_rows = []

    for query_id in query_order:
        qlen = query_lengths[query_id]
        if qlen < args.min_query_bp:
            unassigned_rows.append(
                {
                    "old_id": query_id,
                    "length": qlen,
                    "group": args.unplaced_prefix,
                    "status": "below_min_query_bp",
                    "strand": ".",
                    "best_target": ".",
                    "second_target": ".",
                    "best_bp": 0,
                    "second_bp": 0,
                    "best_ratio": "inf",
                    "best_identity": 0.0,
                    "query_covered": 0.0,
                    "order_key": 10**18,
                }
            )
            continue
        target_stats = paf_stats.get(query_id, {})

        if not target_stats:
            unassigned_rows.append(
                {
                    "old_id": query_id,
                    "length": qlen,
                    "group": args.unplaced_prefix,
                    "status": "no_alignment",
                    "strand": ".",
                        "best_target": ".",
                        "second_target": ".",
                        "best_bp": 0,
                        "second_bp": 0,
                        "best_ratio": "inf",
                    "best_identity": 0.0,
                    "query_covered": 0.0,
                    "order_key": 10**18,
                }
            )
            continue

        ranked = sorted(
            target_stats.items(),
            key=lambda item: (
                merged_bp(item[1]["intervals"]),
                item[1]["bp"],
            ),
            reverse=True,
        )
        best_target, best = ranked[0]
        second_target = ranked[1][0] if len(ranked) > 1 else "."
        second_bp = merged_bp(ranked[1][1]["intervals"]) if len(ranked) > 1 else 0
        best_bp = merged_bp(best["intervals"])
        best_ratio = float("inf") if second_bp == 0 else best_bp / second_bp
        best_identity = (best["nmatch"] / best["alen"]) if best["alen"] > 0 else 0.0
        query_covered = (best_bp / qlen) if qlen > 0 else 0.0
        strand = "+" if best["plus_bp"] >= best["minus_bp"] else "-"
        status = "assigned_lg"

        row = {
            "old_id": query_id,
            "length": qlen,
            "group": best_target,
            "status": status,
            "strand": strand,
            "best_target": best_target,
            "second_target": second_target,
            "best_bp": best_bp,
            "second_bp": second_bp,
            "best_ratio": best_ratio,
            "best_identity": best_identity,
            "query_covered": query_covered,
            "order_key": best["min_tstart"],
        }

        assigned_by_group[best_target].append(row)

    final_rows = []
    order_index = 1

    for group_name in sorted(assigned_by_group.keys(), key=group_sort_key):
        rows = assigned_by_group[group_name]
        rows.sort(key=lambda row: (row["order_key"], -row["best_bp"], row["old_id"]))
        if len(rows) == 1:
            rows[0]["new_id"] = group_name
        else:
            for idx, row in enumerate(rows, start=1):
                row["new_id"] = f"{group_name}p{idx:03d}"
        for row in rows:
            row["order_index"] = order_index
            final_rows.append(row)
            order_index += 1

    unassigned_rows.sort(key=lambda row: (-row["length"], row["old_id"]))
    un_pad = max(3, len(str(max(1, len(unassigned_rows)))))
    for idx, row in enumerate(unassigned_rows, start=1):
        row["new_id"] = f"{args.unplaced_prefix}{idx:0{un_pad}d}"
        row["order_index"] = order_index
        final_rows.append(row)
        order_index += 1

    with open(args.out_map, "w", newline="") as handle:
        for row in final_rows:
            row.pop("order_key", None)
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "old_id",
                "new_id",
                "length",
                "group",
                "status",
                "strand",
                "best_target",
                "second_target",
                "best_bp",
                "second_bp",
                "best_ratio",
                "best_identity",
                "query_covered",
                "order_index",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(final_rows)

    status_counts = Counter()
    status_bp = Counter()
    group_counts = Counter()
    group_bp = Counter()
    for row in final_rows:
        status_counts[row["status"]] += 1
        status_bp[row["status"]] += int(row["length"])
        group_counts[row["group"]] += 1
        group_bp[row["group"]] += int(row["length"])

    with open(args.out_summary, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["kind", "key", "count", "total_bp"],
            delimiter="\t",
        )
        writer.writeheader()
        for key in sorted(status_counts.keys()):
            writer.writerow(
                {
                    "kind": "status",
                    "key": key,
                    "count": status_counts[key],
                    "total_bp": status_bp[key],
                }
            )
        for key in sorted(group_counts.keys()):
            writer.writerow(
                {
                    "kind": "group",
                    "key": key,
                    "count": group_counts[key],
                    "total_bp": group_bp[key],
                }
            )


if __name__ == "__main__":
    main()
