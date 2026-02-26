#!/usr/bin/env python3
import argparse
import statistics


def count_fasta_records(path):
    n = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                n += 1
    return n


def read_cluster_sizes(path):
    sizes = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            rep = line.split("\t", 1)[0]
            sizes[rep] = sizes.get(rep, 0) + 1
    return sizes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fa", required=True)
    parser.add_argument("--clustered-fa", required=True)
    parser.add_argument("--cluster-tsv", required=True)
    parser.add_argument("--out-tsv", required=True)
    args = parser.parse_args()

    n_input = count_fasta_records(args.input_fa)
    n_rep = count_fasta_records(args.clustered_fa)
    n_removed = n_input - n_rep
    pct_retained = (100.0 * n_rep / n_input) if n_input else 0.0
    pct_removed = (100.0 * n_removed / n_input) if n_input else 0.0

    sizes = read_cluster_sizes(args.cluster_tsv)
    size_values = list(sizes.values())

    n_clusters = len(size_values)
    n_singletons = sum(1 for x in size_values if x == 1)
    mean_cluster_size = statistics.mean(size_values) if size_values else 0.0
    median_cluster_size = statistics.median(size_values) if size_values else 0.0
    max_cluster_size = max(size_values) if size_values else 0

    rows = [
        ("n_input_proteins", n_input),
        ("n_representatives", n_rep),
        ("n_removed_as_redundant", n_removed),
        ("pct_retained", f"{pct_retained:.2f}"),
        ("pct_removed", f"{pct_removed:.2f}"),
        ("n_clusters", n_clusters),
        ("n_singleton_clusters", n_singletons),
        ("mean_cluster_size", f"{mean_cluster_size:.3f}"),
        ("median_cluster_size", f"{median_cluster_size:.3f}"),
        ("max_cluster_size", max_cluster_size),
    ]

    with open(args.out_tsv, "w") as out:
        out.write("metric\tvalue\n")
        for metric, value in rows:
            out.write(f"{metric}\t{value}\n")


if __name__ == "__main__":
    main()
