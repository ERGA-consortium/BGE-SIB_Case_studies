#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path


def parse_agat_stats(path: Path) -> dict[str, str]:
    out = {}
    patterns = {
        "genes": re.compile(r"number of genes?\s*[:=]\s*([0-9]+)", re.I),
        "mrna": re.compile(r"number of mrna\s*[:=]\s*([0-9]+)", re.I),
        "single_exon_genes": re.compile(r"number of single exon genes?\s*[:=]\s*([0-9]+)", re.I),
        "single_exon_transcripts": re.compile(r"number of single exon transcripts?\s*[:=]\s*([0-9]+)", re.I),
        "monoexonic_transcripts": re.compile(r"number of monoexonic transcripts?\s*[:=]\s*([0-9]+)", re.I),
    }

    for line in path.read_text().splitlines():
        for key, pat in patterns.items():
            if key in out:
                continue
            m = pat.search(line.strip())
            if m:
                out[key] = m.group(1)
    return out


def parse_busco_summary(busco_dir: Path) -> dict[str, str]:
    out = {}
    files = sorted(busco_dir.rglob("short_summary*.txt"))
    if not files:
        return out

    text = files[0].read_text()
    out["busco_summary_file"] = str(files[0])

    m = re.search(
        r"C:([0-9.]+)%\[S:([0-9.]+)%,D:([0-9.]+)%\],F:([0-9.]+)%,M:([0-9.]+)%,n:([0-9]+)",
        text,
    )
    if m:
        out["busco_C_percent"] = m.group(1)
        out["busco_S_percent"] = m.group(2)
        out["busco_D_percent"] = m.group(3)
        out["busco_F_percent"] = m.group(4)
        out["busco_M_percent"] = m.group(5)
        out["busco_n"] = m.group(6)
        return out

    for line in text.splitlines():
        if "Complete BUSCOs" in line:
            val = re.search(r"([0-9]+)", line)
            if val:
                out["busco_complete_count"] = val.group(1)
        elif "Duplicated BUSCOs" in line:
            val = re.search(r"([0-9]+)", line)
            if val:
                out["busco_duplicated_count"] = val.group(1)
        elif "Fragmented BUSCOs" in line:
            val = re.search(r"([0-9]+)", line)
            if val:
                out["busco_fragmented_count"] = val.group(1)
        elif "Missing BUSCOs" in line:
            val = re.search(r"([0-9]+)", line)
            if val:
                out["busco_missing_count"] = val.group(1)
    return out


def parse_omark_summary(omark_dir: Path) -> dict[str, str]:
    out = {}
    files = sorted(omark_dir.rglob("*.sum"))
    if not files:
        return out
    out["omark_sum_file"] = str(files[0])

    text = files[0].read_text()
    key_patterns = {
        "omark_single": re.compile(r"\bsingle\b[^0-9]*([0-9.]+)", re.I),
        "omark_duplicated": re.compile(r"\bduplicated\b[^0-9]*([0-9.]+)", re.I),
        "omark_missing": re.compile(r"\bmissing\b[^0-9]*([0-9.]+)", re.I),
        "omark_consistent": re.compile(r"\bconsistent\b[^0-9]*([0-9.]+)", re.I),
        "omark_inconsistent": re.compile(r"\binconsistent\b[^0-9]*([0-9.]+)", re.I),
        "omark_contaminants": re.compile(r"\bcontaminants?\b[^0-9]*([0-9.]+)", re.I),
    }
    for key, pat in key_patterns.items():
        m = pat.search(text)
        if m:
            out[key] = m.group(1)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-stats", required=True)
    parser.add_argument("--longest-stats", required=True)
    parser.add_argument("--busco-dir", required=True)
    parser.add_argument("--omark-dir", required=True)
    parser.add_argument("--out-tsv", required=True)
    args = parser.parse_args()

    metrics = {}

    full = parse_agat_stats(Path(args.full_stats))
    for k, v in full.items():
        metrics[f"agat_full_{k}"] = v

    longest = parse_agat_stats(Path(args.longest_stats))
    for k, v in longest.items():
        metrics[f"agat_longest_{k}"] = v

    busco = parse_busco_summary(Path(args.busco_dir))
    for k, v in busco.items():
        metrics[k] = v

    omark = parse_omark_summary(Path(args.omark_dir))
    for k, v in omark.items():
        metrics[k] = v

    with open(args.out_tsv, "w") as out:
        out.write("metric\tvalue\n")
        for key in sorted(metrics):
            out.write(f"{key}\t{metrics[key]}\n")


if __name__ == "__main__":
    main()
