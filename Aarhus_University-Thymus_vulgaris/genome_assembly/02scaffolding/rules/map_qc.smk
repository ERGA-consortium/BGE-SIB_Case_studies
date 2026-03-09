

import math

DOWNSAMPLE_FRAC = float(config.get("downsample", {}).get("frac", 1))
DOWNSAMPLE_TAG = f"frac{int(DOWNSAMPLE_FRAC * 100):03d}"
TMP_BASE = config.get("tmpdir", "/tmp")


rule do_mapping:
    input:
        expand(
            os.path.join(OUTMAIN, "mapping", "qcs", "multiqc_bam", "multiqc_report_" + OUTPRE + ".hap{hap}_bwa_parse_dedup_split.html"),
            hap=HAPS,
        ),




rule downsample_fastq:
    input:
        fq1 = config["fastq"]["p1"],
        fq2 = config["fastq"]["p2"],
    output:
        fq1 = os.path.join(OUTMAIN, "mapping", "input", f"downsample_{DOWNSAMPLE_TAG}.pair1.fq"),
        fq2 = os.path.join(OUTMAIN, "mapping", "input", f"downsample_{DOWNSAMPLE_TAG}.pair2.fq")
    params:
        frac = DOWNSAMPLE_FRAC
    conda:
        "../envs/seqtk.yml"
    resources:
        mem_mb = 30000,
        runtime = 800,
    shell: """
        seqtk sample -s100 {input.fq1} {params.frac} > {output.fq1}
        seqtk sample -s100 {input.fq2} {params.frac} > {output.fq2}
    """


rule do_chroms_file:
    input:
        fai = lambda wc: REF(wc) + ".fai",
    output:
        chroms_file = os.path.join(OUTMAIN, "mapping", "input", OUTPRE + ".hap{hap}.chroms_file.txt")
    resources:
        mem_mb = 5000,
        runtime = 100,
    shell: """
    cut -f1,2 {input.fai} > {output.chroms_file}
    """
 

rule map_hic_bwa:
    input:
        ref = REF,
        fq1 = os.path.join(OUTMAIN, "mapping", "input", f"downsample_{DOWNSAMPLE_TAG}.pair1.fq") if DOWNSAMPLE_FRAC < 1 else config["fastq"]["p1"],
        fq2 = os.path.join(OUTMAIN, "mapping", "input", f"downsample_{DOWNSAMPLE_TAG}.pair2.fq") if DOWNSAMPLE_FRAC < 1 else config["fastq"]["p2"],
        bwaidx = lambda wc: REF(wc) + ".amb",
        faidx = lambda wc: REF(wc) + ".fai"
    output:
        bam = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa.bam"),
    params:
        bwa_threads = lambda wildcards, threads: max(threads // 2 + 2, 1),
        samtools_threads = lambda wildcards, threads: max(threads // 2 - 2, 1),
    conda:
        "../envs/bwa_pairtools_samtools.yml"
    threads: 16
    resources:
        mem_mb = 100000,
        runtime = 4000,
        cpus_per_task = 24
    shell: """
        bwa mem -5SP -T0 -t {threads} {input.ref} {input.fq1} {input.fq2} | \
            samtools view -@8 -b \
                -o {output.bam}
        """



rule parse_pairtools:
    input:
        bam = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa.bam"),
        chroms_file = os.path.join(OUTMAIN, "mapping", "input", OUTPRE + ".hap{hap}.chroms_file.txt")
    output:
        pairs = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa_parse.pairs"),
    threads: 8
    resources:
        mem_mb = 64000,
        runtime = 2000,
        cpus_per_task = 8
    conda:
        "../envs/bwa_pairtools_samtools.yml"
    shell: """
         pairtools parse \
            --min-mapq 40 \
            --walks-policy 5unique \
            --max-inter-align-gap 30 \
            --nproc-in {threads} \
            --nproc-out {threads} \
            --chroms-path {input.chroms_file} \
            {input.bam} > \
            {output.pairs}
    """



rule sort_dedup_pairtools:
    input:
        pairs = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa_parse.pairs"),
    output:
        pairs = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa_parse_dedup.pairs"),
        stats = os.path.join(OUTMAIN, "mapping", "qcs", OUTPRE + ".hap{hap}_bwa_parse_dedup.stats.txt"),
    params:
        tmp_prefix = lambda wc: os.path.join(TMP_BASE, f"{OUTPRE}.hap{wc.hap}.pairtools"),
    threads: 8
    resources:
        mem_mb = 64000,
        runtime = 2000,
        cpus_per_task = 8
    conda:
        "../envs/bwa_pairtools_samtools.yml"
    shell: """
        JOBTAG=${{SLURM_JOB_ID:-local}}
        TMPDIR="{params.tmp_prefix}.${{JOBTAG}}"
        mkdir -p "${{TMPDIR}}"

        pairtools sort \
            --tmpdir="${{TMPDIR}}" \
            --nproc {threads} \
            {input.pairs} | \
        pairtools dedup \
            --nproc-in {threads} \
            --nproc-out {threads} \
            --mark-dups \
            --output-stats {output.stats} > \
            {output.pairs}
"""
        


rule split_bam_sam_pairtools:
    input:
        pairs = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa_parse_dedup.pairs"),
    output:
        pairs = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa_parse_dedup_split.pairs"),
        bam = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa_parse_dedup_split.bam")
    params:
        pairtools_threads = lambda wildcards, threads: threads // 2,
        samtools_threads = lambda wildcards, threads: threads // 2,
        tmp_prefix = lambda wc: os.path.join(TMP_BASE, f"{OUTPRE}.hap{wc.hap}.split"),
    conda:
        "../envs/bwa_pairtools_samtools.yml"
    threads: 16
    resources:
        mem_mb = 64000,
        runtime = 2000,
        cpus_per_task = 16
    shell: """
        JOBTAG=${{SLURM_JOB_ID:-local}}
        TMPDIR="{params.tmp_prefix}.${{JOBTAG}}"
        mkdir -p "${{TMPDIR}}"

        pairtools split \
            --nproc-in {params.pairtools_threads} \
            --nproc-out {params.pairtools_threads} \
            --output-pairs {output.pairs} \
            --output-sam - \
            {input.pairs} | \
        samtools view -bS -@ {params.samtools_threads} | \
        samtools sort \
            -@ {params.samtools_threads} \
            -T "${{TMPDIR}}/samtools_sort" \
            -o {output.bam}

    """



rule namesort_bam:
    input:
        bam = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa_parse_dedup_split.bam")
    output:
        bam = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa_parse_dedup_split_namesort.bam")
    params:
        tmp_prefix = lambda wc: os.path.join(TMP_BASE, f"{OUTPRE}.hap{wc.hap}.namesort"),
    threads: 16
    resources:
        mem_mb = 64000,
        runtime = 1200,
        cpus_per_task = 16
    conda:
        "../envs/bwa_pairtools_samtools.yml"
    shell: """
        JOBTAG=${{SLURM_JOB_ID:-local}}
        TMPDIR="{params.tmp_prefix}.${{JOBTAG}}"
        mkdir -p "${{TMPDIR}}"

        samtools sort -n -@{threads} \
        -T "${{TMPDIR}}/namesort_tmp" \
        -o {output.bam} \
        {input.bam}
    """



rule do_fastqc_bam:
    input:
        bam = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa_parse_dedup_split.bam")
    output:
        directory(os.path.join(OUTMAIN, "mapping", "qcs", "fastq_bams", OUTPRE + ".hap{hap}_bwa_parse_dedup_split"))
    threads: 3
    resources:
        mem_mb = 50000,
        runtime = 100,
        cpus_per_task = 3
    log: os.path.join(OUTMAIN, "mapping", "qcs", "fastq_bams", OUTPRE + ".hap{hap}_bwa_parse_dedup_split.log")
    conda:
        "../envs/fastqc.yml"
    shell: """
    mkdir -p {output}
    fastqc -o {output} --nogroup -f bam -t {threads} {input.bam} &> {log}
    """




rule flagstat_bam:
    input:
        bam = "{prefix}.bam",
    output:
        stats = "{prefix}.flagstats",
    resources:
        mem_mb = 1000,
        runtime = 300
    conda:
        "../envs/bwa_pairtools_samtools.yml"
    shell:
        "samtools flagstat {input.bam} > {output}"


rule stats_bam:
    input:
        bam = "{prefix}.bam",
    output:
        stats = "{prefix}.stats"
    conda:
        "../envs/bwa_pairtools_samtools.yml"
    resources:
        mem_mb = 1000,
        runtime = 300,
        cpus_per_task = 2
    threads: 2
    shell:
        "samtools stats -@ {threads} {input.bam} > {output}"


rule idxstats_bam:
    input:
        bam = "{prefix}.bam",
        bai = "{prefix}.bam.bai",
    output:
        stats = "{prefix}.idxstats"
    conda:
        "../envs/bwa_pairtools_samtools.yml"
    resources:
        mem_mb = 1000,
        runtime = 300
    shell:
        "samtools idxstats {input.bam} > {output}"



rule do_multiqc_bam:
    input:
        flagstats = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa_parse_dedup_split.flagstats"),
        idxstats = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa_parse_dedup_split.idxstats"),
        stats = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa_parse_dedup_split.stats"),
        fastqc_dir = os.path.join(OUTMAIN, "mapping", "qcs", "fastq_bams", OUTPRE + ".hap{hap}_bwa_parse_dedup_split"),
    output:
        f=os.path.join(OUTMAIN, "mapping", "qcs", "multiqc_bam", "multiqc_report_" + OUTPRE + ".hap{hap}_bwa_parse_dedup_split.html")
    params:
        dirname = lambda wildcards, output: os.path.dirname(output.f),
        base =  lambda wildcards, output: os.path.basename(output.f).replace(".html", "")
    resources:
        mem_mb = 100000,
        runtime = 120
    conda:
        "../envs/multiqc.yml"
    shell: """
        multiqc -f -o {params.dirname} -n {params.base} {input}
    """
