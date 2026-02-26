rule fastqc_bam:
    input:
        bam=os.path.join(OUTMAIN, "{assembly}", "bams", "{library}.{assembly}.bam"),
    output:
        d=directory(os.path.join(OUTMAIN, "{assembly}", "qc", "fastqc_bams", "{library}.{assembly}")),
    threads: 3
    resources:
        mem_mb=50000,
        runtime=100,
        cpus_per_task=3,
    log:
        os.path.join(OUTMAIN, "{assembly}", "qc", "fastqc_bams", "{library}.{assembly}.log")
    conda:
        "../envs/fastqc.yml"
    shell:
        """
	mkdir -p {output.d}
        fastqc -o {output.d} --nogroup -f bam -t {threads} {input.bam} &> {log}
        """


rule flagstat_bam:
    input:
        bam="{prefix}.bam",
        bai="{prefix}.bam.bai",
    output:
        txt="{prefix}.flagstats",
    resources:
        mem_mb=1000,
        runtime=300,
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools flagstat {input.bam} > {output.txt}"


rule stats_bam:
    input:
        bam="{prefix}.bam",
        bai="{prefix}.bam.bai",
    output:
        txt="{prefix}.stats",
    threads: 2
    resources:
        mem_mb=1000,
        runtime=300,
        cpus_per_task=2,
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools stats -@ {threads} {input.bam} > {output.txt}"


rule idxstats_bam:
    input:
        bam="{prefix}.bam",
        bai="{prefix}.bam.bai",
    output:
        txt="{prefix}.idxstats",
    resources:
        mem_mb=1000,
        runtime=300,
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools idxstats {input.bam} > {output.txt}"


rule multiqc_bam:
    input:
        lambda wc: (
            expand(
                os.path.join(OUTMAIN, "{assembly}", "bams", "{library}.{assembly}.{ext}"),
                assembly=[str(wc.assembly)],
                library=config["rna_libraries"].keys(),
                ext=["flagstats", "idxstats", "stats"],
            )
            + expand(
                os.path.join(OUTMAIN, "{assembly}", "qc", "fastqc_bams", "{library}.{assembly}"),
                assembly=[str(wc.assembly)],
                library=config["rna_libraries"].keys(),
            )
        ),
    output:
        html=os.path.join(OUTMAIN, "{assembly}", "qc", "multiqc", "multiqc_report.html"),
    params:
        outdir=os.path.join(OUTMAIN, "{assembly}", "qc", "multiqc"),
        name="multiqc_report.html",
    resources:
        mem_mb=20000,
        runtime=120,
    conda:
        "../envs/multiqc.yml"
    shell:
        """
        multiqc -f -o {params.outdir} -n {params.name} {input}
        """
