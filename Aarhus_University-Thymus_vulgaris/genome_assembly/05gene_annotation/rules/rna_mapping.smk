rule star_genome_index:
    input:
        genome=os.path.join(OUTMAIN, "{assembly}", "genome", "{assembly}.softmasked.fa"),
        fai=os.path.join(OUTMAIN, "{assembly}", "genome", "{assembly}.softmasked.fa.fai"),
    output:
        idx=directory(os.path.join(OUTMAIN, "{assembly}", "star", "genomeindex")),
    params:
        nbases=config["star"]["genome_sa_index_nbases"],
    threads: 10
    resources:
        mem_mb=50000,
        runtime=800,
        cpus_per_task=10,
    conda:
        "../envs/star.yml"
    shell:
        """
        STAR --runMode genomeGenerate \
             --genomeSAindexNbases {params.nbases} \
             --genomeDir {output.idx} \
             --genomeFastaFiles {input.genome} \
             --runThreadN {threads}
        """


rule star_pass1:
    input:
        idx=os.path.join(OUTMAIN, "{assembly}", "star", "genomeindex"),
        fq1=lambda wc: config["rna_libraries"][str(wc.library)]["r1"],
        fq2=lambda wc: config["rna_libraries"][str(wc.library)]["r2"],
    output:
        sam=temp(os.path.join(OUTMAIN, "{assembly}", "star", "pass1", "{library}.pass1.Aligned.out.sam")),
        sj=os.path.join(OUTMAIN, "{assembly}", "star", "pass1", "{library}.pass1.SJ.out.tab"),
    params:
        outpre=os.path.join(OUTMAIN, "{assembly}", "star", "pass1", "{library}.pass1."),
        read_cmd=config["star"]["read_files_command"],
        extra=config["star"]["pass1_extra"],
    threads: 6
    resources:
        mem_mb=50000,
        runtime=800,
        cpus_per_task=6,
    conda:
        "../envs/star.yml"
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {input.idx} \
             --readFilesIn {input.fq1} {input.fq2} \
             --readFilesCommand {params.read_cmd} \
             --outFileNamePrefix {params.outpre} \
             {params.extra}
        """


rule star_pass2:
    input:
        idx=os.path.join(OUTMAIN, "{assembly}", "star", "genomeindex"),
        fq1=lambda wc: config["rna_libraries"][str(wc.library)]["r1"],
        fq2=lambda wc: config["rna_libraries"][str(wc.library)]["r2"],
        sj=lambda wc: expand(
            os.path.join(OUTMAIN, "{assembly}", "star", "pass1", "{library}.pass1.SJ.out.tab"),
            assembly=[str(wc.assembly)],
            library=config["rna_libraries"].keys(),
        ),
    output:
        sam=temp(os.path.join(OUTMAIN, "{assembly}", "star", "pass2", "{library}.pass2.Aligned.out.sam")),
    params:
        outpre=os.path.join(OUTMAIN, "{assembly}", "star", "pass2", "{library}.pass2."),
        read_cmd=config["star"]["read_files_command"],
        extra=config["star"]["pass2_extra"],
    threads: 6
    resources:
        mem_mb=50000,
        runtime=800,
        cpus_per_task=6,
    conda:
        "../envs/star.yml"
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {input.idx} \
             --outSAMstrandField intronMotif \
             --sjdbFileChrStartEnd {input.sj} \
             --readFilesIn {input.fq1} {input.fq2} \
             --readFilesCommand {params.read_cmd} \
             --outFileNamePrefix {params.outpre} \
             {params.extra}
        """


rule samtools_sort_bam:
    input:
        sam=os.path.join(OUTMAIN, "{assembly}", "star", "pass2", "{library}.pass2.Aligned.out.sam"),
    output:
        bam=os.path.join(OUTMAIN, "{assembly}", "bams", "{library}.{assembly}.bam"),
    threads: 3
    resources:
        mem_mb=50000,
        runtime=500,
        cpus_per_task=3,
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools sort -@ {threads} -T {output.bam}.tmp -O bam -o {output.bam} {input.sam}
        """


rule index_bam:
    input:
        bam="{prefix}.bam",
    output:
        bai="{prefix}.bam.bai",
    resources:
        mem_mb=1000,
        runtime=300,
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools index {input.bam}"
