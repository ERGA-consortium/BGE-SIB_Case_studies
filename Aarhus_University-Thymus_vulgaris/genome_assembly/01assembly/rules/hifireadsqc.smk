
##############################################
######## TARGET RULES TO RUN PIPELINE ########
##############################################

rule all_hifi_qcs:
    input:
        os.path.join(OUTMAIN, "hifi_reads", "stats", PRENAME + "_raw.stats"),
        os.path.join(OUTMAIN, "hifi_reads", "fastqc", PRENAME + "_raw"),
        os.path.join(OUTMAIN, "hifi_reads", "stats", PRENAME + "_hifis_cleaned_seqkit_cutadapt.stats"),
        os.path.join(OUTMAIN, "hifi_reads", "fastqc", PRENAME + "_hifis_cleaned_seqkit_cutadapt"),
        os.path.join(OUTMAIN, "hifi_reads", "nanoplot", "raw_reads"),
        os.path.join(OUTMAIN, "hifi_reads", "nanoplot", "clean_reads")


rule pretrim_hifi_qcs:
    input:
        os.path.join(OUTMAIN, "hifi_reads", "stats", PRENAME + "_raw.stats"),
        os.path.join(OUTMAIN, "hifi_reads", "fastqc", PRENAME + "_raw"),



###################################
######## RULES TO DO STUFF ########
###################################

rule hifi_stats_raw:
    input:
        fq = config["infq"]
    output:
        stats = os.path.join(OUTMAIN, "hifi_reads", "stats", PRENAME + "_raw.stats")
    threads: 3
    resources:
        mem_mb = 3000,
        runtime = 300,
        cpus_per_task = 3
    conda:
        "../envs/seqkit.yml"
    shell: """
    seqkit stats -a --fq-encoding "sanger" --threads {threads} --out-file {output.stats} {input.fq}
"""


rule fastqc_hifi_raw:
    input:
        fq = config["infq"],
    output:
        directory(os.path.join(OUTMAIN, "hifi_reads", "fastqc", PRENAME + "_raw"))
    threads: 3
    resources:
        mem_mb = 5000,
        runtime = 500,
        cpus_per_task = 3
    log: os.path.join(OUTMAIN, "hifi_reads", "fastqc", PRENAME + "_raw.log")
    conda:
        "../envs/fastqc.yml"
    shell: """
    mkdir -p {output}
    fastqc -o {output} -f fastq -t {threads} {input.fq} &> {log}
    """


rule nanoplot_hifi_raw:
    input:
        fq = config["infq"],
    output:
        outdir = directory(os.path.join(OUTMAIN, "hifi_reads", "nanoplot", "raw_reads"))
    threads: 3
    params:
        prefix = PRENAME + "_rawreads"
    resources:
        mem_mb = 5000,
        runtime = 500,
        cpus_per_task = 3
    conda:
        "../envs/nanoplot.yml"
    shell: """
        NanoPlot --threads {threads} --fastq {input.fq} \
        --outdir {output.outdir} --prefix {params.prefix} \
        --plots dot hex kde --format png --tsv_stats --N50
    """
    
 


rule seqkit_filter:
# remove reads with too low average quality or too short
    input:
        fq = config["infq"]
    output:
        fq = os.path.join(OUTMAIN, "hifi_reads", "filtered_reads", PRENAME + "_hifis_cleaned_seqkit.fq.gz")
    log: os.path.join(OUTMAIN, "hifi_reads", "filtered_reads", PRENAME + "_hifis_cleaned_seqkit.log")
    params:
        minqual = config["cutadapt_params"]["minqual"],
        minlen = config["cutadapt_params"]["minlen"],
    threads: 10
    resources:
        mem_mb = 10000,
        runtime = 600,
        cpus_per_task = 10
    conda:
        "../envs/seqkit.yml"
    shell: """
        seqkit seq --min-qual {params.minqual} --min-len {params.minlen} \
         --out-file {output.fq} {input.fq}
    """




rule cutadapt:
# remove reads with adapters, also remove too short reads
# parameter setting based on VGP pipeline galaxy tutorial:
# https://training.galaxyproject.org/training-material/topics/assembly/tutorials/vgp_genome_assembly/tutorial.html#generation-of-k-mer-spectra-with-meryl
    input:
        fq = os.path.join(OUTMAIN, "hifi_reads", "filtered_reads", PRENAME + "_hifis_cleaned_seqkit.fq.gz")
    output:
        fq = os.path.join(OUTMAIN, "hifi_reads", "filtered_reads", PRENAME + "_hifis_cleaned_seqkit_cutadapt.fq.gz")
    log: os.path.join(OUTMAIN, "hifi_reads", "filtered_reads", PRENAME + "_hifis_cleaned_seqkit_cutadapt.log")
    params:
        error = config["cutadapt_params"]["error"],
        minoverlap = config["cutadapt_params"]["minoverlap"],
        minlen = config["cutadapt_params"]["minlen"],
        adapter1=config["cutadapt_params"]["adapter1"],
        adapter2=config["cutadapt_params"]["adapter2"]
    threads: 10
    resources:
        mem_mb = 10000,
        runtime = 600,
        cpus_per_task = 10
    conda:
        "../envs/cutadapt.yml"
    shell: """
        cutadapt -j {threads} \
         -b {params.adapter1} -b {params.adapter2} \
        --error-rate {params.error} --overlap {params.minoverlap} \
        --minimum-length {params.minlen} \
        --discard-trimmed \
        -o {output.fq} {input.fq} &> {log}
    """


rule hifi_stats_cleaned:
    input:
        fq = os.path.join(OUTMAIN, "hifi_reads", "filtered_reads", PRENAME + "_hifis_cleaned_seqkit_cutadapt.fq.gz")
    output:
        stats = os.path.join(OUTMAIN, "hifi_reads", "stats", PRENAME + "_hifis_cleaned_seqkit_cutadapt.stats")
    threads: 3
    resources:
        mem_mb = 3000,
        runtime = 300,
        cpus_per_task = 3
    conda:
        "../envs/seqkit.yml"
    shell: """
    seqkit stats -a --fq-encoding "sanger" --threads {threads} --out-file {output.stats} {input.fq}
"""


rule fastqc_hifi_cleaned:
    input:
        fq = os.path.join(OUTMAIN, "hifi_reads", "filtered_reads", PRENAME + "_hifis_cleaned_seqkit_cutadapt.fq.gz")
    output:
        directory(os.path.join(OUTMAIN, "hifi_reads", "fastqc", PRENAME + "_hifis_cleaned_seqkit_cutadapt"))
    threads: 3
    resources:
        mem_mb = 5000,
        runtime = 500,
        cpus_per_task = 3
    log: os.path.join(OUTMAIN, "hifi_reads", "fastqc", PRENAME + "_hifis_cleaned_seqkit_cutadapt.log")
    conda:
        "../envs/fastqc.yml"
    shell: """
    mkdir -p {output}
    fastqc -o {output} -f fastq -t {threads} {input.fq} &> {log}
    """


rule nanoplot_hifi_cleaned:
    input:
        fq = os.path.join(OUTMAIN, "hifi_reads", "filtered_reads", PRENAME + "_hifis_cleaned_seqkit_cutadapt.fq.gz")
    output:
        outdir = directory(os.path.join(OUTMAIN, "hifi_reads", "nanoplot", "clean_reads"))
    threads: 3
    params:
        prefix = PRENAME + "_hifis_cleaned_seqkit_cutadapt"
    resources:
        mem_mb = 5000,
        runtime = 500,
        cpus_per_task = 3
    conda:
        "../envs/nanoplot.yml"
    shell: """
        NanoPlot --threads {threads} --fastq {input.fq} \
        --outdir {output.outdir} --prefix {params.prefix} \
        --plots dot hex kde --format png --tsv_stats --N50
    """
    