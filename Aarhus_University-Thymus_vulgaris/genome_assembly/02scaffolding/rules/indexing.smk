
rule bwa_index_fa:
    input:
        "{refpath}"
    output:
        "{refpath}.amb"
    conda:
        "../envs/bwa_pairtools_samtools.yml"
    resources:
        mem_mb = 5000,
        runtime = 300
    shell:
        "bwa index {input}"


rule samtools_index_fa:
    input:
        "{refpath}"
    output:
        "{refpath}.fai"
    conda:
        "../envs/bwa_pairtools_samtools.yml"
    resources:
        mem_mb = 1000,
        runtime = 300
    shell:
        "samtools faidx {input}"


rule index_bam:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    resources:
        mem_mb = 1000,
        runtime = 300
    conda:
        "../envs/bwa_pairtools_samtools.yml"
    shell:  "samtools index {input}"

