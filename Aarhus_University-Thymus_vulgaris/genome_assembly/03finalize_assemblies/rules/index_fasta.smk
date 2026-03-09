rule do_index_all_fasta:
    input:
        expand(
            os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap1_plus_organelle.fa.fai"),
            sample=config["samples"].keys(),
        )
        + expand(
            os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap2_nuclear.fa.fai"),
            sample=config["samples"].keys(),
        ),


rule index_fasta_fai:
    input:
        fa="{fa}",
    output:
        fai="{fa}.fai",
    conda:
        "../envs/samtools.yml"
    resources:
        mem_mb=2000,
        runtime=120,
    shell:
        "samtools faidx {input.fa}"
