rule stage_softmasked_genome:
    input:
        src=lambda wc: config["assemblies"][str(wc.assembly)]["genome"],
    output:
        dst=os.path.join(OUTMAIN, "{assembly}", "genome", "{assembly}.softmasked.fa"),
    resources:
        mem_mb=1000,
        runtime=30,
    shell:
        """
        ln -sfn $(realpath {input.src}) {output.dst}
        """


rule faidx_softmasked_genome:
    input:
        fa=os.path.join(OUTMAIN, "{assembly}", "genome", "{assembly}.softmasked.fa"),
    output:
        fai=os.path.join(OUTMAIN, "{assembly}", "genome", "{assembly}.softmasked.fa.fai"),
    resources:
        mem_mb=2000,
        runtime=120,
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools faidx {input.fa}
        """
