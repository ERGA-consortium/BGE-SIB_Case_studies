
def staged_fa(genome, hap):
    return os.path.join(OUTMAIN, "inputs", "genomes", f"{genome}.{hap}.fasta")


rule stage_genomes:
    input:
        expand(
            os.path.join(OUTMAIN, "inputs", "genomes", "{genome}.{hap}.fasta"),
            zip,
            genome=[g for g, h in GENOME_HAPS],
            hap=[h for g, h in GENOME_HAPS],
        )


rule stage_genome_fa:
    input:
        src=lambda wc: config["genomes"][wc.genome][wc.hap]
    output:
        dst=os.path.join(OUTMAIN, "inputs", "genomes", "{genome}.{hap}.fasta")
    resources:
        mem_mb=1000,
        runtime=10,
    shell: """
        mkdir -p $(dirname {output.dst})
        ln -sfn $(realpath {input.src}) {output.dst}
    """

