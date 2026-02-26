rule run_braker3:
    input:
        genome=os.path.join(OUTMAIN, "{assembly}", "genome", "{assembly}.softmasked.fa"),
        proteins=os.path.join(OUTMAIN, "protein_evidence", "proteins.clustered.fa"),
        bams=lambda wc: expand(
            os.path.join(OUTMAIN, "{assembly}", "bams", "{library}.{assembly}.bam"),
            assembly=[str(wc.assembly)],
            library=config["rna_libraries"].keys(),
        ),
    output:
        gtf=os.path.join(OUTMAIN, "{assembly}", "braker", "braker.gtf"),
        aa=os.path.join(OUTMAIN, "{assembly}", "braker", "braker.aa"),
        cds=os.path.join(OUTMAIN, "{assembly}", "braker", "braker.codingseq"),
    params:
        genemark_bin=config["braker3"]["genemark_bin"],
        prothint_bin=config["braker3"]["prothint_bin"],
        species=config["species"],
        useexisting="--useexisting" if config["braker3"]["useexisting"] else "",
        extra=config["braker3"]["extra"],
        outdir=os.path.join(OUTMAIN, "{assembly}", "braker"),
        bam_csv=lambda wc, input: ",".join(input.bams),
    threads: 20
    resources:
        mem_mb=500000,
        runtime=2880,
        cpus_per_task=20,
    conda:
        "../envs/braker3.yml"
    shell:
        """
        export PATH="{params.genemark_bin}:$PATH"
        export PATH="{params.prothint_bin}:$PATH"

        braker.pl \
            --genome={input.genome} \
            --prot_seq={input.proteins} \
            --bam={params.bam_csv} \
            --threads {threads} \
            --species="{params.species}" \
            {params.useexisting} \
            --workingdir={params.outdir} \
            {params.extra}
        """
