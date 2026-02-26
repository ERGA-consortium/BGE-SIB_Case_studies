
rule stage_ref:
    input:
        src = lambda wc: config["refs"][str(wc.hap)]
    output:
        dst = os.path.join(OUTMAIN, "mapping", "input", OUTPRE + ".hap{hap}.ref.fasta")
    resources:
        mem_mb = 1000,
        runtime = 10
    shell: """
        mkdir -p $(dirname {output.dst})
        ln -sfn $(realpath {input.src}) {output.dst}
    """

