

rule do_hifiasm:
    input:
        expand(os.path.join(OUTMAIN, "assembly", "hifiasm", PRENAME + ".{filt}.asm.bp.p_ctg.gfa"), filt=FILTERS),
        expand(os.path.join(OUTMAIN, "assembly", "hifiasm", PRENAME + ".{filt}.asm.bp.hap1.p_ctg.gfa"), filt=FILTERS),
        expand(os.path.join(OUTMAIN, "assembly", "hifiasm", PRENAME + ".{filt}.asm.bp.hap2.p_ctg.gfa"), filt=FILTERS)


rule hifiasm_assembly:
    """run hifiasm with filter-specific parameter sets."""
    input:
        fq = os.path.join(OUTMAIN, "hifi_reads", "filtered_reads", PRENAME + "_hifis_cleaned_seqkit_cutadapt.fq.gz")
    output:
        gfa_primary = os.path.join(OUTMAIN, "assembly", "hifiasm", PRENAME + ".{filt}.asm.bp.p_ctg.gfa"),
        gfa_haps = expand(os.path.join(OUTMAIN, "assembly", "hifiasm", PRENAME + ".{{filt}}.asm.bp.hap{hap}.p_ctg.gfa"), hap=[1, 2])
    params:
        outprefix = os.path.join(OUTMAIN, "assembly", "hifiasm", PRENAME + ".{filt}.asm"),
        extra_pars = lambda wc: config["filter_settings"][wc.filt]
    threads: 30
    resources:
        mem_mb = 500000,
        runtime = 4000,
        cpus_per_task = 30
    log: os.path.join(OUTMAIN, "assembly", "hifiasm", PRENAME + ".{filt}.asm.log")
    conda:
        "../envs/hifiasm.yml"
    shell: """
    hifiasm -i {params.extra_pars} \
            -o {params.outprefix} \
            -t {threads} {input.fq} 2> {log}
    """
