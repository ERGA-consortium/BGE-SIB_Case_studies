
##############################################
######## TARGET RULES TO RUN PIPELINE ########
##############################################

rule do_merqury_haps:
    input:
        expand(os.path.join(OUTMAIN, "assembly", "evaluations", "merqury", "{filt}.haps." + PRENAME + ".{filt}.asm.bp.hap{hap}.p_ctg.spectra-cn.fl.png"),
               filt=FILTERS, hap=[1, 2]),
        expand(os.path.join(OUTMAIN, "assembly", "evaluations", "merqury", "{filt}.haps.qv"), filt=FILTERS),
        expand(os.path.join(OUTMAIN, "assembly", "evaluations", "merqury", "{filt}.haps.completeness.stats"), filt=FILTERS)


###################################
######## RULES TO DO STUFF ########
###################################

rule meryl_count_hifis:
    input:
        fq = os.path.join(OUTMAIN, "hifi_reads", "filtered_reads", PRENAME + "_hifis_cleaned_seqkit_cutadapt.fq.gz")
    output:
        meryldb = directory(os.path.join(OUTMAIN, "genome_profiling", "kmers", PRENAME + "_hifis.meryl"))
    params:
        k = config["kmer_size"]
    resources:
        mem_mb = 11000,
        runtime = 600,
        cpus_per_task = 10
    threads: 10
    conda:
        "../envs/meryl.yml"
    shell: """
        meryl k={params.k} threads={threads} memory=10 count {input.fq} output {output.meryldb}
    """


rule merqury_haps:
    input:
        meryl_db = os.path.join(OUTMAIN, "genome_profiling", "kmers", PRENAME + "_hifis.meryl"),
        hap1 = os.path.join(OUTMAIN, "assembly", "hifiasm", PRENAME + ".{filt}.asm.bp.hap1.p_ctg.fasta"),
        hap2 = os.path.join(OUTMAIN, "assembly", "hifiasm", PRENAME + ".{filt}.asm.bp.hap2.p_ctg.fasta")
    output:
        kmer_mult_plot_hap1 = os.path.join(OUTMAIN, "assembly", "evaluations", "merqury", "{filt}.haps." + PRENAME + ".{filt}.asm.bp.hap1.p_ctg.spectra-cn.fl.png"),
        kmer_mult_plot_hap2 = os.path.join(OUTMAIN, "assembly", "evaluations", "merqury", "{filt}.haps." + PRENAME + ".{filt}.asm.bp.hap2.p_ctg.spectra-cn.fl.png"),
        qv = os.path.join(OUTMAIN, "assembly", "evaluations", "merqury", "{filt}.haps.qv"),
        completeness = os.path.join(OUTMAIN, "assembly", "evaluations", "merqury", "{filt}.haps.completeness.stats")
    params:
        outdir = os.path.join(OUTMAIN, "assembly", "evaluations", "merqury"),
        outprefix = "{filt}.haps"
    conda:
        "../envs/merqury.yml"
    resources:
        mem_mb = 11000,
        runtime = 400,
        cpus_per_task = 10
    shell: """
        meryl_db=`realpath {input.meryl_db}`
        hap1=`realpath {input.hap1}`
        hap2=`realpath {input.hap2}`
        cd {params.outdir}
        merqury.sh $meryl_db $hap1 $hap2 {params.outprefix}
    """
