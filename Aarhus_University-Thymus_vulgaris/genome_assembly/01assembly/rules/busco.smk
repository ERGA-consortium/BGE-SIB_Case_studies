
########################################
######### RULES TO CALL STUFF ##########
########################################

rule do_busco_haps:
    input:
        expand(os.path.join(OUTMAIN, "assembly", "evaluations", "busco", "plots", "{filt}", "hap{hap}", "busco_figure.png"),
               filt=FILTERS, hap=[1, 2])


#######################################
######### RULES TO RUN STUFF ##########
#######################################

rule busco_hap:
    wildcard_constraints:
        hap = "1|2"
    input:
        fasta = os.path.join(OUTMAIN, "assembly", "hifiasm", PRENAME + ".{filt}.asm.bp.hap{hap}.p_ctg.fasta")
    output:
        outdir = directory(os.path.join(OUTMAIN, "assembly", "evaluations", "busco", "{filt}", "hap{hap}"))
    params:
        busco_db = config["busco_db"]
    resources:
        mem_mb = 50000,
        runtime = 1200,
        cpus_per_task = 10
    threads: 10
    conda:
        "../envs/busco.yml"
    shell: """
        busco --in {input.fasta} --out {output.outdir} --mode genome --lineage_dataset {params.busco_db} --cpu {threads}
    """


rule busco_plot:
    input:
        indir = os.path.join(OUTMAIN, "assembly", "evaluations", "busco", "{dirpath}")
    output:
        os.path.join(OUTMAIN, "assembly", "evaluations", "busco", "plots", "{dirpath}", "busco_figure.png")
    params:
        workdir = os.path.join(OUTMAIN, "assembly", "evaluations", "busco", "plots", "{dirpath}")
    resources:
        mem_mb = 1000,
        runtime = 100
    threads: 10
    conda:
        "../envs/busco.yml"
    shell: """
    cp {input.indir}/short_summary.*.json {params.workdir}/
    busco --plot {params.workdir}
    """
