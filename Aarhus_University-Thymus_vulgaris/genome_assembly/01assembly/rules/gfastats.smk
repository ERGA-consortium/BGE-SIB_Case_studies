
########################################
######### RULES TO CALL STUFF ##########
########################################

rule do_gfastats_haps:
    input:
        expand(os.path.join(OUTMAIN, "assembly", "hifiasm", PRENAME + ".{filt}.asm.bp.hap1.p_ctg.gfa.stats"), filt=FILTERS),
        expand(os.path.join(OUTMAIN, "assembly", "hifiasm", PRENAME + ".{filt}.asm.bp.hap2.p_ctg.gfa.stats"), filt=FILTERS),
        expand(os.path.join(OUTMAIN, "assembly", "hifiasm", PRENAME + ".{filt}.asm.bp.p_ctg.gfa.stats"), filt=FILTERS)


#######################################
######### RULES TO RUN STUFF ##########
#######################################

rule gfa_to_fasta:
    """convert gfa to fasta format."""
    input:
        gfa = "{prefix}.gfa"
    output:
        fa = "{prefix}.fasta"
    resources:
        mem_mb = 5000,
        runtime = 100,
        cpus_per_task = 3
    threads: 3
    conda:
        "../envs/gfastats.yml"
    shell: """
        gfastats --discover-paths --input-sequence {input.gfa} --threads {threads} -o {output.fa}
    """


rule gfa_stats:
    input:
        gfa = "{prefix}.gfa"
    output:
        stats = "{prefix}.gfa.stats"
    resources:
        mem_mb = 5000,
        runtime = 100,
        cpus_per_task = 3
    threads: 3
    conda:
        "../envs/gfastats.yml"
    shell: """
    gfastats --discover-paths {input.gfa} > {output.stats}
    """


rule fasta_stats:
    input:
        fas = "{prefix}.fa"
    output:
        stats = "{prefix}.fa.stats"
    resources:
        mem_mb = 5000,
        runtime = 100,
        cpus_per_task = 3
    threads: 3
    conda:
        "../envs/gfastats.yml"
    shell: """
    gfastats {input.fas} > {output.stats}
    """
