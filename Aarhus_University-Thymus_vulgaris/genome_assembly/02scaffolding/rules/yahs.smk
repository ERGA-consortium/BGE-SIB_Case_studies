

rule do_scaffolding:
    input:
        fa = expand(os.path.join(OUTMAIN, "scaffolding", OUTPRE + ".hap{hap}.yahs_scaffolds_final.fa"), hap=HAPS),
        fai = expand(os.path.join(OUTMAIN, "scaffolding", OUTPRE + ".hap{hap}.yahs_scaffolds_final.fa.fai"), hap=HAPS),
        hic = expand(os.path.join(OUTMAIN, "scaffolding", "contact_maps", OUTPRE + ".hap{hap}_yahs_juicertools.hic"), hap=HAPS)




rule yahs_scaffolding:
    input:
        fa = REF,
        bam = os.path.join(OUTMAIN, "mapping", "output", OUTPRE + ".hap{hap}_bwa_parse_dedup_split_namesort.bam"),
    output:
        agp = os.path.join(OUTMAIN, "scaffolding", OUTPRE + ".hap{hap}.yahs_scaffolds_final.agp"),
        fa = os.path.join(OUTMAIN, "scaffolding", OUTPRE + ".hap{hap}.yahs_scaffolds_final.fa"),
        binfile = os.path.join(OUTMAIN, "scaffolding", OUTPRE + ".hap{hap}.yahs.bin")
    params:
        outpre = os.path.join(OUTMAIN, "scaffolding", OUTPRE + ".hap{hap}.yahs")
    conda:
        "../envs/yahs.yml"
    threads: 8
    resources:
        mem_mb = 100000,
        runtime = 1000,
        cpus_per_task = 8
    shell:
        """
        yahs -o {params.outpre} {input.fa} {input.bam}
        """
