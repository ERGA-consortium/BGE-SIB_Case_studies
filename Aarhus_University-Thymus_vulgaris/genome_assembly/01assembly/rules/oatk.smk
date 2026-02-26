rule do_oatk:
    input:
        mito = os.path.join(OUTMAIN, "organelles", PRENAME + ".mito.ctg.fasta"),
        plstd = os.path.join(OUTMAIN, "organelles", PRENAME + ".pltd.ctg.fasta")


rule run_oatk:
    input:
        fq = os.path.join(OUTMAIN, "hifi_reads", "filtered_reads", PRENAME + "_hifis_cleaned_seqkit_cutadapt.fq.gz"),
        mitodb  = config["oatk"]["mitodb"],
        plstdb = config["oatk"]["plstdb"]
    output:
        gfa = os.path.join(OUTMAIN, "organelles", PRENAME + ".utg.final.gfa"),
        mito = os.path.join(OUTMAIN, "organelles", PRENAME + ".mito.ctg.fasta"),
        plstdb = os.path.join(OUTMAIN, "organelles", PRENAME + ".pltd.ctg.fasta")
    params:
        cov = 150, # coverage threshold, should be 5-10 times nuclear coverage
        syncmer = 1001, # syncmer 
        outprefix = os.path.join(OUTMAIN, "organelles", PRENAME)
    threads: 30
    resources: 
        mem_mb = 200000,
        runtime = 600,
        cpus_per_task = 30
    conda:
        "../envs/oatk.yml"
    shell: """
     oatk -k {params.syncmer} -c {params.cov} -t {threads} \
        -m {input.mitodb} -p {input.plstdb} \
        -o {params.outprefix} {input.fq}
    """
