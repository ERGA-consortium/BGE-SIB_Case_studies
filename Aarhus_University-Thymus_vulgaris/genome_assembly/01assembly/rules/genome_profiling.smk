##############################################
######## TARGET RULES TO RUN PIPELINE ########
##############################################

rule do_genome_profiling:
    input:
        os.path.join(OUTMAIN, "genome_profiling", "genomescope2"),
        os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", PRENAME + "_smudgeplot.png")


###################################
######## RULES TO DO STUFF ########
###################################

rule jellyfish_count:
    input:
        fq = os.path.join(OUTMAIN, "hifi_reads", "filtered_reads", PRENAME + "_hifis_cleaned_seqkit_cutadapt.fq.gz")
    output:
        counts = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", PRENAME + "_hifis_count.jf")
    params:
        k = config["kmer_size"]
    resources:
        mem_mb = 11000,
        runtime = 600,
        cpus_per_task = 10
    threads: 10
    conda:
        "../envs/jellyfish.yml"
    shell: """
        jellyfish count -C -m {params.k} -s 1000000000 -t {threads} <(zcat {input.fq}) -o {output.counts}
    """


rule jellyfish_histogram:
    input:
        counts = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", PRENAME + "_hifis_count.jf")
    output:
        histo = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", PRENAME + "_hifis.histo")
    resources:
        mem_mb = 11000,
        runtime = 300,
        cpus_per_task = 10
    threads: 10
    conda:
        "../envs/jellyfish.yml"
    shell: """
        jellyfish histo -t {threads} {input.counts} > {output.histo}
    """


rule genomescope2:
    input:
        hist = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", PRENAME + "_hifis.histo"),
    output:
        directory(os.path.join(OUTMAIN, "genome_profiling", "genomescope2"))
    params:
        kmer_size = config["kmer_size"],
        prename = PRENAME,
        ploidy = 2
    resources:
        mem_mb = 1000,
        runtime = 100
    conda:
        "../envs/genomescope2.yml"
    shell: """
    genomescope2 -i {input.hist} -o {output} -n {params.prename} -p {params.ploidy} -k {params.kmer_size}
    """


rule get_cutoffs:
    input:
        histo = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", PRENAME + "_hifis.histo"),
    output:
        l_cutoff = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", "l_cutoff.txt"),
        u_cutoff = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", "u_cutoff.txt"),
    resources:
        mem_mb = 1000,
        runtime = 300,
    conda:
        "../envs/smudgeplot.yml"
    shell: """
        smudgeplot cutoff {input.histo} L > {output.l_cutoff}
        smudgeplot cutoff {input.histo} U > {output.u_cutoff}
    """


rule jellyfish_dump:
    input:
        counts = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", PRENAME + "_hifis_count.jf"),
        l_cutoff = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", "l_cutoff.txt"),
        u_cutoff = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", "u_cutoff.txt"),
    output:
        kmers_dump = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", PRENAME + "kmers.dump")
    resources:
        mem_mb = 10000,
        runtime = 300,
    conda:
        "../envs/jellyfish.yml"
    shell: """
        L=$(cat {input.l_cutoff})
        U=$(cat {input.u_cutoff})
        echo $L $U
        jellyfish dump -c -L $L -U $U -o {output.kmers_dump} {input.counts}
    """


rule extract_genomic_kmers:
    input:
        kmers_dump = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", PRENAME + "kmers.dump")
    output:
        kmers_coverages = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", PRENAME + "_kmer_pairs_coverages.tsv")
    resources:
        mem_mb = 250000,
        runtime = 450,
    params:
        outprefix = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", PRENAME + "_kmer_pairs")
    conda:
        "../envs/smudgeplot.yml"
    shell: """
        smudgeplot hetmers -o {params.outprefix} < {input.kmers_dump}
    """


rule plot_smudgeplot:
    input:
        kmers_coverages = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", "kmers", PRENAME + "_kmer_pairs_coverages.tsv")
    output:
        png = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", PRENAME + "_smudgeplot.png")
    resources:
        mem_mb = 5000,
        runtime = 300
    params:
        outprefix = os.path.join(OUTMAIN, "genome_profiling", "smudgeplot", PRENAME)
    conda:
        "../envs/smudgeplot.yml"
    shell: """
        smudgeplot plot {input.kmers_coverages} -o {params.outprefix}
    """
