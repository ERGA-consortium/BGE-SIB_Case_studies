TMP_BASE = config.get("tmpdir", "/tmp")


rule do_hic_contacts:
    input:
        expand(
            os.path.join(OUTMAIN, "scaffolding", "contact_maps", OUTPRE + ".hap{hap}_yahs_juicertools.hic"),
            hap=HAPS,
        )


rule do_jbat_prep:
    input:
        expand(
            os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.hic"),
            hap=HAPS,
        ),
        expand(
            os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.assembly"),
            hap=HAPS,
        ),
        expand(
            os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.liftover.agp"),
            hap=HAPS,
        )


rule do_post_jbat_curation:
    input:
        expand(
            os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.FINAL.agp"),
            hap=HAPS_POST,
        ),
        expand(
            os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.FINAL.fa"),
            hap=HAPS_POST,
        ),
        expand(
            os.path.join(OUTMAIN, "scaffolding", "contact_maps", OUTPRE + ".hap{hap}_curated_yahs_juicertools.hic"),
            hap=HAPS_POST,
        )




rule juicer_pre_contacts:
    input:
        agp = os.path.join(OUTMAIN, "scaffolding", OUTPRE + ".hap{hap}.yahs_scaffolds_final.agp"),
        binfile = os.path.join(OUTMAIN, "scaffolding", OUTPRE + ".hap{hap}.yahs.bin"),
        fai = lambda wc: REF(wc) + ".fai",
    output:
        aln = os.path.join(OUTMAIN, "scaffolding", "contact_maps", OUTPRE + ".hap{hap}_alignments_sorted.txt")
    params:
        tmp_prefix = lambda wc: os.path.join(TMP_BASE, f"{OUTPRE}.hap{wc.hap}.juicer_pre"),
    conda:
        "../envs/yahs.yml"
    threads: 8
    resources:
        mem_mb = 50000,
        runtime = 1000,
        cpus_per_task = 8
    shell: """
        JOBTAG=${{SLURM_JOB_ID:-local}}
        TMPDIR="{params.tmp_prefix}.${{JOBTAG}}"
        mkdir -p "${{TMPDIR}}"

        (juicer pre {input.binfile} {input.agp} {input.fai} | \
            sort -k2,2d -k6,6d -T "${{TMPDIR}}" --parallel={threads} -S32G | \
            awk 'NF' > \
            {output.aln}.part) && \
            (mv {output.aln}.part {output.aln})
    """


rule get_scaffold_sizes_contacts:
    input:
        fa = os.path.join(OUTMAIN, "scaffolding", OUTPRE + ".hap{hap}.yahs_scaffolds_final.fa.fai"),
    output:
        chrom_sizes = os.path.join(OUTMAIN, "scaffolding", "contact_maps", OUTPRE + ".hap{hap}.yahs_scaffolds_final.chrom_sizes.txt")
    resources:
        mem_mb = 5000,
        runtime = 100,
    shell: """
        cut -f1,2 {input.fa} > {output.chrom_sizes}
    """


rule juicertools_pre_contacts:
    input:
        aln = os.path.join(OUTMAIN, "scaffolding", "contact_maps", OUTPRE + ".hap{hap}_alignments_sorted.txt"),
        chrom_sizes = os.path.join(OUTMAIN, "scaffolding", "contact_maps", OUTPRE + ".hap{hap}.yahs_scaffolds_final.chrom_sizes.txt")
    output:
        hic = os.path.join(OUTMAIN, "scaffolding", "contact_maps", OUTPRE + ".hap{hap}_yahs_juicertools.hic")
    params:
        tmp_prefix = lambda wc: os.path.join(TMP_BASE, f"{OUTPRE}.hap{wc.hap}.juicertools_pre_contacts"),
    conda:
        "../envs/juicertools.yml"
    threads: 8
    resources:
        mem_mb = 50000,
        runtime = 1000,
        cpus_per_task = 8
    shell: """
        mkdir -p $(dirname {output.hic})
        JOBTAG=${{SLURM_JOB_ID:-local}}
        TMPDIR="{params.tmp_prefix}.${{JOBTAG}}"
        TMPHIC="${{TMPDIR}}/out.hic.part"
        mkdir -p "${{TMPDIR}}"

        export _JAVA_OPTIONS="-Xms2g -Xmx{resources.mem_mb}m -Djava.io.tmpdir=${{TMPDIR}}"
        (juicer_tools pre -j {threads} {input.aln} "${{TMPHIC}}" {input.chrom_sizes}) && \
            (mv "${{TMPHIC}}" {output.hic})
    """


rule juicer_pre_jbat:
    input:
        agp = os.path.join(OUTMAIN, "scaffolding", OUTPRE + ".hap{hap}.yahs_scaffolds_final.agp"),
        binfile = os.path.join(OUTMAIN, "scaffolding", OUTPRE + ".hap{hap}.yahs.bin"),
        fai = lambda wc: REF(wc) + ".fai",
    output:
        txt = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.txt"),
        assembly = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.assembly"),
        liftover = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.liftover.agp"),
        log = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.log"),
    params:
        outprefix = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT")
    conda:
        "../envs/yahs.yml"
    threads: 8
    resources:
        mem_mb = 50000,
        runtime = 1000,
        cpus_per_task = 8
    shell: """
        mkdir -p $(dirname {params.outprefix})
        juicer pre -a -o {params.outprefix} {input.binfile} {input.agp} {input.fai} > {output.log} 2>&1
    """


rule jbat_pre_csize:
    input:
        log = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.log")
    output:
        chrom_sizes = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.chrom.sizes")
    resources:
        mem_mb = 1000,
        runtime = 100,
    shell: """
        grep PRE_C_SIZE {input.log} | awk '{{print $2" "$3}}' > {output.chrom_sizes}
    """


rule juicertools_pre_jbat:
    input:
        txt = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.txt"),
        chrom_sizes = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.chrom.sizes")
    output:
        hic = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.hic")
    params:
        tmp_prefix = lambda wc: os.path.join(TMP_BASE, f"{OUTPRE}.hap{wc.hap}.juicertools_pre_jbat"),
    conda:
        "../envs/juicertools.yml"
    threads: 2
    resources:
        mem_mb = 50000,
        runtime = 1000,
        cpus_per_task = 2
    shell: """
        mkdir -p $(dirname {output.hic})
        JOBTAG=${{SLURM_JOB_ID:-local}}
        TMPDIR="{params.tmp_prefix}.${{JOBTAG}}"
        TMPHIC="${{TMPDIR}}/out.hic.part"
        mkdir -p "${{TMPDIR}}"

        export _JAVA_OPTIONS="-Xms2g -Xmx{resources.mem_mb}m -Djava.io.tmpdir=${{TMPDIR}}"
        (juicer_tools pre -n -j {threads} {input.txt} "${{TMPHIC}}" {input.chrom_sizes}) && \
            (mv "${{TMPHIC}}" {output.hic})
    """


rule juicer_post_jbat:
    input:
        review_assembly = lambda wc: REVIEW_ASSEMBLY[str(wc.hap)],
        liftover = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.liftover.agp"),
        ref = REF,
    output:
        final_agp = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.FINAL.agp"),
        final_fa = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.FINAL.fa")
    params:
        outprefix = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT")
    conda:
        "../envs/yahs.yml"
    threads: 2
    resources:
        mem_mb = 10000,
        runtime = 500,
        cpus_per_task = 2
    shell: """
        mkdir -p $(dirname {params.outprefix})
        juicer post -o {params.outprefix} {input.review_assembly} {input.liftover} {input.ref}
    """


rule juicer_pre_contacts_curated:
    input:
        agp = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.FINAL.agp"),
        binfile = os.path.join(OUTMAIN, "scaffolding", OUTPRE + ".hap{hap}.yahs.bin"),
        fai = lambda wc: REF(wc) + ".fai",
    output:
        aln = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "curated", OUTPRE + ".hap{hap}_curated_alignments_sorted.txt")
    params:
        tmp_prefix = lambda wc: os.path.join(TMP_BASE, f"{OUTPRE}.hap{wc.hap}.juicer_pre_curated"),
    conda:
        "../envs/yahs.yml"
    threads: 8
    resources:
        mem_mb = 50000,
        runtime = 1000,
        cpus_per_task = 8
    shell: """
        mkdir -p $(dirname {output.aln})
        JOBTAG=${{SLURM_JOB_ID:-local}}
        TMPDIR="{params.tmp_prefix}.${{JOBTAG}}"
        mkdir -p "${{TMPDIR}}"

        (juicer pre {input.binfile} {input.agp} {input.fai} | \
            sort -k2,2d -k6,6d -T "${{TMPDIR}}" --parallel={threads} -S32G | \
            awk 'NF' > \
            {output.aln}.part) && \
            (mv {output.aln}.part {output.aln})
    """


rule get_scaffold_sizes_contacts_curated:
    input:
        agp = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "jbat", OUTPRE + ".hap{hap}.out_JBAT.FINAL.agp"),
    output:
        chrom_sizes = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "curated", OUTPRE + ".hap{hap}.curated.chrom_sizes.txt")
    resources:
        mem_mb = 5000,
        runtime = 100,
    shell: """
        mkdir -p $(dirname {output.chrom_sizes})
        awk 'BEGIN{{OFS="\\t"}} !/^#/ {{if($3>len[$1]) len[$1]=$3}} END {{for(k in len) print k,len[k]}}' {input.agp} | \
            sort -k2,2nr -k1,1 > {output.chrom_sizes}
    """


rule juicertools_pre_contacts_curated:
    input:
        aln = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "curated", OUTPRE + ".hap{hap}_curated_alignments_sorted.txt"),
        chrom_sizes = os.path.join(OUTMAIN, "scaffolding", "contact_maps", "curated", OUTPRE + ".hap{hap}.curated.chrom_sizes.txt")
    output:
        hic = os.path.join(OUTMAIN, "scaffolding", "contact_maps", OUTPRE + ".hap{hap}_curated_yahs_juicertools.hic")
    params:
        tmp_prefix = lambda wc: os.path.join(TMP_BASE, f"{OUTPRE}.hap{wc.hap}.juicertools_pre_contacts_curated"),
    conda:
        "../envs/juicertools.yml"
    threads: 8
    resources:
        mem_mb = 50000,
        runtime = 1000,
        cpus_per_task = 8
    shell: """
        mkdir -p $(dirname {output.hic})
        JOBTAG=${{SLURM_JOB_ID:-local}}
        TMPDIR="{params.tmp_prefix}.${{JOBTAG}}"
        TMPHIC="${{TMPDIR}}/out.hic.part"
        mkdir -p "${{TMPDIR}}"

        export _JAVA_OPTIONS="-Xms2g -Xmx{resources.mem_mb}m -Djava.io.tmpdir=${{TMPDIR}}"
        (juicer_tools pre -j {threads} {input.aln} "${{TMPHIC}}" {input.chrom_sizes}) && \
            (mv "${{TMPHIC}}" {output.hic})
    """
