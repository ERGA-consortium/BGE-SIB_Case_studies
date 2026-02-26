rule agat_keep_longest_isoform:
    input:
        gtf=os.path.join(OUTMAIN, "{assembly}", "braker", "braker.gtf"),
    output:
        gtf=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "{assembly}.longest_isoform.gtf"),
    params:
        tmpdir=config["tmpdir"],
    conda:
        "../envs/agat.yml"
    resources:
        mem_mb=16000,
        runtime=180,
    shadow:
        "shallow"
    shell:
        """
        JOBID="$SLURM_JOB_ID"
        if [ -z "$JOBID" ]; then JOBID="$$"; fi
        JOBTMP="{params.tmpdir}/$JOBID"
        mkdir -p "$JOBTMP"
        export TMPDIR="$JOBTMP"
        agat_sp_keep_longest_isoform.pl \
            --gff {input.gtf} \
            --output {output.gtf}
        """


rule agat_stats_full_annotation:
    input:
        gtf=os.path.join(OUTMAIN, "{assembly}", "braker", "braker.gtf"),
        genome=os.path.join(OUTMAIN, "{assembly}", "genome", "{assembly}.softmasked.fa"),
    output:
        txt=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "{assembly}.agat.full.stats.txt"),
    params:
        tmpdir=config["tmpdir"],
    conda:
        "../envs/agat.yml"
    resources:
        mem_mb=16000,
        runtime=180,
    shadow:
        "shallow"
    shell:
        """
        JOBID="$SLURM_JOB_ID"
        if [ -z "$JOBID" ]; then JOBID="$$"; fi
        JOBTMP="{params.tmpdir}/$JOBID"
        mkdir -p "$JOBTMP"
        export TMPDIR="$JOBTMP"
        agat_sp_statistics.pl \
            --gff {input.gtf} \
            --gs {input.genome} \
            --output {output.txt}
        """


rule agat_stats_longest_annotation:
    input:
        gtf=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "{assembly}.longest_isoform.gtf"),
        genome=os.path.join(OUTMAIN, "{assembly}", "genome", "{assembly}.softmasked.fa"),
    output:
        txt=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "{assembly}.agat.longest.stats.txt"),
    params:
        tmpdir=config["tmpdir"],
    conda:
        "../envs/agat.yml"
    resources:
        mem_mb=16000,
        runtime=180,
    shadow:
        "shallow"
    shell:
        """
        JOBID="$SLURM_JOB_ID"
        if [ -z "$JOBID" ]; then JOBID="$$"; fi
        JOBTMP="{params.tmpdir}/$JOBID"
        mkdir -p "$JOBTMP"
        export TMPDIR="$JOBTMP"
        agat_sp_statistics.pl \
            --gff {input.gtf} \
            --gs {input.genome} \
            --output {output.txt}
        """


rule agat_extract_longest_proteins:
    input:
        gtf=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "{assembly}.longest_isoform.gtf"),
        genome=os.path.join(OUTMAIN, "{assembly}", "genome", "{assembly}.softmasked.fa"),
    output:
        aa=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "{assembly}.longest_isoform.aa"),
    params:
        tmpdir=config["tmpdir"],
    conda:
        "../envs/agat.yml"
    resources:
        mem_mb=16000,
        runtime=180,
    shadow:
        "shallow"
    shell:
        """
        JOBID="$SLURM_JOB_ID"
        if [ -z "$JOBID" ]; then JOBID="$$"; fi
        JOBTMP="{params.tmpdir}/$JOBID"
        mkdir -p "$JOBTMP"
        export TMPDIR="$JOBTMP"
        agat_sp_extract_sequences.pl \
            --gff {input.gtf} \
            --fasta {input.genome} \
            -p \
            --output {output.aa}
        """


rule run_busco_longest_proteins:
    input:
        aa=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "{assembly}.longest_isoform.aa"),
    output:
        outdir=directory(os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "busco", "{assembly}")),
    params:
        lineage=config["annotation_qc"]["busco_lineage"],
        outpath=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "busco"),
    threads: 20
    resources:
        mem_mb=200000,
        runtime=1440,
        cpus_per_task=20,
    conda:
        "../envs/busco.yml"
    shell:
        """
        busco \
            --in {input.aa} \
            --mode proteins \
            --lineage_dataset {params.lineage} \
            --out {wildcards.assembly} \
            --out_path {params.outpath} \
            --cpu {threads}
        """


rule run_omamer_longest_proteins:
    input:
        aa=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "{assembly}.longest_isoform.aa"),
    output:
        omamer=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "omark", "{assembly}.longest_isoform.omamer"),
    params:
        db=config["annotation_qc"]["omamer_db"],
    conda:
        "../envs/omark.yml"
    resources:
        mem_mb=64000,
        runtime=720,
    shell:
        """
        omamer search --db {params.db} --query {input.aa} --out {output.omamer}
        """


rule run_omark_longest_proteins:
    input:
        omamer=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "omark", "{assembly}.longest_isoform.omamer"),
        aa=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "{assembly}.longest_isoform.aa"),
    output:
        outdir=directory(os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "omark", "results")),
    params:
        db=config["annotation_qc"]["omamer_db"],
        taxid=config["annotation_qc"]["omark_taxid"],
        tmpdir=config["tmpdir"],
    conda:
        "../envs/omark.yml"
    resources:
        mem_mb=64000,
        runtime=720,
    shell:
        """
        JOBID="$SLURM_JOB_ID"
        if [ -z "$JOBID" ]; then JOBID="$$"; fi
        JOBHOME="{params.tmpdir}/omark_{wildcards.assembly}_$JOBID"
        mkdir -p "$JOBHOME/.cache"
        export HOME="$JOBHOME"
        export XDG_CACHE_HOME="$JOBHOME/.cache"
        omark -f {input.omamer} -d {params.db} -o {output.outdir} -t {params.taxid} -of {input.aa}
        """


rule summarize_annotation_qc:
    input:
        full_stats=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "{assembly}.agat.full.stats.txt"),
        longest_stats=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "{assembly}.agat.longest.stats.txt"),
        busco_dir=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "busco", "{assembly}"),
        omark_dir=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "omark", "results"),
    output:
        tsv=os.path.join(OUTMAIN, "{assembly}", "annotation_qc", "{assembly}.annotation_qc.summary.tsv"),
    conda:
        "../envs/python.yml"
    shell:
        """
        python scripts/summarize_annotation_qc.py \
            --full-stats {input.full_stats} \
            --longest-stats {input.longest_stats} \
            --busco-dir {input.busco_dir} \
            --omark-dir {input.omark_dir} \
            --out-tsv {output.tsv}
        """
