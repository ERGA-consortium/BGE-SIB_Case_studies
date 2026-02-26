


rule concat_protein_fastas:
    input:
        config["protein_evidence"]["input_fastas"],
    output:
        fa=os.path.join(OUTMAIN, "protein_evidence", "proteins.concat.fa"),
    resources:
        mem_mb=8000,
        runtime=120,
    shell:
        """
        cat {input} > {output.fa}
        """


rule cluster_protein_fastas:
    input:
        fa=os.path.join(OUTMAIN, "protein_evidence", "proteins.concat.fa"),
    output:
        cluster_tsv=os.path.join(OUTMAIN, "protein_evidence", "proteins_cluster.tsv"),
        fa=os.path.join(OUTMAIN, "protein_evidence", "proteins.clustered.fa"),
    params:
        outprefix=os.path.join(OUTMAIN, "protein_evidence", "proteins"),
        tmpdir=config["tmpdir"],
        min_seq_id=config["protein_evidence"]["cluster"]["min_seq_id"],
        coverage=config["protein_evidence"]["cluster"]["coverage"],
        cov_mode=config["protein_evidence"]["cluster"]["cov_mode"],
    threads: config["protein_evidence"]["cluster"]["threads"]
    resources:
        mem_mb=64000,
        runtime=600,
        cpus_per_task=config["protein_evidence"]["cluster"]["threads"],
    conda:
        "../envs/mmseqs.yml"
    shell:
        """
        mmseqs easy-linclust {input.fa} {params.outprefix} {params.tmpdir}/${{SLURM_JOB_ID}} \
            --min-seq-id {params.min_seq_id} -c {params.coverage} \
            --cov-mode {params.cov_mode} \
            --threads {threads}
        cp {params.outprefix}_rep_seq.fasta {output.fa}
        """


rule summarize_protein_clustering:
    input:
        input_fa=os.path.join(OUTMAIN, "protein_evidence", "proteins.concat.fa"),
        clustered_fa=os.path.join(OUTMAIN, "protein_evidence", "proteins.clustered.fa"),
        cluster_tsv=os.path.join(OUTMAIN, "protein_evidence", "proteins_cluster.tsv"),
    output:
        tsv=os.path.join(OUTMAIN, "protein_evidence", "proteins.cluster_summary.tsv"),
    params:
        script=os.path.join(workflow.basedir, "scripts", "summarize_protein_clustering.py"),
    resources:
        mem_mb=4000,
        runtime=60,
        cpus_per_task=1,
    conda:
        "../envs/python.yml"
    shell:
        """
        python {params.script} \
            --input-fa {input.input_fa} \
            --clustered-fa {input.clustered_fa} \
            --cluster-tsv {input.cluster_tsv} \
            --out-tsv {output.tsv}
        """
