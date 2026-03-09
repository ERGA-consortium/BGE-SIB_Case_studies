rule do_contamination_all:
    input:
        expand(
            os.path.join(
                config["outmain"],
                "{sample}",
                "{hap}",
                "final",
                "{sample}.{hap}.nuclear.clean.fa.fai",
            ),
            sample=config["samples"].keys(),
            hap=["hap1", "hap2"],
        ),


rule run_fcs_gx:
    input:
        fa=os.path.join(config["outmain"], "{sample}", "{hap}", "final", "{sample}.{hap}.nuclear.no_organelle.fa"),
    output:
        report=os.path.join(
            config["outmain"],
            "{sample}",
            "{hap}",
            "contamination",
            "fcs_gx",
            "{sample}.{hap}.fcs_gx_report.txt",
        ),
        taxonomy=os.path.join(
            config["outmain"],
            "{sample}",
            "{hap}",
            "contamination",
            "fcs_gx",
            "{sample}.{hap}.taxonomy.rpt",
        ),
    params:
        gx_db=config["contamination"]["fcs_gx_db_prefix"],
        tax_id=config["contamination"]["tax_id"],
        outdir=os.path.join(config["outmain"], "{sample}", "{hap}", "contamination", "fcs_gx"),
        split_fasta="true",
    conda:
        "../envs/fcs-gx.yml"
    threads: 20
    resources:
        mem_mb=220000,
        runtime=1440,
        cpus_per_task=20,
    shell: """
        GX_NUM_CORES={threads} run_gx.py \
            --fasta {input.fa} \
            --tax-id {params.tax_id} \
            --gx-db {params.gx_db} \
            --split-fasta {params.split_fasta} \
            --action-report true \
            --out-dir {params.outdir} \
            --out-basename {wildcards.sample}.{wildcards.hap}
    """


rule select_bacterial_contigs_from_fcs:
    input:
        report=os.path.join(config["outmain"], "{sample}", "{hap}", "contamination", "fcs_gx", "{sample}.{hap}.fcs_gx_report.txt"),
        fa=os.path.join(config["outmain"], "{sample}", "{hap}", "final", "{sample}.{hap}.nuclear.no_organelle.fa"),
    output:
        remove_list=os.path.join(config["outmain"], "{sample}", "{hap}", "contamination", "{sample}.{hap}.contamination.remove.txt"),
        summary=os.path.join(config["outmain"], "{sample}", "{hap}", "contamination", "{sample}.{hap}.contamination.summary.tsv"),
    params:
        script=os.path.join(workflow.basedir, "scripts", "select_fcs_bacterial_contigs.py"),
    conda:
        "../envs/python.yml"
    resources:
        mem_mb=4000,
        runtime=120,
    shell: """
        python {params.script} \
            --report {input.report} \
            --input-fa {input.fa} \
            --out-remove {output.remove_list} \
            --out-summary {output.summary} \
            --remove-action EXCLUDE
    """


rule write_contamination_clean_fasta:
    input:
        fa=os.path.join(config["outmain"], "{sample}", "{hap}", "final", "{sample}.{hap}.nuclear.no_organelle.fa"),
        remove_list=os.path.join(config["outmain"], "{sample}", "{hap}", "contamination", "{sample}.{hap}.contamination.remove.txt"),
    output:
        fa=os.path.join(config["outmain"], "{sample}", "{hap}", "final", "{sample}.{hap}.nuclear.clean.fa"),
    params:
        script=os.path.join(workflow.basedir, "scripts", "filter_fasta_by_remove_list.py"),
    conda:
        "../envs/python.yml"
    resources:
        mem_mb=4000,
        runtime=120,
    shell: """
        python {params.script} \
            --input-fa {input.fa} \
            --remove-list {input.remove_list} \
            --output-fa {output.fa} \
            --mode remove_contigs
    """
