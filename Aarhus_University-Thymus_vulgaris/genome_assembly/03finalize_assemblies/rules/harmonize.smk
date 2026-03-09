rule do_harmonize_all:
    input:
        os.path.join(
            config["outmain"],
            "reference",
            f"{config['reference']['sample']}_hap1.lg_reference.fa.fai",
        ),
        expand(
            os.path.join(
                config["outmain"],
                "{sample}",
                "{hap}",
                "harmonized",
                "{sample}.{hap}.harmonized.fa.fai",
            ),
            zip,
            sample=[
                str(sample)
                for sample, sample_cfg in config["samples"].items()
                for hap in sample_cfg["haps"].keys()
                if not (
                    str(sample) == str(config["reference"]["sample"])
                    and str(hap) == "hap1"
                )
            ],
            hap=[
                str(hap)
                for sample, sample_cfg in config["samples"].items()
                for hap in sample_cfg["haps"].keys()
                if not (
                    str(sample) == str(config["reference"]["sample"])
                    and str(hap) == "hap1"
                )
            ],
        ),


rule build_reference_lg:
    input:
        fa=lambda wc: os.path.join(
            config["outmain"],
            "input",
            str(config["reference"]["sample"]),
            "hap1.fa",
        ),
    output:
        lg_ref_fa=os.path.join(
            config["outmain"],
            "reference",
            f"{config['reference']['sample']}_hap1.lg_reference.fa",
        ),
        map_tsv=os.path.join(
            config["outmain"],
            "reference",
            f"{config['reference']['sample']}_hap1.rename_map.tsv",
        ),
        unplaced_tsv=os.path.join(
            config["outmain"],
            "reference",
            f"{config['reference']['sample']}_hap1.unplaced.tsv",
        ),
    params:
        script=os.path.join(workflow.basedir, "scripts", "build_lg_reference.py"),
        min_bp=lambda wc: int(config.get("harmonize", {}).get("lg_min_bp", 10000000)),
        lg_prefix=lambda wc: config.get("harmonize", {}).get("lg_prefix", "chr"),
        unplaced_prefix=lambda wc: config.get("harmonize", {}).get("unplaced_prefix", "unplaced"),
    conda:
        "../envs/python.yml"
    resources:
        mem_mb=8000,
        runtime=200,
    shell: """
        python {params.script} \
            --input-fa {input.fa} \
            --output-lg-fa {output.lg_ref_fa} \
            --output-map {output.map_tsv} \
            --output-unplaced {output.unplaced_tsv} \
            --lg-min-bp {params.min_bp} \
            --lg-prefix {params.lg_prefix} \
            --unplaced-prefix {params.unplaced_prefix}
    """


rule map_to_reference_lg:
    input:
        query=os.path.join(config["outmain"], "input", "{sample}", "{hap}.fa"),
        ref=os.path.join(
            config["outmain"],
            "reference",
            f"{config['reference']['sample']}_hap1.lg_reference.fa",
        ),
    output:
        paf=os.path.join(config["outmain"], "{sample}", "{hap}", "alignments", "{sample}.{hap}.to_ref.paf"),
    params:
        preset=lambda wc: config.get("harmonize", {}).get("synteny", {}).get("minimap2_preset", "asm5"),
    conda:
        "../envs/minimap2.yml"
    threads: 10
    resources:
        mem_mb=220000,
        runtime=1400,
        cpus_per_task=10,
    shell: """
        minimap2 -x {params.preset} -t {threads} {input.ref} {input.query} > {output.paf}
    """


rule assign_lg_from_paf:
    input:
        paf=os.path.join(config["outmain"], "{sample}", "{hap}", "alignments", "{sample}.{hap}.to_ref.paf"),
        query=os.path.join(config["outmain"], "input", "{sample}", "{hap}.fa"),
    output:
        map_tsv=os.path.join(config["outmain"], "{sample}", "{hap}", "mapping", "{sample}.{hap}.to_ref.mapping.tsv"),
        summary_tsv=os.path.join(config["outmain"], "{sample}", "{hap}", "mapping", "{sample}.{hap}.to_ref.summary.tsv"),
    params:
        script=os.path.join(workflow.basedir, "scripts", "assign_lg_from_paf.py"),
        min_query_bp=lambda wc: int(config.get("harmonize", {}).get("lg_min_bp", 10000000)),
        unplaced_prefix=lambda wc: config.get("harmonize", {}).get("unplaced_prefix", "unplaced"),
    conda:
        "../envs/python.yml"
    resources:
        mem_mb=8000,
        runtime=300,
    shell: """
        python {params.script} \
            --paf {input.paf} \
            --query-fa {input.query} \
            --out-map {output.map_tsv} \
            --out-summary {output.summary_tsv} \
            --min-query-bp {params.min_query_bp} \
            --unplaced-prefix {params.unplaced_prefix}
    """


rule rewrite_harmonized_fasta:
    input:
        fa=os.path.join(config["outmain"], "input", "{sample}", "{hap}.fa"),
        map_tsv=lambda wc: (
            os.path.join(
                config["outmain"],
                "reference",
                f"{config['reference']['sample']}_hap1.rename_map.tsv",
            )
            if (
                str(wc.sample) == str(config["reference"]["sample"])
                and str(wc.hap) == "hap1"
            )
            else os.path.join(
                config["outmain"],
                str(wc.sample),
                str(wc.hap),
                "mapping",
                f"{wc.sample}.{wc.hap}.to_ref.mapping.tsv",
            )
        ),
    output:
        fa=os.path.join(config["outmain"], "{sample}", "{hap}", "harmonized", "{sample}.{hap}.harmonized.fa"),
    params:
        script=os.path.join(workflow.basedir, "scripts", "rewrite_fasta_by_mapping.py"),
    conda:
        "../envs/python.yml"
    resources:
        mem_mb=4000,
        runtime=120,
    shell: """
        python {params.script} \
            --input-fa {input.fa} \
            --map-tsv {input.map_tsv} \
            --output-fa {output.fa}
    """
