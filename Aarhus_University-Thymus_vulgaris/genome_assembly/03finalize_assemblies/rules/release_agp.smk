rule do_release_agp_all:
    input:
        expand(
            os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap1_nuclear.agp"),
            sample=config["samples"].keys(),
        ),
        expand(
            os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap1_organelle.agp"),
            sample=config["samples"].keys(),
        ),
        expand(
            os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap1_plus_organelle.agp"),
            sample=config["samples"].keys(),
        ),
        expand(
            os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap2_nuclear.agp"),
            sample=config["samples"].keys(),
        ),


rule build_release_agp_hap1_nuclear:
    input:
        agp=lambda wc: config["samples"][str(wc.sample)]["agp_refs"]["hap1"],
        map_tsv=lambda wc: (
            os.path.join(
                config["outmain"],
                "reference",
                f"{config['reference']['sample']}_hap1.rename_map.tsv",
            )
            if str(wc.sample) == str(config["reference"]["sample"])
            else os.path.join(
                config["outmain"],
                str(wc.sample),
                "hap1",
                "mapping",
                f"{wc.sample}.hap1.to_ref.mapping.tsv",
            )
        ),
        keep_fa=os.path.join(config["outmain"], "{sample}", "hap1", "final", "{sample}.hap1.nuclear.clean.fa"),
    output:
        agp=os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap1_nuclear.agp"),
    params:
        script=os.path.join(workflow.basedir, "scripts", "rewrite_release_agp.py"),
    conda:
        "../envs/python.yml"
    resources:
        mem_mb=4000,
        runtime=120,
    shell: """
        python {params.script} \
            --input-agp {input.agp} \
            --map-tsv {input.map_tsv} \
            --keep-fa {input.keep_fa} \
            --output-agp {output.agp}
    """


rule build_release_agp_hap1_organelle:
    input:
        fai=os.path.join(config["outmain"], "{sample}", "organelle", "{sample}.organelle.fa.fai"),
    output:
        agp=os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap1_organelle.agp"),
    resources:
        mem_mb=1000,
        runtime=30,
    shell: """
        awk 'BEGIN{{OFS="\\t"}} {{print $1,1,$2,1,"W",$1,1,$2,"+"}}' {input.fai} > {output.agp}
    """


rule build_release_agp_hap1_plus_organelle:
    input:
        nuclear=os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap1_nuclear.agp"),
        organelle=os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap1_organelle.agp"),
    output:
        agp=os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap1_plus_organelle.agp"),
    resources:
        mem_mb=1000,
        runtime=20,
    shell: """
        cat {input.nuclear} {input.organelle} > {output.agp}
    """


rule build_release_agp_hap2:
    input:
        agp=lambda wc: config["samples"][str(wc.sample)]["agp_refs"]["hap2"],
        map_tsv=os.path.join(config["outmain"], "{sample}", "hap2", "mapping", "{sample}.hap2.to_ref.mapping.tsv"),
        keep_fa=os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap2_nuclear.fa"),
    output:
        agp=os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap2_nuclear.agp"),
    params:
        script=os.path.join(workflow.basedir, "scripts", "rewrite_release_agp.py"),
    conda:
        "../envs/python.yml"
    resources:
        mem_mb=4000,
        runtime=120,
    shell: """
        python {params.script} \
            --input-agp {input.agp} \
            --map-tsv {input.map_tsv} \
            --keep-fa {input.keep_fa} \
            --output-agp {output.agp}
    """
