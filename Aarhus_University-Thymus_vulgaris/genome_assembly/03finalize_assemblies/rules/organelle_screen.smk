rule do_organelle_screen_all:
    input:
        expand(
            os.path.join(config["outmain"], "{sample}", "{hap}", "final", "{sample}.{hap}.nuclear.no_organelle.fa.fai"),
            sample=config["samples"].keys(),
            hap=["hap1", "hap2"],
        ),


rule do_finalize_all:
    input:
        expand(
            os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap1_plus_organelle.fa.fai"),
            sample=config["samples"].keys(),
        ),
        expand(
            os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap2_nuclear.fa.fai"),
            sample=config["samples"].keys(),
        ),
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


rule build_organelle_reference:
    input:
        plastid=lambda wc: config["samples"][str(wc.sample)]["organelle_refs"]["plastid"],
        mitochondria=lambda wc: config["samples"][str(wc.sample)]["organelle_refs"]["mitochondria"],
    output:
        fa=os.path.join(config["outmain"], "{sample}", "organelle", "{sample}.organelle.fa"),
    params:
        script=os.path.join(workflow.basedir, "scripts", "build_organelle_reference.py"),
        plastid_complete=lambda wc: config.get("organelle_screen", {}).get("naming", {}).get("plastid_complete_name", "chrC"),
        mito_complete=lambda wc: config.get("organelle_screen", {}).get("naming", {}).get("mitochondria_complete_name", "chrM"),
        plastid_fragment_prefix=lambda wc: config.get("organelle_screen", {}).get("naming", {}).get("plastid_fragment_prefix", "pt"),
        mito_fragment_prefix=lambda wc: config.get("organelle_screen", {}).get("naming", {}).get("mitochondria_fragment_prefix", "mt"),
        plastid_min_complete_bp=lambda wc: int(config.get("organelle_screen", {}).get("completeness", {}).get("plastid_min_complete_bp", 120000)),
        mito_min_complete_bp=lambda wc: int(config.get("organelle_screen", {}).get("completeness", {}).get("mitochondria_min_complete_bp", 100000)),
    conda:
        "../envs/python.yml"
    resources:
        mem_mb=2000,
        runtime=60,
    shell: """
        python {params.script} \
            --plastid-fa {input.plastid} \
            --mitochondria-fa {input.mitochondria} \
            --output-fa {output.fa} \
            --plastid-complete-name {params.plastid_complete} \
            --mitochondria-complete-name {params.mito_complete} \
            --plastid-fragment-prefix {params.plastid_fragment_prefix} \
            --mitochondria-fragment-prefix {params.mito_fragment_prefix} \
            --plastid-min-complete-bp {params.plastid_min_complete_bp} \
            --mitochondria-min-complete-bp {params.mito_min_complete_bp}
    """


rule align_organelle_candidates:
    input:
        query=os.path.join(config["outmain"], "{sample}", "{hap}", "harmonized", "{sample}.{hap}.harmonized.fa"),
        ref=os.path.join(config["outmain"], "{sample}", "organelle", "{sample}.organelle.fa"),
    output:
        paf=os.path.join(config["outmain"], "{sample}", "{hap}", "organelle_screen", "alignments", "{sample}.{hap}.vs_organelle.paf"),
    params:
        preset=lambda wc: config.get("organelle_screen", {}).get("minimap2_preset", "asm5"),
    conda:
        "../envs/minimap2.yml"
    threads: 10
    resources:
        mem_mb=220000,
        runtime=1000,
        cpus_per_task=10,
    shell: """
        minimap2 -x {params.preset} -t {threads} {input.ref} {input.query} > {output.paf}
    """


rule summarize_organelle_hits:
    input:
        paf=os.path.join(config["outmain"], "{sample}", "{hap}", "organelle_screen", "alignments", "{sample}.{hap}.vs_organelle.paf"),
        query=os.path.join(config["outmain"], "{sample}", "{hap}", "harmonized", "{sample}.{hap}.harmonized.fa"),
    output:
        remove_list=os.path.join(config["outmain"], "{sample}", "{hap}", "organelle_screen", "{sample}.{hap}.contigs_to_remove.txt"),
        bed=os.path.join(config["outmain"], "{sample}", "{hap}", "organelle_screen", "{sample}.{hap}.organelle_like_segments.bed"),
        summary=os.path.join(config["outmain"], "{sample}", "{hap}", "organelle_screen", "{sample}.{hap}.summary.tsv"),
    params:
        script=os.path.join(workflow.basedir, "scripts", "screen_organelle_hits.py"),
        min_identity=lambda wc: float(config.get("organelle_screen", {}).get("min_identity", 0.98)),
        min_aln_bp=lambda wc: int(config.get("organelle_screen", {}).get("min_aln_bp", 5000)),
        min_mapq=lambda wc: int(config.get("organelle_screen", {}).get("min_mapq", 20)),
        min_cov_frac=lambda wc: float(config.get("organelle_screen", {}).get("min_cov_frac", 0.90)),
        min_cov_bp=lambda wc: int(config.get("organelle_screen", {}).get("min_cov_bp", 50000)),
        plastid_complete=lambda wc: config.get("organelle_screen", {}).get("naming", {}).get("plastid_complete_name", "chrC"),
        mito_complete=lambda wc: config.get("organelle_screen", {}).get("naming", {}).get("mitochondria_complete_name", "chrM"),
        plastid_fragment_prefix=lambda wc: config.get("organelle_screen", {}).get("naming", {}).get("plastid_fragment_prefix", "pt"),
        mito_fragment_prefix=lambda wc: config.get("organelle_screen", {}).get("naming", {}).get("mitochondria_fragment_prefix", "mt"),
    conda:
        "../envs/python.yml"
    resources:
        mem_mb=4000,
        runtime=120,
    shell: """
        python {params.script} \
            --paf {input.paf} \
            --query-fa {input.query} \
            --out-remove {output.remove_list} \
            --out-bed {output.bed} \
            --out-summary {output.summary} \
            --min-identity {params.min_identity} \
            --min-aln-bp {params.min_aln_bp} \
            --min-mapq {params.min_mapq} \
            --min-cov-frac {params.min_cov_frac} \
            --min-cov-bp {params.min_cov_bp} \
            --plastid-complete-name {params.plastid_complete} \
            --mitochondria-complete-name {params.mito_complete} \
            --plastid-fragment-prefix {params.plastid_fragment_prefix} \
            --mitochondria-fragment-prefix {params.mito_fragment_prefix}
    """


rule write_organelle_curated_fasta:
    input:
        fa=os.path.join(config["outmain"], "{sample}", "{hap}", "harmonized", "{sample}.{hap}.harmonized.fa"),
        remove_list=os.path.join(config["outmain"], "{sample}", "{hap}", "organelle_screen", "{sample}.{hap}.contigs_to_remove.txt"),
    output:
        fa=os.path.join(config["outmain"], "{sample}", "{hap}", "final", "{sample}.{hap}.nuclear.no_organelle.fa"),
    params:
        script=os.path.join(workflow.basedir, "scripts", "filter_fasta_by_remove_list.py"),
        mode=lambda wc: config.get("organelle_screen", {}).get("mode", "remove_contigs"),
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
            --mode {params.mode}
    """


rule release_hap1_with_organelle:
    input:
        nuclear=os.path.join(config["outmain"], "{sample}", "hap1", "final", "{sample}.hap1.nuclear.clean.fa"),
        organelle=os.path.join(config["outmain"], "{sample}", "organelle", "{sample}.organelle.fa"),
    output:
        fa=os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap1_plus_organelle.fa"),
    resources:
        mem_mb=1000,
        runtime=30,
    shell: """
        cat {input.nuclear} {input.organelle} > {output.fa}
    """


rule release_hap2_nuclear:
    input:
        nuclear=os.path.join(config["outmain"], "{sample}", "hap2", "final", "{sample}.hap2.nuclear.clean.fa"),
    output:
        fa=os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap2_nuclear.fa"),
    resources:
        mem_mb=1000,
        runtime=20,
    shell: """
        ln -sfn $(realpath {input.nuclear}) {output.fa}
    """
