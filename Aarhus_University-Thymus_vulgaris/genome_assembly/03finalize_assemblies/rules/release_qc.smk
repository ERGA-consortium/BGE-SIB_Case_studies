rule do_release_qc_all:
    input:
        expand(
            os.path.join(config["outmain"], "{sample}", "release_qc", "stats", "{sample}.{release}.fa.stats"),
            sample=config["samples"].keys(),
            release=["hap1_plus_organelle", "hap2_nuclear"],
        ),
        expand(
            os.path.join(config["outmain"], "{sample}", "release_qc", "busco", "{sample}.{release}.short_summary.txt"),
            sample=config["samples"].keys(),
            release=["hap1_plus_organelle", "hap2_nuclear"],
        ),
        expand(
            os.path.join(config["outmain"], "{sample}", "release_qc", "merqury", "{sample}.release.qv"),
            sample=config["samples"].keys(),
        ),
        expand(
            os.path.join(config["outmain"], "{sample}", "release_qc", "merqury", "{sample}.release.completeness.stats"),
            sample=config["samples"].keys(),
        ),


rule release_fasta_stats:
    wildcard_constraints:
        release="hap1_plus_organelle|hap2_nuclear"
    input:
        fa=os.path.join(config["outmain"], "{sample}", "release", "{sample}.{release}.fa"),
    output:
        stats=os.path.join(config["outmain"], "{sample}", "release_qc", "stats", "{sample}.{release}.fa.stats"),
    threads: 3
    resources:
        mem_mb=5000,
        runtime=100,
        cpus_per_task=3,
    conda:
        "../envs/gfastats.yml"
    shell: """
        gfastats {input.fa} > {output.stats}
    """


rule busco_release:
    wildcard_constraints:
        release="hap1_plus_organelle|hap2_nuclear"
    input:
        fa=os.path.join(config["outmain"], "{sample}", "release", "{sample}.{release}.fa"),
    output:
        summary=os.path.join(config["outmain"], "{sample}", "release_qc", "busco", "{sample}.{release}.short_summary.txt"),
    params:
        busco_db=config["busco_db"],
        out_path=os.path.join(config["outmain"], "{sample}", "release_qc", "busco"),
        out_name=lambda wc: f"{wc.sample}.{wc.release}",
    resources:
        mem_mb=50000,
        runtime=1200,
        cpus_per_task=10,
    threads: 10
    conda:
        "../envs/busco.yml"
    shell: """
        busco --in {input.fa} \
              --out {params.out_name} \
              --out_path {params.out_path} \
              --mode genome \
              --lineage_dataset {params.busco_db} \
              --cpu {threads}
        cp {params.out_path}/{params.out_name}/short_summary*.txt {output.summary}
    """


rule merqury_release:
    input:
        meryl_db=lambda wc: config["release_qc"]["merqury"]["meryl_db"][str(wc.sample)],
        hap1=os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap1_plus_organelle.fa"),
        hap2=os.path.join(config["outmain"], "{sample}", "release", "{sample}.hap2_nuclear.fa"),
    output:
        qv=os.path.join(config["outmain"], "{sample}", "release_qc", "merqury", "{sample}.release.qv"),
        completeness=os.path.join(config["outmain"], "{sample}", "release_qc", "merqury", "{sample}.release.completeness.stats"),
    params:
        outdir=os.path.join(config["outmain"], "{sample}", "release_qc", "merqury"),
        outprefix=lambda wc: f"{wc.sample}.release",
    conda:
        "../envs/merqury.yml"
    resources:
        mem_mb=11000,
        runtime=400,
        cpus_per_task=10,
    shell: """
        meryl_db=`realpath {input.meryl_db}`
        hap1=`realpath {input.hap1}`
        hap2=`realpath {input.hap2}`
        cd {params.outdir}
        merqury.sh $meryl_db $hap1 $hap2 {params.outprefix}
    """
