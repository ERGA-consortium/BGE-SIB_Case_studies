
############################################
############ RULES TO CALL STUFF ###########
############################################

rule do_models:
    input:
        os.path.join(OUTMAIN, "repeatmodels", DB_NAME + "-families.fa")


rule do_masks:
    input:
        expand(
            os.path.join(OUTMAIN, "repeatmasks", DB_NAME, "{genome}", "{hap}", "{genome}.{hap}.fasta.masked"),
            zip,
            genome=[g for g, h in GENOME_HAPS],
            hap=[h for g, h in GENOME_HAPS],
        )


rule do_release_all:
    input:
        expand(
            os.path.join(OUTMAIN, "release", "{genome}", "{hap}", "{genome}.{hap}.softmasked.fa.fai"),
            zip,
            genome=[g for g, h in GENOME_HAPS],
            hap=[h for g, h in GENOME_HAPS],
        ),
        expand(
            os.path.join(OUTMAIN, "release", "{genome}", "{hap}", "{genome}.{hap}.repeats.bed"),
            zip,
            genome=[g for g, h in GENOME_HAPS],
            hap=[h for g, h in GENOME_HAPS],
        )


###########################################
############ RULES TO RUN STUFF ###########
###########################################

rule build_database:
    input:
        staged_fastas = rules.stage_genomes.input,
    output:
        done = os.path.join(OUTMAIN, "repeatmodels", "done_builddb.txt")
    params:
        work_dir = os.path.join(OUTMAIN, "repeatmodels"),
        db_name = DB_NAME,
        staged_dir = os.path.join(OUTMAIN, "inputs", "genomes"),
    resources:
        mem_mb = 50000,
        runtime = 600,
    conda:
        "../envs/repeatmodeler.yml"
    shell: """
        mkdir -p {params.work_dir}
        staged_dir=`realpath {params.staged_dir}`
        cd {params.work_dir}
        BuildDatabase -name {params.db_name} -dir $staged_dir
        outdone=`basename {output.done}`
        echo "done building database" > ${{outdone}}
    """


rule model_repeats:
    input:
        done = os.path.join(OUTMAIN, "repeatmodels", "done_builddb.txt")
    output:
        families_fa = os.path.join(OUTMAIN, "repeatmodels", DB_NAME + "-families.fa")
    params:
        work_dir = os.path.join(OUTMAIN, "repeatmodels"),
        db_name = DB_NAME,
        flags = " ".join([
            flag_if(REPEATMODELER, "ltrstruct", "-LTRStruct", default=True),
            extra_flags(REPEATMODELER),
        ]).strip(),
    threads: int(REPEATMODELER.get("threads", 20))
    resources:
        mem_mb = 50000,
        runtime = 4000,
        cpus_per_task = int(REPEATMODELER.get("threads", 20)),
    conda:
        "../envs/repeatmodeler.yml"
    shell: """
        cd {params.work_dir}
        RepeatModeler -database {params.db_name} -threads {threads} {params.flags}
    """


rule repeat_masker:
    input:
        fa = os.path.join(OUTMAIN, "inputs", "genomes", "{genome}.{hap}.fasta"),
        families_fa = os.path.join(OUTMAIN, "repeatmodels", DB_NAME + "-families.fa"),
    output:
        masked = os.path.join(OUTMAIN, "repeatmasks", DB_NAME, "{genome}", "{hap}", "{genome}.{hap}.fasta.masked"),
        rm_out = os.path.join(OUTMAIN, "repeatmasks", DB_NAME, "{genome}", "{hap}", "{genome}.{hap}.fasta.out"),
    params:
        outdir = os.path.join(OUTMAIN, "repeatmasks", DB_NAME, "{genome}", "{hap}"),
        flags = " ".join([
            flag_if(REPEATMASKER, "sensitive", "-s", default=True),
            flag_if(REPEATMASKER, "gff", "-gff", default=True),
            flag_if(REPEATMASKER, "xsmall", "-xsmall", default=True),
            extra_flags(REPEATMASKER),
        ]).strip(),
    threads: int(REPEATMASKER.get("threads", 20))
    resources:
        mem_mb = 50000,
        runtime = 4000,
        cpus_per_task = int(REPEATMASKER.get("threads", 20)),
    conda:
        "../envs/repeatmasker.yml"
    shell: """
        RepeatMasker -lib {input.families_fa} -pa {threads} -dir {params.outdir} {params.flags} {input.fa}
    """


rule release_softmasked_fasta:
    input:
        masked = os.path.join(OUTMAIN, "repeatmasks", DB_NAME, "{genome}", "{hap}", "{genome}.{hap}.fasta.masked"),
    output:
        fa = os.path.join(OUTMAIN, "release", "{genome}", "{hap}", "{genome}.{hap}.softmasked.fa"),
    resources:
        mem_mb=1000,
        runtime=20,
    shell: """
        ln -sfn $(realpath {input.masked}) {output.fa}
    """


rule index_release_softmasked_fasta:
    input:
        fa = os.path.join(OUTMAIN, "release", "{genome}", "{hap}", "{genome}.{hap}.softmasked.fa"),
    output:
        fai = os.path.join(OUTMAIN, "release", "{genome}", "{hap}", "{genome}.{hap}.softmasked.fa.fai"),
    conda:
        "../envs/samtools.yml"
    resources:
        mem_mb=2000,
        runtime=60,
    shell:
        "samtools faidx {input.fa}"


rule repeatmasker_out_to_bed:
    input:
        rm_out = os.path.join(OUTMAIN, "repeatmasks", DB_NAME, "{genome}", "{hap}", "{genome}.{hap}.fasta.out"),
    output:
        bed = os.path.join(OUTMAIN, "release", "{genome}", "{hap}", "{genome}.{hap}.repeats.bed"),
    resources:
        mem_mb=1000,
        runtime=20,
    shell: """
        awk 'BEGIN{{OFS="\\t"}} $1 ~ /^[0-9]+$/ {{s=$6-1; if (s<0) s=0; print $5, s, $7, $11}}' {input.rm_out} > {output.bed}
    """
