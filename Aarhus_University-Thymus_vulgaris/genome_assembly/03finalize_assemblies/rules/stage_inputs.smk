_STAGE_PAIRS = [
    (str(sample), str(hap))
    for sample, sample_cfg in config["samples"].items()
    for hap in sample_cfg["haps"].keys()
]
_STAGE_SAMPLES = [sample for sample, _ in _STAGE_PAIRS]
_STAGE_HAPS = [hap for _, hap in _STAGE_PAIRS]


rule do_stage_inputs:
    input:
        expand(
            os.path.join(config["outmain"], "input", "{sample}", "{hap}.fa"),
            zip,
            sample=_STAGE_SAMPLES,
            hap=_STAGE_HAPS,
        ),


rule stage_assembly:
    input:
        src=lambda wc: config["samples"][str(wc.sample)]["haps"][str(wc.hap)],
    output:
        dst=os.path.join(config["outmain"], "input", "{sample}", "{hap}.fa"),
    resources:
        mem_mb=1000,
        runtime=20,
    shell: """
        ln -sfn $(realpath {input.src}) {output.dst}
    """
