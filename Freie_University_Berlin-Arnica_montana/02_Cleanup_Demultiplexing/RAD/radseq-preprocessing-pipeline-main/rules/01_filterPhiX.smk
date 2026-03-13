print('')
print("This run is:", runNumber)
print('')
print("Your pooled Library names for this run are:", Library)
print('')
print("Your individual sample IDs for this run are:", inID)

forwardSuffix = config['suffixForwardRead']
reverseSuffix = config['suffixReverseRead']

FORMAT = ""
if config['gzipped'] in affirmativeUserResponse and config['fqORfastq'] == "fastq":
	FORMAT = "fastq.gz"
elif config['gzipped'] in affirmativeUserResponse and config['fqORfastq'] == "fq":
	FORMAT = "fq.gz"
elif config['gzipped'] not in affirmativeUserResponse and config['fqORfastq'] == "fastq":
	FORMAT = "fastq"
elif config['gzipped'] not in affirmativeUserResponse and config['fqORfastq'] == "fq":
	FORMAT = "fq"

rule bowtie2:
	input:
		ref = config["Phix_index"],
		read1 = os.path.join(config['RawDataFolder'], "{Library}" + forwardSuffix + "." + FORMAT),
		read2 = os.path.join(config['RawDataFolder'], "{Library}" + reverseSuffix + "." + FORMAT)
	output:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Library}_phiX_bowtie2.sam")
	log:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/logs/{Library}_bowtie2phiXmapping.log")
	threads:
		config['Bowtie2Threads']
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_1.yaml")
	shell:
		"(bowtie2 --no-unal -x {input.ref} -1 {input.read1} -2 {input.read2} -p {threads} -S {output}) 2> {log}"

rule filterUnmapped:
	input:
		sam = rules.bowtie2.output,
		read1 = os.path.join(config['RawDataFolder'], "{Library}" + forwardSuffix + "." + FORMAT),
		read2 = os.path.join(config['RawDataFolder'], "{Library}" + reverseSuffix + "." + FORMAT)
	params:
		script = os.path.join(workflow.basedir, "scripts/filterFQ/filterFA_FQ_with_BLASTorSAM_PE.py")
	output:
		R1out = os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Library}_noPhiX" + forwardSuffix + ".fastq"),
		R2out = os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Library}_noPhiX" + reverseSuffix + ".fastq")
	log:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/logs/{Library}_getReadsFromSAM.log")
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_2.yaml")
	shell:
		"(python {params.script} {input.read1} {input.read2} {input.sam} {output.R1out} {output.R2out}) &> {log}"
