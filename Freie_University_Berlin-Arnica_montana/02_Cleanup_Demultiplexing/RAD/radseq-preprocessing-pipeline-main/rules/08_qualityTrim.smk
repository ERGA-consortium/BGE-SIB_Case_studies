
rule trimmomatic:
	input:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[2] + "/fastq/{inID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteFiltered.fastq")
	output:
		fastq = os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[3] + "/fastq/{inID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteFiltered.Q" + trimmomaticQvalue + ".fastq")
	log:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[3] + "/logs/{inID}_trimmomaticQ" + trimmomaticQvalue + "-" + EnzymeStringUnderFixed + ".log")
	threads:
		16
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_1.yaml")
	shell:
		"(trimmomatic SE -threads {threads} -phred33 {input} {output.fastq} AVGQUAL:{trimQvalue}) 2> {log}"
