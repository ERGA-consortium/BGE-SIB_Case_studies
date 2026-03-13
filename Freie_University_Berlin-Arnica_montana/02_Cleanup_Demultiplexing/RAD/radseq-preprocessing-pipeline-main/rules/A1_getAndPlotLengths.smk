rule getLengths:
	input:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/fastq/{inID}_" + runNumber + filenaming + ".assembled.fastq")
	output:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/Lengths/{inID}_" + runNumber + filenaming + ".Lengths")
	shell:
		"cat {input} | awk '{{if(NR%4==2) print length($1)}}' | sort -n | uniq -c > {output}"

rule plotLengths:
	input:
		rules.getLengths.output
	params:
		os.path.join(workflow.basedir, "scripts/plotFragSizes_args_fastq.R")
	output:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/Lengths/{inID}_" + runNumber + filenaming + ".Lengths.svg")
	shell:
		"Rscript {params} {input} {output}"


rule getLengths2:
	input:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[3] + "/fastq/{inID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteFiltered.Q" + trimmomaticQvalue + ".fastq")
	output:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[3] + "/Lengths/{inID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteFiltered.Q" + trimmomaticQvalue + ".Lengths")
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_1.yaml")
	shell:
		"cat {input} | awk '{{if(NR%4==2) print length($1)}}' | sort -n | uniq -c > {output}"

rule plotLengths2:
	input:
		rules.getLengths2.output
	params:
		script = os.path.join(workflow.basedir, "scripts/plotFragSizes_args_fastq.R")
	output:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[3] + "/Lengths/{inID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteFiltered.Q" + trimmomaticQvalue + ".Lengths.svg")
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_1.yaml")
	shell:
		"Rscript {params.script} {input} {output}"
