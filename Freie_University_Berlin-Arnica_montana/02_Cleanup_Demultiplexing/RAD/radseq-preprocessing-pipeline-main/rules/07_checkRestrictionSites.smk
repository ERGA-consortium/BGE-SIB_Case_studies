rule checkRestriction:
	input:
		fastqIn = os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/fastq/{inID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + ".fastq")
	params:
		script = os.path.join(workflow.basedir, "scripts/digestion/checkRestrictionSites.py"),
		firstSequence = firstEnzyme,
		secondSequence = secondEnzyme
	output:
		filtered = os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[2] + "/fastq/{inID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteFiltered.fastq"),
		errors = os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[2] + "/fastq/{inID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteErrors.tsv")
	log:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[2] + "/logs/{inID}_checkRestriction-" + EnzymeStringUnderFixed + ".log")
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_2.yaml")
	shell:
		"(python {params.script} {input.fastqIn} {output.filtered} {output.errors} {params.firstSequence} {params.secondSequence}) 2> {log}"

rule moveFilescheckRestriction:
	input:
		check1= expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[2] + "/logs/{indiID}_checkRestriction-" + EnzymeStringUnderFixed + ".log"), indiID=inID)
	params:
		errorsFile=os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[2] + "/fastq/*siteErrors.tsv"),
		errorsFolder=os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[2] + "/siteErrors/")
	output:
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[2] + "/siteErrors/{indiID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteErrors.tsv"), indiID=inID),
	shell:
		"""
		mv {params.errorsFile} {params.errorsFolder}
		"""
