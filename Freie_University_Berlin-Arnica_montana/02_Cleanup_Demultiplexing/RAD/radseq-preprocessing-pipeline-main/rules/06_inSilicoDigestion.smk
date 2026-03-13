rule runDigest:
	input:
		files = os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/fastq/{inID}_" + runNumber + filenaming + ".assembled.fastq")
	params:
		enzyme = EnzString,
		script = os.path.join(workflow.basedir, "scripts/digestion/RAD_digestion_v2.0.py"),
		outpath = os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/fastq/{inID}_" + runNumber + filenaming + ".assembled")
	output:
		fastq = os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/fastq/{inID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + ".fastq")
	log:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/logs/{inID}_runDigest-" + EnzymeStringUnderFixed + ".log")
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_2.yaml")
	shell:
		"(python {params.script} --g {input.files} --e {params.enzyme} --o {params.outpath} --dd --fq --rad) 2> {log}"


rule moveFilesRunDigest:
	input:
		check1= expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/logs/{indiID}_runDigest-" + EnzymeStringUnderFixed + ".log"), indiID=inID)
	params:
		fragmentFile=os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/fastq/*fragmentSizes.tsv"),
		statsFile=os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/fastq/*statistics.tsv"),
		fragmentFolder=os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/fragmentSizes/"),
		statsFolder=os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/stats/")
	output:
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/fragmentSizes/{indiID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderhyphen + "_fragmentSizes.tsv"), indiID=inID),
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/stats/{indiID}_" + runNumber + filenaming + ".assembled_digestion_statistics.tsv"), indiID=inID)
	shell:
		"""
		mv {params.fragmentFile} {params.fragmentFolder}
		mv {params.statsFile} {params.statsFolder}
		"""
