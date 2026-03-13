rule pear:
	input:
		concat1=os.path.join(config['Results'], runNumber + "_configuration/{inID}_renameR1.done"),
		concat2=os.path.join(config['Results'], runNumber + "_configuration/{inID}_renameR2.done"),
		read1 = os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[4] + "/concat/{inID}_" + runNumber + "_1_conc" + filenaming + ".fastq"),
		read2 = os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[4] + "/concat/{inID}_" + runNumber + "_2_conc" + filenaming + ".fastq")
	params:
		prefix = os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/fastq/{inID}_" + runNumber + filenaming)
	output:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/fastq/{inID}_" + runNumber + filenaming + ".assembled.fastq")
	log:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/logs/{inID}_pearMerge.log")
	threads:
		config['PEARthreads']
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_1.yaml")
	shell:
		"(pear -j {threads} -n 30 -m 600 -v 30 -f {input.read1} -r {input.read2} -o {params.prefix}) &> {log} "
