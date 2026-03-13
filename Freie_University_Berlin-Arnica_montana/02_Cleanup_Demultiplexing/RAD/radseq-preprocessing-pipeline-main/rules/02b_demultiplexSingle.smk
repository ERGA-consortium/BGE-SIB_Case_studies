rule demultiplexSingleBarcode:
	input:
		R1 = os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Library}_noPhiX" + forwardSuffix + ".fastq"),
		R2 = os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Library}_noPhiX" + reverseSuffix + ".fastq")
	params:
		barcodes = flexbarFileOne
		outPath = os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{Library}")
	output:
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{{Library}}_barcode_{barcodez}_1.fastq"), barcodez=barcodeID),
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{{Library}}_barcode_{barcodez}_2.fastq"), barcodez=barcodeID)
	log:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/logs/{Library}_RunFlexbar.log")
	threads:
		config['DemultiplexThreads']
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_1.yaml")
	shell:
		"flexbar --threads {threads} -t {params.outPath} -r {input.R1} -p {input.R2} -b {params.barcodes} -bt LTAIL -u 3 -O {log}"
