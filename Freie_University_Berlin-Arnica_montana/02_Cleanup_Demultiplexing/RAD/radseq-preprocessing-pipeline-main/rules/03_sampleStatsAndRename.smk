rule getReadCounts1:
	input:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/logs/{Library}_RunFlexbar.log")
	output:
		counts = temp(os.path.join(config['Results'], runNumber + "_configuration/{Library}_read_count.tsv")),
		totalProcessed = temp(os.path.join(config['Results'], runNumber + "_configuration/{Library}_totalproccessed.tsv"))
	shell:
		"""
		cat {input} | awk '/written reads/ {{ print $3 }}' | awk '{{ if(NR%2==1) print $1 }}' > {output.counts}
		cat {input} | awk '/Processed reads/ {{ print $3 }}' > {output.totalProcessed}
		cat {input} | awk '/Discarded reads overall/ {{ print $4 }}' >> {output.totalProcessed}
		cat {input} | awk '/Remaining reads/ {{ print $3 }}' >> {output.totalProcessed}
		"""


rule cleanBarcodeList:
	input:
		bcodez = os.path.join(workflow.basedir, "configuration/barcodesList.tsv")
	output:
		temp(os.path.join(config['Results'], runNumber + "_configuration/barcodesList2Use.tsv"))
	shell:
		"""
		sed 's/\r//g' < {input} >  {output}
		rm {input}
		"""


rule sampleAssignmentsPerLibraryColumn1:
	input:
		"configuration/FlexbarSheet.tsv"
	output:
		temp(os.path.join(config['Results'], runNumber + "_configuration/{Library}_assignedTogether.tsv"))
	shell:
		"cat {input} | awk '/{wildcards.Library}/,/^\s*$/' | awk 'FNR>2 {{ print $0 }}' | sed '/^$/q' > {output}"


rule myPandasScript:
	input:
		os.path.join(config['Results'], runNumber + "_configuration/{Library}_read_count.tsv"),
		os.path.join(config['Results'], runNumber + "_configuration/{Library}_totalproccessed.tsv"),
		os.path.join(workflow.basedir, "configuration/renamingFileLocation.tsv"),
		os.path.join(config['Results'], runNumber + "_configuration/{Library}_assignedTogether.tsv"),
		os.path.join(config['Results'], runNumber + "_configuration/barcodesList2Use.tsv")
	output:
		temp(os.path.join(config['Results'], runNumber + "_configuration/{Library}_renameR1.tsv")),
		temp(os.path.join(config['Results'], runNumber + "_configuration/{Library}_renameR2.tsv")),
		temp(os.path.join(config['Results'], runNumber + "_configuration/{Library}_FlexbarWithCountsAndLibraryName.tsv"))
	script:
		os.path.join(workflow.basedir, "scripts/pandass.py")


rule concatFlexbarSheetsWithCounts:
	input:
		individualSheets = expand(os.path.join(config['Results'], runNumber + "_configuration/{Librarys}_FlexbarWithCountsAndLibraryName.tsv"), Librarys=Library),
		oldSheet = "configuration/FlexbarSheet.tsv"
	output:
		temp(os.path.join(config['Results'], runNumber + "_configuration/FlexbarSheetWithCountsNoRun.tsv"))
	shell:
		"awk 'FNR==1{{print \"\"}}1' {input.individualSheets} > {output}"


rule addRunNumber:
	input:
		rules.concatFlexbarSheetsWithCounts.output
	params:
		run = runNumber
	output:
		os.path.join(config['Results'], runNumber + "_configuration/FlexbarSheetWithCounts.tsv")
	shell:
		"sed '1s/^/#{params.run}\\n/' {input} > {output}"


rule renameR1:
	input:
		check2 = os.path.join(config['Results'], runNumber + "_configuration/FlexbarSheetWithCounts.tsv"),
		check = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{{Library}}_barcode_{barcodez}_1.fastq"), barcodez=barcodeID),
		file = expand(os.path.join(config['Results'], runNumber + "_configuration/{{Library}}_renameR1.tsv"))
	output:
		one = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{{Library}}_{indiID}_" + runNumber + "_1.fastq"), indiID=inID),
		two = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{{Library}}_none_" + runNumber + "_1.fastq"))
	shell:
		"touch {output.one} {output.two} | awk -F'\t' 'system(\"mv \" $1 \" \" $2)' {input.file}"


rule renameR2:
	input:
		check2 = os.path.join(config['Results'], runNumber + "_configuration/FlexbarSheetWithCounts.tsv"),
		check = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{{Library}}_barcode_{barcodez}_2.fastq"), barcodez=barcodeID),
		file = expand(os.path.join(config['Results'], runNumber + "_configuration/{{Library}}_renameR2.tsv"))
	output:
		one = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{{Library}}_{indiID}_" + runNumber + "_2.fastq"), indiID=inID),
		two = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{{Library}}_none_" + runNumber + "_2.fastq"))
	shell:
		"touch {output.one} {output.two} | awk -F'\t' 'system(\"mv \" $1 \" \" $2)' {input.file}"


rule rmNoneRenameDone:
	input:
		renameR1Done = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{Librarys}_{indiID}_" + runNumber + "_1.fastq"), Librarys=Library, indiID=inID),
		renameR2Done = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{Librarys}_{indiID}_" + runNumber + "_2.fastq"), Librarys=Library, indiID=inID),
		remove1 = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{Librarys}_none_" + runNumber + "_1.fastq"), Librarys=Library),
		remove2 = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{Librarys}_none_" + runNumber + "_2.fastq"), Librarys=Library),
		remove3 = os.path.join(config['Results'], runNumber + "_configuration/barcodesList2Use.tsv"),
		remove4 = os.path.join(workflow.basedir, "configuration/renamingFileLocation.tsv"),
	output:
		temp(os.path.join(config['Results'], runNumber + "_configuration/renameNonConcat.done"))
	shell:
		"rm {input.remove1} {input.remove2} {input.remove3} {input.remove4} | touch {output}"
