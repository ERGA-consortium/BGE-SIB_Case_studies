
rule compress01:
	input:
		check1= expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/logs/{Librarys}_RunFlexbar.log"), Librarys=Library),
		check2 =expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Librarys}_phiX_bowtie2.sam"), Librarys=Library),
		check3 =expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Librarys}_noPhiX" + forwardSuffix + ".fastq"), Librarys=Library),
		check4 =expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Librarys}_noPhiX" + reverseSuffix + ".fastq"), Librarys=Library)
	params:
		fastq = os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/")
	output:
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Librarys}_phiX_bowtie2.sam.gz"), Librarys=Library),
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Librarys}_noPhiX" + forwardSuffix + ".fastq.gz"), Librarys=Library),
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Librarys}_noPhiX" + reverseSuffix + ".fastq.gz"), Librarys=Library)
	threads:
		16
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_3.yaml")
	shell:
		"""
		pigz -p {threads} -r {params.fastq}
		"""

rule compress03b:
	input:
		check = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/logs/{indiID}_pearMerge.log"), indiID=inID),
		three = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{Librarys}_{indiID}_" + runNumber + "_1.fastq"), Librarys=Library,indiID=inID),
		four = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{Librarys}_{indiID}_" + runNumber + "_2.fastq"), Librarys=Library, indiID=inID)
	params:
		zip = os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/")
	output:
		three = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{Librarys}_{indiID}_" + runNumber + "_1.fastq.gz"), Librarys=Library,indiID=inID),
		four = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{Librarys}_{indiID}_" + runNumber + "_2.fastq.gz"), Librarys=Library, indiID=inID)
	threads:
		16
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_3.yaml")
	shell:
		"""
		pigz -p {threads} -r {params.zip}
		"""

rule compress04:
	input:
		check1= expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/logs/{indiID}_runDigest-" + EnzymeStringUnderFixed + ".log"), indiID=inID),
		check2 =expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/Lengths/{indiID}_" + runNumber + filenaming + ".Lengths.svg"), indiID=inID)
	params:
		fastq = os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/fastq/")
	output:
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/fastq/{indiID}_" + runNumber + filenaming + ".assembled.fastq.gz"), indiID=inID),
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/fastq/{indiID}_" + runNumber + filenaming + ".discarded.fastq.gz"), indiID=inID),
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/fastq/{indiID}_" + runNumber + filenaming + ".unassembled.forward.fastq.gz"), indiID=inID),
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/fastq/{indiID}_" + runNumber + filenaming + ".unassembled.reverse.fastq.gz"), indiID=inID)
	threads:
		16
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_3.yaml")
	shell:
		"""
		pigz -p {threads} -r {params.fastq}
		"""

rule compress05:
	input:
		check1= expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/fragmentSizes/{indiID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderhyphen + "_fragmentSizes.tsv"), indiID=inID),
		check2=expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/stats/{indiID}_" + runNumber + filenaming + ".assembled_digestion_statistics.tsv"), indiID=inID),
		check = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[2] + "/logs/{indiID}_checkRestriction-" + EnzymeStringUnderFixed + ".log"), indiID=inID)
	params:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/fastq/")
	output:
		fastq = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/fastq/{indiID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + ".fastq.gz"), indiID=inID)
	threads:
		16
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_3.yaml")
	shell:
		"""
		pigz -p {threads} -r {params}
		"""

rule compress06:
	input:
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[3] + "/logs/{indiID}_trimmomaticQ" + trimmomaticQvalue + "-" + EnzymeStringUnderFixed + ".log"), indiID=inID),
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[2] + "/siteErrors/{indiID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteErrors.tsv"), indiID=inID)
	params:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[2] + "/fastq/")
	output:
		fastq = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[2] + "/fastq/{indiID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteFiltered.fastq.gz"), indiID=inID)
	threads:
		16
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_3.yaml")
	shell:
		"""
		pigz -p {threads} -r {params}
		"""

rule compress07_TrimmedFiltered:
	input:
		check = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[3] + "/Lengths/{indiID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteFiltered.Q" + trimmomaticQvalue + ".Lengths.svg"), indiID=inID)
	params:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[3] + "/fastq/")
	output:
		fastq = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[3] + "/fastq/{indiID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteFiltered.Q" + trimmomaticQvalue + ".fastq.gz"), indiID=inID)
	threads:
		16
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_3.yaml")
	shell:
		"""
		pigz -p {threads} -r {params}
		"""

rule doneAll:
	input:
		three = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{Librarys}_{indiID}_" + runNumber + "_1.fastq.gz"), Librarys=Library,indiID=inID),
		four = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{Librarys}_{indiID}_" + runNumber + "_2.fastq.gz"), Librarys=Library, indiID=inID),
		compress01_a=expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Librarys}_phiX_bowtie2.sam.gz"), Librarys=Library),
		compress01_b=expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Librarys}_noPhiX" + forwardSuffix + ".fastq.gz"), Librarys=Library),
		compress01_c=expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/01_filterPhix/{Librarys}_noPhiX" + reverseSuffix + ".fastq.gz"), Librarys=Library),
		mvPcr_a=expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/nonConcatCLUSTERS/{Librarys}_{indiID}_" + runNumber + "_1" + filenaming + "_CLUSTERS.fastq.gz"), Librarys=Library, indiID=inID),
		mvPcr_b=expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/nonConcatCLUSTERS/{Librarys}_{indiID}_" + runNumber + "_2" + filenaming + "_CLUSTERS.fastq.gz"), Librarys=Library, indiID=inID),
		mvPcr_c=expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/nonConcatStats/{Librarys}_{indiID}_" + runNumber + "_filterPCRdups_Stats.tsv"), Librarys=Library, indiID=inID),
		compress04_1=expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/fastq/{indiID}_" + runNumber + filenaming + ".assembled.fastq.gz"), indiID=inID),
		compress04_2=expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/fastq/{indiID}_" + runNumber + filenaming + ".discarded.fastq.gz"), indiID=inID),
		compress04_3=expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/fastq/{indiID}_" + runNumber + filenaming + ".unassembled.forward.fastq.gz"), indiID=inID),
		compress04_4=expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[0] + "/fastq/{indiID}_" + runNumber + filenaming + ".unassembled.reverse.fastq.gz"), indiID=inID),
		ok = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[3] + "/Lengths/{indiID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteFiltered.Q" + trimmomaticQvalue + ".Lengths.svg"), indiID=inID),
		fastq1 = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[2] + "/fastq/{indiID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteFiltered.fastq.gz"), indiID=inID),
		fastq2 = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[1] + "/fastq/{indiID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + ".fastq.gz"), indiID=inID),
		fastq = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[3] + "/fastq/{indiID}_" + runNumber + filenaming + ".assembled_" + EnzymeStringUnderFixed + "_siteFiltered.Q" + trimmomaticQvalue + ".fastq.gz"), indiID=inID)
	output:
		temp(os.path.join(config['Results'], runNumber + "_configuration/preProcessing.done"))
	shell:
		"touch {output} | echo Execution is complete! thank you for using this Pre-Processing pipeline"


rule moveFinal:
	input:
		os.path.join(config['Results'], runNumber + "_configuration/preProcessing.done")
	params:
		configOG="configuration/config.yaml",
		FlexbarOG="configuration/FlexbarSheet.tsv"
	output:
		configNEW= os.path.join(config['Results'], runNumber + "_configuration/config.yaml"),
		FlexbarNEW= os.path.join(config['Results'], runNumber + "_configuration/FlexbarSheet.tsv")
	shell:
		"""
		cp {params.configOG} {output.configNEW}
		cp {params.FlexbarOG} {output.FlexbarNEW}
		"""
