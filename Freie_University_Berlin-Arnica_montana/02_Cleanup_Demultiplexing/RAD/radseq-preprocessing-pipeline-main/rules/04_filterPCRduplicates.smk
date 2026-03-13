rule removePCRDuplicates:
	input:
		renamedDone = os.path.join(config['Results'], runNumber + "_configuration/renameNonConcat.done"),
		R1 = os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{Library}_{inID}_" + runNumber + "_1.fastq"),
		R2 = os.path.join(config['Results'], runNumber + "_0_PreProcessing/02_individualDemultiplex/beforeConcat/{Library}_{inID}_" + runNumber + "_2.fastq")
	params:
		script = os.path.join(workflow.basedir, "scripts/filterPCRdups/filterPCRdups_CM.py"),
		outpath = os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/beforeConcat/")
	output:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/beforeConcat/{Library}_{inID}_" + runNumber + "_1" + filenaming + "_CLUSTERS.fastq"),
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/beforeConcat/{Library}_{inID}_" + runNumber + "_1" + filenaming + ".fastq"),
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/beforeConcat/{Library}_{inID}_" + runNumber + "_2" + filenaming + "_CLUSTERS.fastq"),
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/beforeConcat/{Library}_{inID}_" + runNumber + "_2" + filenaming + ".fastq"),
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/beforeConcat/{Library}_{inID}_" + runNumber + "_filterPCRdups_Stats.tsv")
	log:
		os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/logs/{Library}_{inID}_filterPCRdups.log")
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_2.yaml")
	shell:
		"(python {params.script} -r1 {input.R1} -r2 {input.R2} -o {params.outpath}) &> {log}"

rule move03a_filterDups:
	input:
		check1= expand(os.path.join(config['Results'], runNumber + "_configuration/{indiID}_renameR1.done"), indiID=inID),
		check2= expand(os.path.join(config['Results'], runNumber + "_configuration/{indiID}_renameR2.done"), indiID=inID)
	params:
		clustersFile=os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/beforeConcat/*CLUSTERS.fastq"),
		statsFile=os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/beforeConcat/*Stats.tsv"),
		clustersFolder=os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/nonConcatCLUSTERS/"),
		statsFolder=os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/nonConcatStats/"),
		rmBeforeConcat=os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/beforeConcat/")
	output:
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/nonConcatCLUSTERS/{Librarys}_{indiID}_" + runNumber + "_1" + filenaming + "_CLUSTERS.fastq.gz"), Librarys=Library, indiID=inID),
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/nonConcatCLUSTERS/{Librarys}_{indiID}_" + runNumber + "_2" + filenaming + "_CLUSTERS.fastq.gz"), Librarys=Library, indiID=inID),
		expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/03_filterPCRduplicates/nonConcatStats/{Librarys}_{indiID}_" + runNumber + "_filterPCRdups_Stats.tsv"), Librarys=Library, indiID=inID)
	threads:
		16
	conda:
		os.path.join(workflow.basedir, "envs/preProcessing_3.yaml")
	shell:
		"""
		mv {params.clustersFile} {params.clustersFolder}
		mv {params.statsFile} {params.statsFolder}
		rm -r {params.rmBeforeConcat}
		pigz -p {threads} -r {params.clustersFolder}
		"""
