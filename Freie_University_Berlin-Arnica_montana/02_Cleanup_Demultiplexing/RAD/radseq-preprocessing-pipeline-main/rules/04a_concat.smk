rule concateReplicatesR1:
	input:
		renameR1Done = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[4] + "/beforeConcat/{Librarys}_{{inID}}_" + runNumber + "_1" + filenaming + ".fastq"), Librarys=Library)
	output:
		concatenatedR1 = os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[4] + "/concat/{inID}_" + runNumber + "_1_conc" + filenaming + ".fastq"),
		renameDoneR1 = temp(os.path.join(config['Results'], runNumber + "_configuration/{inID}_renameR1.done"))
	shell:
		"touch {output.renameDoneR1} | cat {input.renameR1Done} > {output.concatenatedR1}"


rule concateReplicatesR2:
	input:
		renameR2Done = expand(os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[4] + "/beforeConcat/{Librarys}_{{inID}}_" + runNumber + "_2" + filenaming + ".fastq"), Librarys=Library),
	output:
		concatenatedR2 = os.path.join(config['Results'], runNumber + "_0_PreProcessing/" + folderNaming[4] + "/concat/{inID}_" + runNumber + "_2_conc" + filenaming + ".fastq"),
		renameDoneR2 = temp(os.path.join(config['Results'], runNumber + "_configuration/{inID}_renameR2.done"))
	shell:
		"touch {output.renameDoneR2} | cat {input.renameR2Done} > {output.concatenatedR2}"
