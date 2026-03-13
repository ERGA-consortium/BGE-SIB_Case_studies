#!/usr/bin/python

##################################################################################
#### property of BeGenDiv(https://begendiv.de/)								 #####
#### @authors: Harald Detering, Marie Jeschek, Maximilian Driller            #####
#### contact: drillermax@gmail.com                                           #####
##################################################################################


import argparse, os, sys, re
from Bio import SeqIO
from Bio.Seq import Seq

def mkOutDir(path):
	try:
		os.makedirs(os.path.dirname(path))
	except OSError:
		print("Output directory already exists. Already existing files will be overwritten.")
	if os.path.basename(path):
		return(path+"_")
	else: return(path)

def parseEnzymeList(eFile):
	enzymes = []	
	for line in eFile:
		data = line.strip().split("\t")
		thisEnzymes = []
		
		for i in range(0,len(data),2):
			enzymeName=data[i]
			try:
				enzymePattern=data[i+1]
				thisEnzymes.append(myRestrictionEnzyme(enzymeName, enzymePattern))
			except IndexError:
				print("ERROR: Wrong format for enzyme file. No restriction pattern found for enzyme %s." % enzymeName)
				print("The current line will be ignored: %s" % line)
				thisEnzymes=None
				break
		
		if thisEnzymes: enzymes.append(thisEnzymes)
	
	return(enzymes)

class myRestrictionEnzyme:
	def __init__(self, name, pattern):
		self.name=name
		self.pattern=pattern
		self.parsePattern()
	def __repr__(self):
		return("%s: %s" % (self.name, self.pattern))
	def parsePattern(self):
		self.pPattern=self.pattern.replace('^', '').upper()
		self.pPattern=self.pPattern.replace('_', '').upper() #for snakemake preprocessing pipeline
		self.offset=self.pattern.find('^')
		if self.offset > len(self.pattern)/2: # needed to flip the offset to correct site
			self.offset = len(self.pattern)-1-self.offset

		#print(self.name, self.pPattern, self.offset)



# parse arguments
parser = argparse.ArgumentParser(description="In-silico RAD digestion.", formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--g", "-genome", required=True, metavar="genomeFile", type=str,
	help="Sequences to digest in fasta or fastq format (f.e. a genome file).")
parser.add_argument("--e", "-enzymes", required=True, metavar="enzymeFile", type=argparse.FileType('r'),
	help="List with enzymes/enzyme combinations including restriction pattern. Tab-separated.\nFormat:\nenzyme1\tpattern1\tenzyme2\tpattern2 -- the first enzyme is the forward enzyme, the second one in the list is reverse (format is as shown in the file downloaded by the get_rad_enzymes.py script)")
parser.add_argument("--dd", action='store_true',
	help="If argument is set, only fragments cut by two different enzymes will be returned.")
parser.add_argument("--o", "-outdir", default=sys.stdout, metavar="outputFile", type=mkOutDir,
	help="Output file prefix. [Default: stdout]")
parser.add_argument("--min", type=int,
	help="Minimum size for fragments.")
parser.add_argument("--max", type=int,
	help="Maximum size for fragments.")
parser.add_argument("--rad", "-radseq", action='store_true',
	help="Allow fragments to have Start(=first/forward enzyme in list) and End(=second/reverse enzyme in list) of a scaffold as a restriction side. Enable if digesting already digested sequences.")
parser.add_argument("--fq", "-fastq", action='store_true',
	help="enable if genome is in fastq format e.g. reads redigested")
args = parser.parse_args()


enzymeCombinations = parseEnzymeList(args.e) # parse enzyme list



# iterate over contigs (from genome or reads)
outFiles={}
#dictionary for number of fragments per enzyme combintaion
fragDict={}
#dictionary for number of fragments per length 
sizeDICT={}
#variable to keep track of the longest fragment --> needed for output later
maxFrag=0

if args.g.endswith(".gz"):
	import gzip
	contigParser=gzip.open(args.g, "rt")
else:
	contigParser=open(args.g)


if args.fq:
	infmt="fastq"
else:
	infmt="fasta"


for contig in SeqIO.parse(contigParser, infmt):
	for enzymeComb in enzymeCombinations: # for each enzyme combination:

		# find cutting sites of enzymes ...
		cuttingSites = {}
		
		enzyme1=enzymeComb[0]
		enzyme2=enzymeComb[1]

		#print(enzyme1.name, enzyme1.pPattern, enzyme1.offset)
		#print(enzyme2.name, enzyme2.pPattern, enzyme2.offset)

		#if args.rad: # ignore first and last bp in rad mode to avoid finding a restriction site there
	#		cutPositions_Enzyme1 = [(m.start()+enzyme1.offset, m.end()-enzyme1.offset, enzyme1.name) for m in re.finditer(enzyme1.pPattern, str(contig.seq[1:len(contig.seq)-2]).upper())] # cutting sites for enzyme1
#			cutPositions_Enzyme2 = [(m.start()+enzyme2.offset, m.end()-enzyme2.offset, enzyme2.name) for m in re.finditer(enzyme2.pPattern, str(contig.seq[1:len(contig.seq)-2]).upper())] # cutting sites for enzyme2
#		else:
		cutPositions_Enzyme1 = [(m.start()+enzyme1.offset, m.end()-enzyme1.offset, enzyme1.name) for m in re.finditer(enzyme1.pPattern, str(contig.seq).upper())] # cutting sites for enzyme1
		cutPositions_Enzyme2 = [(m.start()+enzyme2.offset, m.end()-enzyme2.offset, enzyme2.name) for m in re.finditer(enzyme2.pPattern, str(contig.seq).upper())] # cutting sites for enzyme2
		
		cutPositions_Joined = cutPositions_Enzyme1 + cutPositions_Enzyme2 # join cutting sites
		cutPositions_Joined_sorted = sorted(cutPositions_Joined, key=lambda tup:tup[0]) # sort list by cut index
		#print(cutPositions_Joined_sorted)



		if args.rad: #digest rad loci --> keep start and end e.g. add artifical start end cutting positions
			if cutPositions_Joined_sorted != []:
				if cutPositions_Joined_sorted[0][0]==0:
					cutPositions_Joined_sorted.pop(0)
			if cutPositions_Joined_sorted != []:
				if cutPositions_Joined_sorted[-1][1]==len(contig.seq):
					cutPositions_Joined_sorted.pop(-1)


			startPos = (0, 0, 'start')
			endPos = (len(contig.seq), len(contig.seq), 'end')
			cutPositions_Joined_sorted.insert(0, startPos)
			cutPositions_Joined_sorted.append(endPos)

		nFrags = 1
		for i in range(0, len(cutPositions_Joined_sorted)-1):
			curPos = cutPositions_Joined_sorted[i]
			nextPos = cutPositions_Joined_sorted[i+1]
			#print(curPos)
			#print(nextPos)

			if not args.dd or (curPos[2] != nextPos[2]): # if not dd rad mode otherwise needs to be cut by two different enzymes
				
				if args.rad and (curPos[2] == 'start'and nextPos[2] == enzyme1.name):
					continue
				elif args.rad and (curPos[2] == enzyme2.name and nextPos[2] == 'end'):
					continue

				newFrag = contig.seq[curPos[0]:nextPos[1]]
				#print(newFrag)
				if args.fq:
					quals = contig.letter_annotations["phred_quality"][curPos[0]:nextPos[1]]
					#print(quals)
				
				if curPos[2] == enzyme2.name: # case where fragment needs to be reverse complemented
					newFrag = newFrag.reverse_complement()
					if args.fq:
						quals.reverse() # reverse quality if seq is reverse complemented
					newName = "%s|frag_%i|%i-%i|%s-%s"%(contig.id,nFrags,nextPos[0],curPos[1],nextPos[2],curPos[2])
				else:
					newName = "%s|frag_%i|%i-%i|%s-%s"%(contig.id,nFrags,curPos[0],nextPos[1],curPos[2],nextPos[2])
				#print(newFrag)
	
				################################################################
				################################################################
				#### creating outputs ##########################################
				outFileKey = "_".join([x.name for x in enzymeComb])
				if outFileKey not in outFiles:
					try:
						if args.min and args.max:
							outFiles[outFileKey] = open(args.o+"%s_%i-%i.%s" % ("_".join([x.name for x in enzymeComb]), args.min, args.max, infmt),"w")
						elif args.min:
							outFiles[outFileKey] = open(args.o+"%s_%i-.%s" % ("_".join([x.name for x in enzymeComb]), args.min, infmt),"w")
						elif args.max:
							outFiles[outFileKey] = open(args.o+"%s_-%i.%s" % ("_".join([x.name for x in enzymeComb]), args.max, infmt),"w")
						else:
							outFiles[outFileKey] = open(args.o+"%s.%s" % ("_".join([x.name for x in enzymeComb]), infmt),"w")
					except AttributeError:
						outFiles[outFileKey] = sys.stdout
				################################################################
				################################################################
				################################################################

				removed = False

				fragLen = len(newFrag)
				
				#checking for 3rd enzyme
				if len(enzymeComb) > 2: #true if we have a third enzyme
					enzyme3 = enzymeComb[2]
					if newFrag.upper().find(enzyme3.pPattern)!=-1:
						removed = True

				if args.min: # checking min len
					if fragLen < args.min:
						removed=True

				if args.max: # checking min len
					if fragLen > args.max:
						removed=True				

				if removed == False:
					#################################################
					################ add statistics #################
					#################################################
					enzNames = enzymeComb[0].name +"-"+ enzymeComb[1].name
					if not enzNames in fragDict:
						fragDict[enzNames] = 1
					else:
						fragDict[enzNames] += 1

					if sizeDICT.get(enzNames, "no") == "no": # init
						sizeDICT[enzNames] = {}
					if not fragLen in sizeDICT[enzNames]: # store num fragments per length
						sizeDICT[enzNames][fragLen] = 1
					else:
						sizeDICT[enzNames][fragLen] += 1
					if fragLen > maxFrag:
						maxFrag = fragLen
					#################################################
					#################################################
					#################################################
					if args.fq:
						qual = ""
						for bqual in quals:
							qual += chr(bqual+33)
						outFiles[outFileKey].write("@%s\n%s\n+\n%s\n" % (newName,newFrag,qual))
					else:
						outFiles[outFileKey].write(">%s\n%s\n" % (newName,newFrag))

					nFrags += 1
#########################
##### write ouputs ######
#########################
try:
	#output a file containing the number of fragments for each enzyme combination
	overviewWriter = open(args.o + "digestion_statistics.tsv", "w")
	

	for key in fragDict:
		overviewWriter.write("%s\t%i\n" % (key, fragDict[key]))


	for eComb in sizeDICT:
		oWriter = open(args.o + "%s_fragmentSizes.tsv"%(eComb), "w")
		for i in range(0, maxFrag+1):
		
			if sizeDICT[eComb].get(i, "nope") == "nope":
				oWriter.write("%i\t%i\n"%(i, 0))
			else:
				oWriter.write("%i\t%i\n"%(i, sizeDICT[eComb][i]))
	
except AttributeError:
	overviewWriter = sys.stdout

finally:
	overviewWriter.close()
