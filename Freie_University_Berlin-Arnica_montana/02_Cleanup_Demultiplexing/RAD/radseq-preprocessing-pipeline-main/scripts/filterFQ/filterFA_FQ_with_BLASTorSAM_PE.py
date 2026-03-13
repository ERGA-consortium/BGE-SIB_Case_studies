from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
import argparse


parser = argparse.ArgumentParser(description="")
parser.add_argument("r1", help="path to fastq fwd file as an input (can also be zipped)", type=str)
parser.add_argument("r2", help="path to fastq rev file as an input (can also be zipped)", type=str)
parser.add_argument("sam", help="sam file only unmapped reads will be written out", type=str)
parser.add_argument("out1", help="path and name of output file which will be created", type=str)
parser.add_argument("out2", help="path and name of output file which will be created", type=str)
parser.add_argument("--m", "-mapped", help="put put reads that MAPPED/BLASTED against something and dont output the ones with no hits", action="store_true")
parser.add_argument("--b", "-blast", help="if is done for a blast table not sam file BLAST ouput must be in format 6 or 7 default columns", action="store_true")
parser.add_argument("--fa", "-fasta", help="if file in fasta not fastq format", action="store_true")

args = parser.parse_args()


if args.r1.endswith(".gz"):
	try:
		import gzip 
		handle1 = gzip.open(args.r1, mode="rt")
	except IOError:
		print("Error - could not open %s" % (args.r1))

else:
	try:
		handle1 = open(args.r1, mode="rt")
	except IOError:
		print("Error - could not open %s" % (args.r1))

if args.r2.endswith(".gz"):
	try:
		import gzip 
		handle2 = gzip.open(args.r2, mode="rt")
	except IOError:
		print("Error - could not open %s" % (args.r2))

else:
	try:
		handle2 = open(args.r2, mode="rt")
	except IOError:
		print("Error - could not open %s" % (args.r2))



outhandle1 = open(args.out1, "w")
outhandle2 = open(args.out2, "w")
samReader=open(args.sam)


mappedDICT={}

for line in samReader:
	
	#checking SAM
	if not args.b:
		if not line.startswith("@"):
			splitted=line.strip().split("\t")


			rID = splitted[0].strip(";")

			if mappedDICT.get(rID, "nope") == "nope":
				mappedDICT[rID] = 1
			#else:
				#print("Weird rID %s is already in dictionary!!!"%(rID))
	

	#checking BLAST
	#can have more than one hit per read
	else:
		if not line.startswith("#"):
			splitted=line.strip().split("\t")
			rID = splitted[0].strip(";")
			#print(rID)			

			if mappedDICT.get(rID, "nope") == "nope":
				mappedDICT[rID] = 1
			#else:
			#	print("Weird rID %s is already in dictionary!!!"%(rID))



samReader.close()

seqsMapped=len(mappedDICT)
totalcount = 0
seqcount = 0


if args.fa:
	inFormat = "fasta"
else:
	inFormat = "fastq"

#import io
#seq_file = handle2
#byte_str = seq_file.read()
#text_obj = byte_str.decode('UTF-8')

r2Iterator=SeqIO.parse(handle2, inFormat)

for read1 in SeqIO.parse(handle1, inFormat):
	read2=next(r2Iterator)


	#print(read1.id)
	totalcount+=1
	
	if not args.m:

		curID1 = read1.id.strip(";")
		curID2 = read2.id.strip(";")
		#print(curID1)


		if mappedDICT.get(curID1, "nope") == "nope":
			SeqIO.write(read1, outhandle1, inFormat)
			SeqIO.write(read2, outhandle2, inFormat)
			seqcount+=1
		elif mappedDICT.get(curID1, "nope") == "nope":
			SeqIO.write(read1, outhandle1, inFormat)
			SeqIO.write(read2, outhandle2, inFormat)
			seqcount+=1
	#ouput reads in SAM/BLAST file
	else:
		curID1 = read1.id.strip(";")
		curID2 = read2.id.strip(";")
		#print(curID)

		if mappedDICT.get(curID1, "nope") != "nope":
			SeqIO.write(read1, outhandle1, inFormat)
			SeqIO.write(read2, outhandle2, inFormat)
			seqcount+=1
		elif mappedDICT.get(curID1, "nope") != "nope":
			SeqIO.write(read1, outhandle1, inFormat)
			SeqIO.write(read2, outhandle2, inFormat)
			seqcount+=1




print("All done %i of %i sequences survived number of seqs filtered out: %i." % (seqcount, totalcount, seqsMapped))

outhandle1.close()
outhandle2.close()
handle1.close()
handle2.close()
