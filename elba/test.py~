from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Seq
import sys

def doBlast (fastafile):
	handle = open(fastafile)
	
	for record in SeqIO.parse(handle, "fasta") :
		result = NCBIWWW.qblast("blastp", "nr", record.seq)
		blastfile = open("blast.xml", "w")
		blastfile.write(result.read())
		blastfile.close()	
	handle.close()

def selectProt(blastxml):
	result = open(blastxml, "r")
	out = open("candidates.blast", "w")
	for n in NCBIXML.parse(result):
		for alignment in n.alignments:
			for hsp in alignment.hsps:
	#			if hsp.expect < 0.1:
					out.write("****Alignment****\n")
					out.write("Sequence: "+str(alignment.title)+"\n")
					out.write("Length: "+str(alignment.length)+"\n")
					out.write("E-value: "+str(hsp.expect)+"\n")
					out.write(str(hsp.query[:50]) + "...\n")
					out.write(str(hsp.match[:50]) + "...\n")
					out.write(str(hsp.sbjct[:50]) + "...\n")
	result.close()
#doBlast ("thrombin.fa")
selectProt("blast.xml")


#BLAST running locally --> output 
#from output blast: selectsequencies with less than 1e-5 e-value. Extract outputfile: id/seq/evalue/homology/coverage/
#filter outputfile: maxium 15 sequences present in both outputfiles and save each hit  in two files for each protein
#do clustalW for both files
#construct phylogenetic tree: tree representation  
#construct matrix
#r pearson correlation 





