from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Seq

def doBlast (fastafile):
	handle = open(fastafile)
	
	for record in SeqIO.parse(handle, "fasta") :
		result = NCBIWWW.qblast("blastp", "nr", record.seq)
		blastfile = open("blast.xml", "w")
		blastfile.write(result.read())
		blastfile.close()	
	handle.close()

result = open("blast.xml", "r")
for n in NCBIXML.parse(result):
	for alignment in n.alignments:
		for hsp in alignment.hsps:
		#			if hsp.expect < 0.1:
						print ("****Alignment****")
						print ("sequence:", alignment.title)
						print ("length:", alignment.length)
						print ("E-value:", hsp.expect)
						print (hsp.query[:50] + "...")
						print (hsp.match[:50] + "...")
						print (hsp.sbjct[:50] + "...")
result.close()	

	


doBlast ("thrombin.fa")



#BLAST running locally --> output 
#from output blast: selectsequencies with less than 1e-5 e-value. Extract outputfile: id/seq/evalue/homology/coverage/
#filter outputfile: maxium 15 sequences present in both outputfiles and save each hit  in two files for each protein
#do clustalW for both files
#construct phylogenetic tree: tree representation  
#construct matrix
#r pearson correlation 




