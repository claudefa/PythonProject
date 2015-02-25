from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Seq
import sys
import re

def doBlast (fastafile):
	
	handle = open(fastafile)
	# count = 0
	
	for record in SeqIO.parse(handle, "fasta") :
		print ("Doing blast %s" %record.id)
		try:
			result = NCBIWWW.qblast("blastp", "swissprot", record.seq)
		except:
			print ("Impossible to do blast")
		print("Blast %s done..." %record.id)
		
		blastfile = open("%s.xml" %(record.id[:4]), "w")
		blastfile.write(result.read())
		blastfile.close()
		# count += 1
	
	# if count ==2:
	# 	print("correct file")
	# else:
	# 	print ("please introduce a fasta with two sequences")
	# 	sys.exit()

	# blastfile.close()	
	# blastfile2.close()	

	handle.close()

# def doBlast (fastafile):
# 	handle = open(fastafile)
# 	for record in SeqIO.parse(handle, "fasta") :
# 		result = NCBIWWW.qblast("blastp", "nr", record.seq)
# 		blastfile = open("blast.xml", "w")
# 		blastfile.write(result.read())
# 		blastfile.close()	
# 	handle.close()


def selectProt(blastxml):
	result = open(blastxml, "r")
	out = open("candidates.blast", "w")
	evalue = 0.00001
	p = re.compile( "Precursor[^\[]*\[([^\]]*)\]")
	for n in NCBIXML.parse(result):
		for alignment in n.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < evalue:
					out.write("****Alignment****\n")
					out.write("Hit id: %s \n" %alignment.hit_id)
					m = p.findall(alignment.hit_def)
					out.write ("Hit specie: %s \n" %m)
					out.write("Def: %s \n" %(alignment.title[0:50]))
					out.write("Length: %s \n" %(alignment.length) )
					out.write("E-value: %s \n" %(hsp.expect))
					out.write("Score: %s \n" %(hsp.score))
# 					out.write(str(hsp.query)+"\n")
# 					out.write(str(hsp.match) + "...\n")
					out.write("Hit Sequence: %s \n" %(hsp.sbjct) )
	
	result.close()
doBlast ("thrombin.fa")
selectProt("blast.xml")


#BLAST running locally --> output 
#from output blast: selectsequencies with less than 1e-5 e-value. Extract outputfile: id/seq/evalue/homology/coverage/
#filter outputfile: maxium 15 sequences present in both outputfiles and save each hit  in two files for each protein
#do clustalW for both files
#construct phylogenetic tree: tree representation  
#construct matrix
#r pearson correlation 





