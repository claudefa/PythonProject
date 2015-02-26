from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Seq
import sys
import re

def doBlast (fastafile):
	"""
	Function to perform blast online. Internet connection needed. Returns a xml file for each protein with blast output.
	"""
	handle = open(fastafile)

	#Do blast for each fasta sequence in file. Not handeling right now two input files, only one!
	
	for record in SeqIO.parse(handle, "fasta"):
		sys.stderr.write("Doing blast %s ..." %(record.id[:4]))
		try:
			result = NCBIWWW.qblast("blastp", "swissprot", record.seq)
			# result = record.seq
		except:
			sys.stderr.write("Impossible to do blast, check your internet connection") # check error type
			sys.exit()
		sys.stderr.write("Blast %s done!" %(record.id[:4]))
		
		blastfile = open("%s.xml" %(record.id[:4]), "w")
		blastfile.write(result.read())
		# blastfile.write(str(result))
		blastfile.close()

	handle.close()



def selectProt(blastxml):
	"""
	Function to extract sequencies with selected e-value form xml blast output. 
	"""
	result = open(blastxml, "r")
	out = open("%s.out.blast" %(blastxml[:4]), "w")
	evalue = 0.00001
	coverage=0.0
	sys.stderr.write("Extracting sequencies with e-value: %s " %(evalue))
	p = re.compile( "Precursor[^\[]*\[([^\]]*)\]") #this is the pattern to 
	for n in NCBIXML.parse(result):
		for alignment in n.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < evalue:
					out.write("*********Alignment**********\n")
					out.write("Hit id: \t %s \n" %alignment.hit_id)
					m = p.findall(alignment.hit_def)
					out.write ("Hit specie: \t %s \n" %m)
					# out.write("Def: %s \n" %(alignment.title[0:50]))
					out.write("Length: \t %s \n" %(alignment.length) )
					out.write("E-value: \t %s \n" %(hsp.expect))
					out.write("Score: \t %s \n" %(hsp.score))
					# coverage=(int(hsp.align_length)/int(hsp.query_length)*100)
					# out.write("Coverage: \t %s \n" %(coverage)) #elbaaaa help no trobo el coverage
					out.write("Identity: \t %s \n" %hsp.identities)
 					# out.write(str(hsp.query)+"\n")
 					# out.write(str(hsp.match) + "...\n")
					out.write("Hit Sequence: \t %s \n" %(hsp.sbjct) )
	
	result.close()

#doBlast ("fasta.fa")
selectProt("1COW.xml")
selectsequencies("3D49.xml")


#BLAST running locally --> output 
#from output blast: selectsequencies with less than 1e-5 e-value. Extract outputfile: id/seq/evalue/homology/coverage/
#filter outputfile: maxium 15 sequences present in both outputfiles and save each hit  in two files for each protein
#do clustalW for both files
#construct phylogenetic tree: tree representation  
#construct matrix
#r pearson correlation 





