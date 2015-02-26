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
	sys.stderr.write("Extracting sequencies with e-value: %s \n" %(evalue))
	p = re.compile( "[^\[]*\[([^\]]*)\]") #this is the pattern to 
	for n in NCBIXML.parse(result):
		for alignment in n.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < evalue:
					out.write("#"*10 + "Alignment" + "#"*10 +"\n")
					m = p.findall(alignment.hit_def)
					out.write ("Hit specie: %s \n" %m)
					out.write("Hit id: %s \n" %alignment.hit_id)
					# out.write("Def: %s \n" %(alignment.title[0:50]))
					out.write("Length: %s \n" %(alignment.length) )
					out.write("E-value: %s \n" %(hsp.expect))
					out.write("Score: %s \n" %(hsp.score))
					# coverage=(int(hsp.align_length)/int(hsp.query_length)*100)
					# out.write("Coverage: \t %s \n" %(coverage)) #elbaaaa help no trobo el coverage
					out.write("Identity: %s \n" %hsp.identities)
 					# out.write(str(hsp.query)+"\n")
 					# out.write(str(hsp.match) + "...\n")
					out.write("Hit Sequence:\n%s \n" %(hsp.sbjct) )
	
	result.close()



def comparefiles (file1, file2):
	"""
	Function to compare both blast out and select homologous protein in the same species for both proteins. 
	"""
	blast1 = open(file1, "r")
	blast2 = open(file1, "r")
	out1 = open("out1.txt", "w")
	out2 = open("out2.txt", "w")
	set1 = set()
	set2 = set()

	for line in blast1:
		line = line.strip()
		if line.startswith("Hit specie"):
			set1.add(line[11:])

	for line in blast2:
		line = line.strip()
		if line.startswith("Hit specie"):
			set2.add(line[11:])
	
	intersect = set()
	# print (set1, "\n")
	# print (set2, "\n")

	intersect = set1.intersection(set2)
	# print (intersect)

	blast1.close()

	blast1 = open(file1, "r")

	##########ELBAAAAAAAAA###########
	#como puedo, una vez encontrada la especie (esta en el intersect set) cojer las filas de hit id y sequence????
	

	for line in blast1:
		line = line.strip()
		for element in intersect:
			if element in line:
				pass

	blast1.close()
	blast2.close()
	out1.close()
	out2.close()

#doBlast ("fasta.fa")
# selectProt("1COW.xml")
# selectProt("3D49.xml")
comparefiles("1COW.out.blast","3D49.out.blast")


#BLAST running locally --> output 
#from output blast: selectsequencies with less than 1e-5 e-value. Extract outputfile: id/seq/evalue/homology/coverage/
#filter outputfile: maxium 15 sequences present in both outputfiles and save each hit  in two files for each protein
#do clustalW for both files
#construct phylogenetic tree: tree representation  
#construct matrix
#r pearson correlation 





