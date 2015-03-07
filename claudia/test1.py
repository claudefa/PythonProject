from modules import *

def doBlast (fastafile):
	"""
	This function performs a BLAST search using the QBLAST server at NCBI.
	Returns a xml file for each protein with blast output.
	"""
	try:
		handle = open(fastafile, "r")
	except IOError:
		sys.stderr.write("Impossible to open this file. It does not exist!\nAborting program.\n")
		sys.exit()
	
	blastlist = []
	#Do blast for each fasta sequence in file. Not handeling two input files, only one!
	
	for record in SeqIO.parse(handle, "fasta"):
		sys.stderr.write("Doing blast %s ...\n" %(record.id[:4])) #the name of each fasta must be revised (now we suposo that it is a pdb but what if not??)
		try:
			result = NCBIWWW.qblast("blastp", "swissprot", record.seq)
			# pass # ONLY TO UNCOMMEND WHEN YOU DON'T WANT TO RUN BLAST BUT CHECK THE SCRIPT FLOW. OTHERWISE BLAST KICKS YOU OUT 
		except:
			sys.stderr.write("Impossible to do blast, check your internet connection\n") # check error type
			sys.exit()
		
		sys.stderr.write("Blast %s done!\n" %(record.id[:4]))

		blastfile = open("%s.xml" %(record.id[:4]), "w")
		blastfile.write(result.read()) 
		# blastfile.write(str(result))# ONLY TO UNCOMMEND WHEN YOU DON'T WANT TO RUN BLAST BUT CHECK THE SCRIPT FLOW. OTHERWISE BLAST KICKS YOU OUT 
		blastfile.close()

		blastlist.append(str(record.id[:4])+".xml")
		sys.stderr.write("Blast output with extension '%s.xml'\n\n" %(record.id[:4]))

	return blastlist 
	# return ["1COW.xml","3D49.xml"] # ONLY TO UNCOMMEND WHEN YOU DON'T WANT TO RUN BLAST BUT CHECK THE SCRIPT FLOW. OTHERWISE BLAST KICKS YOU OUT 
	handle.close()



def selectProt(blastxml, evalue):
	"""
	Function to extract sequencies with selected e-value form xml blast output. 
	"""
	result = open(blastxml, "r")
	out = open("%s.out.blast" %(blastxml[:4]), "w")
	# coverage=0.0 IS IT IMPORTANT THE COVERAGE


	p = re.compile( "[^\[]*\[([^\]]*)\]") #this is the pattern to select the specie, unbelievable 
	for n in NCBIXML.parse(result):
		for alignment in n.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < evalue:
					out.write("#"*10 + "Alignment" + "#"*10 +"\n")
					m = p.findall(alignment.hit_def)
					out.write ("Hit_specie: %s \n" %m)
					out.write("Hit_id: %s \n" %alignment.hit_id)
					# out.write("Def: %s \n" %(alignment.title[0:50]))
					out.write("Length: %s \n" %(alignment.length) )
					out.write("E-value: %s \n" %(hsp.expect))
					out.write("Score: %s \n" %(hsp.score))
					# coverage=(int(hsp.align_length)/int(hsp.query_length)*100)
					# out.write("Coverage: \t %s \n" %(coverage)) #elbaaaa help no trobo el coverage
					out.write("Identity: %s \n" %hsp.identities)
 					# out.write(str(hsp.query)+"\n")
 					# out.write(str(hsp.match) + "...\n")
					out.write("Hit_Sequence: %s \n" %(hsp.sbjct) )
	
	name = str(blastxml[:4]) + ".out.blast"

	result.close()
	return name



class Protein(object):

	def __init__(self, seq_id, seq, sp, e_val, score, identity):
		self.seq_id=seq_id
		self.seq=seq
		self.sp=sp
		self.e_val=e_val
		self.score=score
		self.identity=identity
    
	def __str__(self):
		return ("%s\n%s\n%s\n%s\n%s\n%s\n" %(self.seq_id, self.seq, self.sp, self.e_val, self.score, self.identity))
	def get_specie (self):
		return str(self.sp)
	def get_id (self):
		return self.seq_id
	def get_seq (self):
		return self.seq

def Protein_creator(filename):
	fd = open(filename, "r")
	info = []
	p = re.compile("\[(.*)\]")
	for field in fd:
		if "#" in field:
			if not info == []:
				yield Protein(info[1], info[6], info[0], info[3], info[4], info[5])
				info = []
		else:
			if "[" in field:
				sp = p.search(field).group(1)
				info.append(sp)
			else:
				line = field.split()
				info.append(line[1])			
	yield Protein(info[1], info[6], info[0], info[3], info[4], info[5])
	fd.close()

def species_selector(intersect, filename, outfile):
	"""
	This function selects only the higher scoring protein of each specie.
	Write a fasta file with selected proteins.
	"""
	sp_set = set()	
	out = open (outfile, "w")

	for protein in Protein_creator(filename):
		for specie in intersect:
			if protein.get_specie() == specie:
				if protein.get_specie() not in sp_set:
					out.write(">"+str(protein.get_id())+"\n"+str(protein.get_seq()).replace("-","")+"\n")
					sp_set.add(protein.get_specie())
	out.close()
	return()


def comparefiles (file1, file2):
	"""
	Function to compare both blast out and select homologous protein in the same species for both proteins. Return set.
	"""
	(set1,set2) = (set(),set())

	for protein in Protein_creator(file1):
		set1.add(protein.get_specie())

	for protein in Protein_creator(file2):
		set2.add(protein.get_specie())

	intersect = set1.intersection(set2) #Get species shared in both files

	species_selector(intersect, file1, "multifasta1.fa")
	species_selector(intersect, file2, "multifasta2.fa")
		
	return()





# cline1 = ClustalwCommandline("clustalw", infile="fasta1.fa")
# cline2 = ClustalwCommandline("clustalw", infile="fasta2.fa")
# cline1()
# cline2()
# tree1 = Phylo.read("fasta1.dnd", "newick")
# tree2 = Phylo.read("fasta2.dnd", "newick")
# Phylo.draw_ascii(tree1)
# Phylo.draw_ascii(tree2)



#Un cop tenim el fasta: 1) quedar-nos amb la sequencia que tingui el e-value millor/millor identitat
 						# 2) treure els gaps 


#EXECUTE FUNCTIONS


#doBlast ("fasta.fa")
#selectProt("1COW.xml", 0.00001)
#selectProt("3D49.xml", 0.00001)
#comparefiles("1COW.out.blast","3D49.out.blast")


#BLAST running locally --> output 
#from output blast: selectsequencies with less than 1e-5 e-value. Extract outputfile: id/seq/evalue/homology/coverage/
#filter outputfile: maxium 15 sequences present in both outputfiles and save each hit  in two files for each protein
#do clustalW for both files
#construct phylogenetic tree: tree representation  
#construct matrix
#r pearson correlation 


