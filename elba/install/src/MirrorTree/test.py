#!/usr/bin/python3

from modules import *

def doBlast (fastafile):
	"""
	This function performs a BLAST search using the QBLAST server at NCBI.
	Returns a xml file for each protein with blast output.
	"""
	try:
		handle = open(fastafile, "r")
	except IOError:
		sys.stderr.write("Impossible to open this file. It does not exist!\nAborting...\n")
		sys.exit()
	
	blastlist = []
	#Do blast for each fasta sequence in file. Not handeling two input files, only one!
	
	for record in SeqIO.parse(handle, "fasta"):
		sys.stderr.write("Doing blast %s ...\n" %(record.id[:4])) #the name of each fasta must be revised (now we suposo that it is a pdb but what if not??)
		try:
			result = NCBIWWW.qblast("blastp", "swissprot", record.seq)
			# pass # ONLY TO UNCOMMEND WHEN YOU DON'T WANT TO RUN BLAST BUT CHECK THE SCRIPT FLOW. OTHERWISE BLAST KICKS YOU OUT 
		except:
			sys.stderr.write("Impossible to perform BLAST!\nAborting...\n") # check error type
			sys.exit()
		
		sys.stderr.write("Blast %s done!\n" %(record.id[:4]))

		blastfile = open("%s.xml" %(record.id[:4]), "w")
		blastfile.write(result.read()) 
		# blastfile.write(str(result))# ONLY TO UNCOMMEND WHEN YOU DON'T WANT TO RUN BLAST BUT CHECK THE SCRIPT FLOW. OTHERWISE BLAST KICKS YOU OUT 
		blastfile.close()

		blastlist.append(str(record.id[:4])+".xml")
		sys.stderr.write("Blast output with extension '%s.xml'\n\n" %(record.id[:4]))

	return blastlist 
	#return ["1AIE.xml","2J0I.xml"] # ONLY TO UNCOMMEND WHEN YOU DON'T WANT TO RUN BLAST BUT CHECK THE SCRIPT FLOW. OTHERWISE BLAST KICKS YOU OUT 
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
	p = re.compile("\[(.*)\]") # Pattern for extract the specie
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
	out = open(outfile, "w")

	for protein in Protein_creator(filename):
		for specie in intersect:
			if protein.get_specie() == specie:
				if protein.get_specie() not in sp_set:
					out.write(">"+str(protein.get_id())+"\n"+str(protein.get_seq()).replace("-","")+"\n")
					sp_set.add(protein.get_specie())
	out.close()
	return()


def comparefiles (file_list):
	"""
	Function to compare both blast out and select homologous protein in the same species for both proteins. Return set.
	"""
	(set1,set2) = (set(),set())

	for protein in Protein_creator(file_list[0]):
		set1.add(protein.get_specie())

	for protein in Protein_creator(file_list[1]):
		set2.add(protein.get_specie())

	intersect = set1.intersection(set2) #Get species shared in both files

	species_selector(intersect, file_list[0], "multifasta1.fa") #CANVIAR NOM PER TAL Q SIGUI VARIABLE
	species_selector(intersect, file_list[1], "multifasta2.fa")
		
	return ["multifasta1.fa","multifasta2.fa"]


def doClustalW (multifastafile):
	"""
	Given a multifasta file peform an alignment using ClustalW. Return two files: .aln and .dnd
	"""
	path_clustal = "/usr/bin/clustalw" 
	cline1 = ClustalwCommandline(path_clustal, infile=multifastafile)
	cline1()


def read_matrix(matrix):
	"""
	Traversing throught matrix and return values and average
	"""
	values = []
	for element in matrix:
		for i in element:
			values.append(i)
	average = sum(values)/len(values)
	return (values,average)

def diff(matrix):
	diff = []
	average = read_matrix(matrix)[1]
	for element in read_matrix(matrix)[0]:
		diff.append(element-average)
	return diff

def listmatrix (matrix): #change, it is the same as readmatrix but without average
	values = []
	for element in matrix:
		for i in element:
			values.append(i)
	return values

def compute_r(matrix_list):
	# (numerator,r_square,s_square) = (0,0,0)
	# difference = list(zip(diff(matrix1),diff(matrix2)))
	# for element in difference:
 # 		numerator += element[0]*element[1]
 # 		r_square += element[0]**2
 # 		s_square += element[1]**2
	# return numerator/(math.sqrt(r_square*s_square))

	return numpy.corrcoef(listmatrix(matrix_list[0]), listmatrix(matrix_list[1]))[0, 1]


def plotData(matrix1,matrix2):
	x = listmatrix(matrix1)
	y = listmatrix(matrix2)
	plt.scatter(x,y)
	fit =numpy.polyfit(x,y,1)
	p = numpy.poly1d(fit)
	plt.plot(x, p(x), '--g')
	title('Linear regression')
	plt.xlabel('Matrix1') #canviar noms per ser variables
	plt.ylabel('Matrix2')
	savefig("plot.png") 
