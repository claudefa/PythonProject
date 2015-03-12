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
	
	#Do blast for each fasta sequence in file.
	n = 1
	for record in SeqIO.parse(handle, "fasta"):
	
		try:
			result = NCBIWWW.qblast("blastp", "swissprot", record.seq)
		except:
			sys.stderr.write("Impossible to perform BLAST!\nAborting...\n") 
			sys.exit()
	
		blastfile = open("blast%s.xml" %(n), "w")
		blastfile.write(result.read()) 
		blastfile.close()
		blastlist.append("blast"+str(n)+".xml")
		n += 1
	handle.close()
	return blastlist 
	



def selectProt(blastxml, evalue, identity):
	"""
	Function to extract sequencies with selected e-value form xml blast output. 
	"""
	#Open blast file and open output file
	result = open(blastxml, "r")
	out = open("%s.out.blast" %(blastxml[:-4]), "w")

	p = re.compile( "[^\[]*\[([^\]]*)\]") #this is the pattern to select the specie, unbelievable 
	
	#for each alignment extrat most relevant fields
	for n in NCBIXML.parse(result):
		for alignment in n.alignments:
			for hsp in alignment.hsps:
				if hsp.expect < evalue: #and hsp.identities >= identity:
					out.write("#"*10 + "Alignment" + "#"*10 +"\n")
					m = p.findall(alignment.hit_def)
					out.write ("Hit_specie: %s \n" %m)
					out.write("Hit_id: %s \n" %alignment.hit_id)
					out.write("Length: %s \n" %(alignment.length) )
					out.write("E-value: %s \n" %(hsp.expect))
					out.write("Score: %s \n" %(hsp.score))
					out.write("Identity: %s \n" %hsp.identities)
					out.write("Hit_Sequence: %s \n" %(hsp.sbjct) )
	
	name = str(blastxml[:-4]) + ".out.blast"

	result.close()
	out.close()
	return name



class Protein(object):
	"""
	Protein object related to the file of selected proteins from blast output. 
	"""
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
		"""Get specie."""
		return str(self.sp)
	def get_id (self):
		"""Get identification."""
		return self.seq_id
	def get_seq (self):
		"""Get sequence."""
		return self.seq

class NEHomologous(ValueError):
	"""
	Raises error when Not Enough Homologous found in both files to perform an accurate mirror tree. 
	"""
	def __init__ (self, number):
		self.number = number
	
	def __str__(self):
		return "Not enough homologs after blast results filtering. Sorry! :(\n"

def Protein_creator(filename):
	"""
	Function to create Protein Objectes from the selected proteins after blast. 
	"""
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

def querySequence(galiza):
	percebe = []
	percebeiro = open(galiza,"r")
	for mejillon in SeqIO.parse(percebeiro, "fasta"):
		percebe.append((mejillon.id,mejillon.seq))
	return percebe


def species_selector(intersect, filename, outfile, query):
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
	out.write(">"+str(query[0])+"\n"+str(query[1])+"\n")
	out.close()
	return()


def comparefiles (file_list, filename):
	"""
	Function to compare both blast out and select homologous protein in the same species for both proteins. Return set.
	"""
	(set1,set2) = (set(),set())

	for protein in Protein_creator(file_list[0]):
		set1.add(protein.get_specie())
	for protein in Protein_creator(file_list[1]):
		set2.add(protein.get_specie())

	intersect = set1.intersection(set2) #Get species shared in both files
	
	if len(intersect) >= 3:
		query = querySequence(filename)
		species_selector(intersect, file_list[0], "multifasta1.fa",query[0])  
		species_selector(intersect, file_list[1], "multifasta2.fa",query[1])
		return ["multifasta1.fa","multifasta2.fa"]

	else:
		sys.stderr.write("Not enough homologous after blast results filtering. Sorry! :(\n")
		return()
		sys.exit()


def doClustalW (multifastafile, path_clustal):
	"""
	Given a multifasta file peform an alignment using ClustalW. Return two files: .aln and .dnd
	"""
	cline1 = ClustalwCommandline(path_clustal, infile=multifastafile)
	cline1()

def listmatrix (matrix): 
	"""
	Traversign throught matrix and return values.
	"""
	values = []
	for element in matrix:
		for i in element:
			values.append(i)
	return values

def average(matrix):
	"""
	Traversing throught matrix and return values and average.
	"""
	values = listmatrix(matrix)
	average = sum(values)/len(values)
	return (values,average)


def compute_r(matrix_list):
	"""
	Function to compute Pearson correlation coefficient.
	"""
	return numpy.corrcoef(listmatrix(matrix_list[0]), listmatrix(matrix_list[1]))[0, 1]


def plotData(matrix_list):
	"""
	Function to Plot both distance matrix and see linear regression.
	"""
	x = listmatrix(matrix_list[0])
	y = listmatrix(matrix_list[1])
	plt.scatter(x,y)
	fit =numpy.polyfit(x,y,1)
	p = numpy.poly1d(fit)
	plt.plot(x, p(x), '--g')
	title('Linear regression for distance matrices')
	plt.xlabel('Distance Family 1') 
	plt.ylabel('Distance Family 2')
	savefig("plot.png")
	return "plot.png" 


def cleaningWorkspace(folder, files, input_list):
	p = re.compile("(.*)\..*")
	m = ""
	for element in input_list:
		m += p.search(element).group(1)+"_"
	if not os.path.isdir(m[:-1]):
		os.system("mkdir MirrorTree_%s" %(m[:-1]))
	os.system("mkdir MirrorTree_%s/%s" %(m[:-1],folder))
	for element in files:
		os.system("mv %s ./MirrorTree_%s/%s" %(element, m[:-1], folder))
	return()


