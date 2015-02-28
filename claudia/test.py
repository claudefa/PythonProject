from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Seq
import sys
import re
# from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo

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

def Protein_creator(fin):
    fd=open(fin, "r")
    i=0
    for line in fd:
        if i==0:
            i+=1
        elif i==1:
            sp=line[14:-4]
            i+=1
        elif i==2:
            seq_id=line.strip()[8:]
            i+=1
        elif i==3:
            i+=1
        elif i==4:
            e_val=line.strip()[9:]
            i+=1
        elif i==5:
            score=line.strip()[7:]
            i+=1
        elif i==6:
            identity=line.strip()[10:]
            i+=1
        elif i==7:
            i+=1
        elif i==8:
            seq=line.strip()
            yield Protein(seq_id, seq, sp, e_val, score, identity)
            i=0
    fd.close()


def comparefiles (file1, file2):
	"""
	Function to compare both blast out and select homologous protein in the same species for both proteins. Return set.
	"""
	set1 = set()
	set2 = set()
	prot1 = []
	prot2 = []
	for protein in Protein_creator(file1):
		prot1.append(protein)
	
	for protein in Protein_creator(file2):
		prot2.append(protein)

	for protein in prot1:
		set1.add(protein.get_specie())

	for protein in prot2:
		set2.add(protein.get_specie())

	intersect = set()
	intersect = set1.intersection(set2) #not handleing finding two proteins of the same species. but when traversing 
	
	  									#in file, since the better hits are at the beginning you must keep them. 
	out1 = open("fasta1.fa","w")
	out2 = open ("fasta2.fa", "w")

	sp_list=[]
      

	for protein in Protein_creator(file1):
		for element in intersect:
			if protein.sp == element:
				if protein.sp not in sp_list:
					out1.write(">"+str(protein.get_id())+"\n"+str(protein.get_seq()).replace("-","")+"\n")
					sp_list.append(protein.sp)

	sp_list=[]
	
	for protein in Protein_creator(file2):
		for element in intersect:
			if protein.sp == element:
				if protein.sp not in sp_list:
					out2.write(">"+str(protein.get_id())+"\n"+str(protein.get_seq()).replace("-","")+"\n")
					sp_list.append(protein.sp)
		
comparefiles("1COW.out.blast","3D49.out.blast")


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
# selectProt("1COW.xml")
# selectProt("3D49.xml")



#BLAST running locally --> output 
#from output blast: selectsequencies with less than 1e-5 e-value. Extract outputfile: id/seq/evalue/homology/coverage/
#filter outputfile: maxium 15 sequences present in both outputfiles and save each hit  in two files for each protein
#do clustalW for both files
#construct phylogenetic tree: tree representation  
#construct matrix
#r pearson correlation 





