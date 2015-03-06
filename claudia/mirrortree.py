from test import * 
import sys
import os
import os.path
import argparse
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio.Phylo.Consensus import *
import numpy
from numpy import corrcoef, arange
import math
from pylab import pcolor, show, colorbar, xticks, yticks
import scipy.signal
import matplotlib.pyplot as plt

###################################################################################
###################################################################################
#						 HOW TO RUN THIS SCRIPT 								  #
#		put in the command line: python3 mirrortree.py -i fastafile.fa -v         #
# 																				  #
#       for the moment only starting at the beginning with one fasta infput file  #
#		To run this script you need:  										      #
#			- internet connection   											  #
#			- clustalw intalled		
#			- need NumPy 	
#			- need Treeconstruction module									      #
###################################################################################
###################################################################################

parser = argparse.ArgumentParser(description="This program manages user arguments")

parser.add_argument('-i', '--input',
            dest = "infile",
            action = "store", 
            default = False,
            help = "Input file with both proteins in FASTA format (.fa/.fasta)") 

parser.add_argument('-a1', '--align1',
            dest = "alignment1",
            action = "store", 
            default = False,
            help = "Input first CLUSTALW alignment (.aln)") 

parser.add_argument('-a2', '--align2',
            dest = "alignment2",
            action = "store", 
            default = False,
            help = "Input second CLUSTALW alignment (.aln)") 


# parser.add_argument('-o', '--output',
#                 dest = "output",
#                 action = "store",
#                 default = "./",
#                 help = "Print output in an output file")

parser.add_argument('-v', '--verbose',
            dest = "verbose",
            action = "store_true",
            default = False,
            help = "Print log in stderr")

parser.add_argument('-e', '--evalue',
            dest = "evalue",
            action = "store",
            default = 0.00001,
            type= int,
            help = "Specify e-value thereshold for selection of hits in BLAST output")	

options = parser.parse_args()

input_file = options.infile
input_aln1 = options.alignment1
input_aln2 = options.alignment2

verbose = options.verbose

#CAPTURING THE INPUT FILE 
try:
	if os.path.isfile(input_file) and input_file.endswith(".fa" or ".fasta"):
		if verbose:
			sys.stderr.write("You have selected '%s' file to execute this awesome mirrorTree script from the beginning. \n" %(input_file))
	elif os.path.isfile(input_aln1) and os.path.isfile(input_aln2) and input_aln1.endswith (".aln") and input_aln2.endswith(".aln"):
		sys.stderr.write ("You have selected '%s' and '%s' file. Starting from alignment directly. \n" %(input_aln1, input_aln2))
	elif os.path.isdir(input_file):
		sys.exit()
	else:
		if verbose:
			sys.stderr.write("Not a correct file to run this file. You may want to see documentation.\n")
		sys.exit()
except:
	sys.stderr.write("Sorry.\n")
	sys.exit() 

#EXECUTE BLAST WITH INPUT FILE
if input_file:
	if verbose:
		sys.stderr.write("Connecting to BLAST... \n\n")

		blastlist = doBlast(input_file)

	if verbose:
		sys.stderr.write("BLAST has finished.\n\n")


	#EXTRACT BEST HITS FROM BLAST XML FILE 
	evalue = options.evalue

	if verbose:
		sys.stderr.write("Extracting sequencies with e-value: '%s' from %s file after BLAST...\n" %(evalue, blastlist[0]))

	protfile1 = selectProt(blastlist[0], evalue)

	if verbose:
		sys.stderr.write("Done!\n")
		sys.stderr.write("Extracting sequencies with e-value: '%s' from %s file after BLAST...\n" %(evalue, blastlist[1]))

	protfile2 = selectProt (blastlist[1], evalue)

	if verbose:
		sys.stderr.write("Done!\n\n")
		sys.stderr.write("Comparing both BLAST output files to extract hit sequences (BEWARE! They are the same species for both proteins)\n")

	#COMPARE BOTH FILES TO ONLY SELECT HITS PRESENT IN BOTH FILES (IN PRESENTS OF PARALOGS IGNORE THEM)

	multifastafiles = comparefiles(protfile1, protfile2)

	# multifastafiles = comparefiles("1COW.out.blast","3D49.out.blast")


	#DO CLUSTAL ALIGNMENT

	# if verbose:
	# 	sys.stderr.write("Done!\n\n")
	# 	sys.stderr.write("Doing ClustalW alignment from multifasta for both proteins \n") #We should start program from here also

	path_clustal = "/Volumes/clustalw-2.1-macosx/clustalw-2.1-macosx/clustalw2" 

	cline1 = ClustalwCommandline(path_clustal, infile=multifastafiles[0])
	cline2 = ClustalwCommandline(path_clustal, infile=multifastafiles[1])
	cline1()
	cline2()


#DO FILOGENETIC TREE

# # if verbose:
# # 	sys.stderr.write("ClustalW done!\n")
# # 	sys.stderr.write("Drawing phylogenetic trees...\n")
# tree1 = Phylo.read("%s.dnd"%(multifastafiles[0][:-3]), "newick")
# tree2 = Phylo.read("%s.dnd"%(multifastafiles[1][:-3]), "newick")

# tree1.ladderize()
# tree2.ladderize()

# Phylo.draw_ascii(tree1) #IT WOULD BE INTERESTING TO HAVE GOOD NAMES IN EACH BRANCH.
# Phylo.draw_ascii(tree2) #Estos no los queremos. Queremos los de despues de tener matrix

# # if verbose:
# # 	sys.stderr.write("Phylogenetic tree done!\n")
# # 	sys.stderr.write("Saving the tree in pdf format")

# # Phylo.draw(tree1) # NEED TO OBTAIN A CORRECT FORMAT TO OUTPUT
# # Phylo.draw(tree2)

# #we must construct a tree NJ
					 
#DISTANCE MATRIX
if input_file:
	aln1 = AlignIO.read("%s.aln" %(multifastafiles[0][:-3]), 'clustal')
	aln2 = AlignIO.read("%s.aln" %(multifastafiles[1][:-3]), 'clustal')

elif input_aln1 and input_aln2:
	aln1 = AlignIO.read(input_aln1, 'clustal')
	aln2 = AlignIO.read(input_aln2, 'clustal')

calculator = DistanceCalculator("blosum62") #using identity, you can use also blosum62/identity

print (aln1)
dm1 = calculator.get_distance(aln1)
print (dm1)


print (aln2)
dm2 = calculator.get_distance(aln2)
print (dm2)


constructor = DistanceTreeConstructor(calculator, 'nj')
tree3 = constructor.build_tree(aln1)
# print(tree3)

tree4 = constructor.build_tree(aln2)
# print(tree4)


#BOOSTRAP

# msa1 = AlignIO.read("%s.aln" %(multifastafiles[0][:-3]), 'clustal')
# msas1 = bootstrap(msa1, 100)

# msa2 = AlignIO.read("%s.aln" %(multifastafiles[1][:-3]), 'clustal')
# msas2 = bootstrap(msa2, 100)

# constructor_boot = DistanceTreeConstructor(calculator, 'nj')
# # trees = bootstrap_trees(msa1, 100, constructor_boot)

# consensus_tree = bootstrap_consensus(msa1, 100, constructor_boot, majority_consensus)
# print(tree3)
# Phylo.draw_ascii(tree3)
# print(consensus_tree)

# Phylo.draw_ascii(consensus_tree)


print()

def read_matrix(matrix):
	values = []
	for element in matrix:
		for i in element:
			values.append(i)
	average = sum(values)/len(values)
	return (values,average)

#print(read_matrix(dm1))
# def formula(value1,average1,value2,average2):


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

def compute_r(matrix1,matrix2):
	(numerator,r_square,s_square) = (0,0,0)
	difference = list(zip(diff(matrix1),diff(matrix2)))
	for element in difference:
 		numerator += element[0]*element[1]
 		r_square += element[0]**2
 		s_square += element[1]**2
	# return numerator/(math.sqrt(r_square*s_square))

	return numpy.corrcoef(listmatrix(matrix1), listmatrix(matrix2))[0, 1]

print(compute_r(dm1,dm2))




def plotData(matrix1,matrix2):
	llista = list(zip(listmatrix(matrix1),listmatrix(matrix2)))
	plt.scatter(listmatrix(matrix1), listmatrix(matrix2))
	fit =numpy.polyfit(listmatrix(matrix1), listmatrix(matrix2),1)
	p = numpy.poly1d(fit)
	plt.plot(listmatrix(matrix1), p(listmatrix(matrix1)), '--g')
	plt.show() #no .show sino .save para guardarlo



plotData(dm1, dm2)





