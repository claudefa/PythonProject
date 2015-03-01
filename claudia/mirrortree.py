from test import * 
import sys
import os
import os.path
import argparse
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
# from TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
import numpy

###################################################################################
###################################################################################
#						 HOW TO RUN THIS SCRIPT 								  #
#		put in the command line: python3 mirrortree.py -i fastafile.fa -v         #
# 																				  #
#       for the moment only starting at the beginning with one fasta infput file  #
#		To run this script you need:  										      #
#			- internet connection   											  #
#			- clustalw intalled		
#			- need NumPy 													      #
###################################################################################
###################################################################################

# parser = argparse.ArgumentParser(description="This program manages user arguments")

# parser.add_argument('-i', '--input',
#             dest = "infile",
#             required = True,
#             action = "store", 
#             default = False,
#             help = "Input file with both proteins in FASTA format") #canviar per tal de poder entrar des d'altres punts del programa

# # parser.add_argument('-o', '--output',
# #                 dest = "output",
# #                 action = "store",
# #                 default = "./",
# #                 help = "Print output in an output file")

# parser.add_argument('-v', '--verbose',
#             dest = "verbose",
#             action = "store_true",
#             default = False,
#             help = "Print log in stderr")

# parser.add_argument('-e', '--evalue',
#             dest = "evalue",
#             action = "store",
#             default = 0.00001,
#             type= int,
#             help = "Specify e-value thereshold for selection of hits in BLAST output")	

# options = parser.parse_args()

# input_file = options.infile
# verbose = options.verbose

# #CAPTURING THE INPUT FILE 
# try:
# 	if os.path.isfile(input_file) and input_file.endswith(".fa" or ".fasta"):
# 		if verbose:
# 			sys.stderr.write("You have selected '%s' file to execute this awesome mirrorTree script from the beginning. \n" %(input_file))
# 	elif os.path.isfile(input_file) and input_file.endswith (".aln"):
# 		sys.stderr.write ("You have selected '%s' file. Starting from alignment directly. \n" %(input_file))
# 	elif os.path.isdir(input_file):
# 		sys.exit()
# 	else:
# 		if verbose:
# 			sys.stderr.write("Not a correct file to run this file. You may want to see documentation.\n")
# 		sys.exit()
# except:
# 	sys.stderr.write("Sorry.\n")
# 	sys.exit() 

# #EXECUTE BLAST WITH INPUT FILE
# if verbose:
# 	sys.stderr.write("Connecting to BLAST... \n\n")

# blastlist = doBlast(input_file)

# if verbose:
# 	sys.stderr.write("BLAST has finished.\n\n")


# #EXTRACT BEST HITS FROM BLAST XML FILE 
# evalue = options.evalue

# if verbose:
# 	sys.stderr.write("Extracting sequencies with e-value: '%s' from %s file after BLAST...\n" %(evalue, blastlist[0]))

# protfile1 = selectProt(blastlist[0], evalue)

# if verbose:
# 	sys.stderr.write("Done!\n")
# 	sys.stderr.write("Extracting sequencies with e-value: '%s' from %s file after BLAST...\n" %(evalue, blastlist[1]))

# protfile2 = selectProt (blastlist[1], evalue)

# if verbose:
# 	sys.stderr.write("Done!\n\n")
# 	sys.stderr.write("Comparing both BLAST output files to extract hit sequences (BEWARE! They are the same species for both proteins)\n")

# #COMPARE BOTH FILES TO ONLY SELECT HITS PRESENT IN BOTH FILES (IN PRESENTS OF PARALOGS IGNORE THEM)

# multifastafiles = comparefiles(protfile1, protfile2)

multifastafiles = comparefiles("1COW.out.blast","3D49.out.blast")


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

# if verbose:
# 	sys.stderr.write("ClustalW done!\n")
# 	sys.stderr.write("Drawing phylogenetic trees...\n")
tree1 = Phylo.read("%s.dnd"%(multifastafiles[0][:-3]), "newick")
tree2 = Phylo.read("%s.dnd"%(multifastafiles[1][:-3]), "newick")

tree1.ladderize()
tree2.ladderize()

Phylo.draw_ascii(tree1) #IT WOULD BE INTERESTING TO HAVE GOOD NAMES IN EACH BRANCH.
Phylo.draw_ascii(tree2)

# if verbose:
# 	sys.stderr.write("Phylogenetic tree done!\n")
# 	sys.stderr.write("Saving the tree in pdf format")

# Phylo.draw(tree1) # NEED TO OBTAIN A CORRECT FORMAT TO OUTPUT
# Phylo.draw(tree2)

#we must construct a tree NJ
 
#DISTANCE MATRIX
aln1 = AlignIO.read("msa.phy", "phylip") #no se m'obra be
calculator = DistanceCalculator("identity") #using identity, you can use also blosum62

print (aln1)
dm = calculator.get_distance(aln)
print (dm)





# constructor = DistanceTreeConstructor(calculator, 'nj')

# tree = constructor.build_tree(aln)




# #MATRIX DISTANCES - NOT PROPERLY DONE (TRYING THE ONE ABOVE)
# # if verbose:
# # 	sys.stderr ("Computing distances between clades...\n")

# matrix1 = to_distance_matrix(tree1)

# matrix2 = to_distance_matrix(tree2)

# for element in matrix1:
# 	print (element)

# for element in matrix2:
# 	print (element)



