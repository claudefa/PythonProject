from test import * 
import sys
import os
import os.path
import argparse
###################################################################################
###################################################################################
#						 HOW TO RUN THIS SCRIPT 								  #
#		put in the command line: python3 mirrortree.py -i fastafile.fa -v         #
# 																				  #
#       for the moment only starting at the beginning with one fasta infput file  #
#																			      #
###################################################################################
###################################################################################

parser = argparse.ArgumentParser(description="This program manages user arguments")

parser.add_argument('-i', '--input',
            dest = "infile",
            required = True,
            action = "store", 
            default = False,
            help = "Input file with both proteins in FASTA format") #canviar per tal de poder entrar des d'altres punts del programa

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
verbose = options.verbose

#CAPTURING THE INPUT FILE 
try:
	if os.path.isfile(input_file) and input_file.endswith(".fa" or ".fasta"):
		if verbose:
			sys.stderr.write("You have selected '%s' file to execute this awesome mirrorTree script from the beginning. \n" %(input_file))
	elif os.path.isfile(input_file) and input_file.endswith (".aln"):
		sys.stderr.write ("You have selected '%s' file. Starting from alignment directly. \n" %(input_file))
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
if verbose:
	sys.stderr.write("Connecting to BLAST... \n\n")

xmlblast = doBlast(input_file)

if verbose:
	sys.stderr.write ("BLAST has finished.\n\n")


# READING THE SEQUENCE FILES AND STORE THEM IN A LIST

# protein_sequence_objects_list = []

# for filename in list_files:
# 	protein_sequence_objects_list.extend( (d.translate()for d in FASTA_iterator(filename) ) )
# 	if verbose:
# 		sys.stderr.write("%s finished.\n" %filename)

# if verbose:
# 	sys.stderr.write("%s sequences found.\n" %len(protein_sequence_objects_list))


# # SELECT OUTPUT
# outputfile = options.output

# if outputfile == False :	
# 	out_fd = sys.stdout

# else:	
# 	if outputfile.endswith (".gz"):
# 		out_fd = gzip.open(outputfile, "wt") 
# 	else:		
# 		out_fd = open(outputfile, "w")

# #SELECT PROTEINS WITH A PATTERN
# pattern = options.pattern
# m=[]
# if pattern == False:
# 	if verbose:
# 			sys.stderr.write("No pattern selected\n")
# 	m = protein_sequence_objects_list
# else:
# 	p = re.compile(pattern, re.IGNORECASE) 
# 	for protein in protein_sequence_objects_list:
# 		if p.search(protein.get_sequence()) != None:
# 			m.append(protein) 
	

# #SELECT NUMBER OF SEQUENCES TO BE PRINTED
# randomnum = options.random 
# if randomnum:
# 	try:
# 		if verbose:
# 			sys.stderr.write("Selecting a sample with size %s.\n" %randomnum)
# 		num = random.sample(m, randomnum)
# 	except ValueError:
# 		sys.stderr.write ("%s is greater than the number of sequences! Try again (less than %s)\n" %(randomnum, len(m)))
# 		sys.exit()

# else:
# 	if verbose:
# 		sys.stderr.write("No sampling done. \n")
# 	num = m

# # SORT THE SEQUENCES BY THE LENGTH
# if verbose:
# 	sys.stderr.write("Sorting the sequences...\n")

# num.sort(reverse=True)

# if verbose:
# 	sys.stderr.write("Sort process finished.\n")

# # WRITE OUTPUT
# for protein in num:
# 	out_fd.write("%s\t%s\t%s\n" %(	protein.get_identifier(), len(protein), protein.get_mw()))


# out_fd.close()

# if verbose:
# 	sys.stderr.write("Program finished correctly.\n")

