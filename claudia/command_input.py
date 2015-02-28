import test 
import sys
import os
import os.path
import sequence_data
import argparse


parser = argparse.ArgumentParser(description="This program manages user arguments")

parser.add_argument('-i', '--input',
            dest = "infile",
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

#CAPTURING THE INPUT FILE(s)
while True:
	try os.path.isfile(input_file):
		if verbose:
			sys.stderr.write("You have selected '%s' file to execute this awsome mirrorTree script" %(input_file))
		
	except FileNotFoundError:
        sys.stderr.write("This file does not exist. Check your spelling or choose another file")
        sys.exit() 




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

