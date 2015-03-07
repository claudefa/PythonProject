import argparse, os, sys

parser = argparse.ArgumentParser(description="This program manages user arguments")

parser.add_argument('-i', '--input',
	dest = "infile",
	action = "store", 
	default = None,
#	type = argparse.FileType('r'), # Open each argument as a file for reading
	nargs = '+', # One or more files as input
	help = "Input: Choose one FASTA file (.fa/.fasta) or two CLUSTALW alignment files (.aln)")

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
            type = int,
            help = "Specify e-value thereshold for selection of hits in BLAST output")	

options = parser.parse_args()

input_list = options.infile
verbose = options.verbose

#CAPTURING THE INPUT FILE

if len(input_list) < 3:
	for input_file in input_list:

		if os.path.isdir(input_file):
			sys.stderr.write("This program cannot handle directories as input.\nPlease, use files.\nAborting...\n")
			sys.exit()

		if os.path.isfile(input_file) and (".fa" or ".fasta") in input_file:
			if len(input_list) is 1:
				if verbose:
					sys.stderr.write("You have selected '%s' file to execute this awesome MirrorTree script.\n" %(input_file))
			else:
				sys.stderr.write("Impossible to handle this input request. Please, check the documentation.\nAborting...\n")
				sys.exit()

		elif len(input_list) is 2:
			if os.path.isfile(input_list[0]) and (".aln") in input_list[0] and os.path.isfile(input_list[1]) and (".aln") in input_list[1]:
				if verbose:
					sys.stderr.write ("You have selected '%s' and '%s' file. Starting from alignment directly. \n"
					%(input_list[0], input_list[1],))
				break
			else:
				sys.stderr.write("Impossible to handle this input request. Please, check the documentation.\nAborting...\n")
				sys.exit()
		else:
			sys.stderr.write("Impossible to handle this input request. Please, check the documentation.\nAborting...\n")
			sys.exit()
else:
	sys.stderr.write("Only able to handle one or two files!\nAborting...\n")
	sys.exit()
