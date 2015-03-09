
from test import *

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
evalue = options.evalue

#CAPTURING THE INPUT FILE(S)

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

#EXECUTE BLAST WITH INPUT FILE

#IF INPUT IS ONE FASTA FILE

if len(input_list) == 1:
	if verbose:
		sys.stderr.write("Connecting to BLAST... \n\n")

	#blastlist = doBlast(input_list[0])
	blastlist = ["1AIE.xml","2J0I.xml"]

	if verbose:
		sys.stderr.write("BLAST has finished.\n\n")


	#EXTRACT BEST HITS FROM BLAST XML FILE 

	(protfile_list, aln) = ([],[])

	for xml in blastlist:
		if verbose:
			sys.stderr.write("Extracting sequencies with e-value: '%s' from %s file after BLAST...\n" %(evalue, xml))

		protfile_list.append(selectProt(xml, evalue))

		if verbose:
			sys.stderr.write("Done!\n")
	if verbose:
		sys.stderr.write("Comparing both BLAST output files to extract hit sequences (BEWARE! They are the same species for both proteins)...\nPerforming ClustalW alignment...\n")

	for element in comparefiles(protfile_list):
		doClustalW(element) #PERFORMING CLUSTALW ALIGNMENT
		aln.append(AlignIO.read("%s.aln" %(element[:-3]), 'clustal')) #BUILDING DISTANCE MATRIX

	if verbose:
		sys.stderr.write("Comparison and ClustalW for both files done!\n")
		sys.stderr.write("Obtaining distance matrices from the alignments...\n")


# IF INPUT ARE TWO .aln FILES

else:
	aln = []
	for element in input_list:
		aln.append(AlignIO.read(element, 'clustal')) #BUILDING DISTANCE MATRIX
	if verbose:
		sys.stderr.write("Obtaining distance matrices from the alignments...\n")

#BUILDING FILOGENETIC TREE

calculator = DistanceCalculator("blosum62") # You can use blosum62/identity
constructor = DistanceTreeConstructor(calculator, 'nj') #Neighbour Joining = 'nj'

(dm,tree) = ([],[])
for element in aln:
	dm.append(calculator.get_distance(element))
#dm2 = calculator.get_distance(aln2)
	tree.append(constructor.build_tree(element))
#tree4 = constructor.build_tree(aln2)

if verbose:
	sys.stderr.write("Phylogenetic tree done!\n")

for element in tree:
	Phylo.draw_ascii(element)
	print()

#Phylo.draw_ascii(tree4)
#print()

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

#COMPUTE R CORRELATION 
if verbose:
	sys.stderr.write("Phylogenetic tree done!\n")

sys.stdout.write("This is the correlation for both proteins: %.3f \n"%(compute_r(dm)))

if verbose:
	sys.stderr.write("Plotting linear regression. Saved as 'plot.png'\n")

#plotData(dm1, dm2)


