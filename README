# MIRRORTREE 1.0 

This is the PYT and SBI project by [**Elba Raimúndez**](https://github.com/elbaraim), [**Clàudia Fontserè**](https://github.com/claudefa) and [**Lucas Michel**]!

- [Description](#description)
- [Requirements before installing](#requirements)
- [Installing](#installing)
- [Organization of the code](#codeorganization)
    * [mirrorTree](#mirrortree)
    * [functions.py](#functions)
    * [modules.py](#modules)
- [How to execute](#execute)

##Description
MirrorTree 1.0 is a program to predict protein-protin interaction using similarity between two protein families (Pazos and Valencia, 2001).

It can initiated from 2 protein sequences (in the same fasta file with fa/fasta extension) and from 2 multiple sequence alignments (two files with aln extension).

The workflow of the mirrorTree is the following:
    - Input a fasta file with two proteins (fasta format).
    - Find orthologs from each protein. This is done by connecting to BLAST (online) and comparing to swissprot/uniprot database. 
        * The parameters controlling the BLAST search are % identity (according with the BLAST alignment) and e-value. By default these parameters are set to ≥30% and ≤1e-5 and ≥60% respectively, and they can be modified by the user. 
        * Retain sequences that belong to the same species in both families. And only 1 sequence per organism, the one with highest homology with the query is chosen. 
    - The resulting sequences are aligned using ClustalW (Thompson JD, 1994) with default parameters. This step is performed locally. 
    As a result of this process, we obtain a multiple sequence alignment of the orthologs for our two query sequences.
    - (The user may start the execution here)
    - Afterwards, distance matrices for both alignments are computed.
    Phylogenetic trees are obtained from these alignments with the neighbor-joining (NJ) algorithm implemented in ClustalW (Chenna, et al., 2003) using bootstrap (100 repetitions).
    - Finally, The tree similarity between the two families is calculated as the correlation between their distance matrices using Pearson correaltion (Pazos and Valencia, 2001)


##Requirements before installing
- Beware that this script is written in Python 3.4.2. Check your version before executing mirrorTree.
- Required python modules
    * **Biopython**
    * **Numpy**
    * **pylab**
    * **matplotlib**
- ClustalW locally installed is required to run this script. Remember to change the clustalw path to your local path in [modules.py](https://github.com/claudefa ... ).  
- Internet connection is needed to perform Blast online. 

##Installing
To install the mirrorTree...

INSTRUCTIONS FOR INSTALLING MIRRORTREE SCRIPT:
'sudo ./install.sh' --> root privileges are needed


##Organization of the code
This program is split in different modules and scripts. Here, you can find the reference of what is contained in each file.

#####[mirrorTree](https://github.com/claudefa ...)
The workflow of the program. In this script all functions needed are called. It also contains the argument parser. 

#####[functions.py](https://github.com/claudefa ...)
All the core, helper functions and classes are found here. The documentation of each function is available in this file. 

#####[modules.py](https://github.com/claudefa ...)
In this files you can find all modules that the program needs to run. From this file they are imported to the others. 

##How to execute mirrorTree

This script is executed in the command line as following:

FOR RUNNING FROM THE COMMAND LINE:
'mirrorTree -h' --> This will show you how to execute the script



