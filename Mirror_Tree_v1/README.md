# MIRRORTREE 1.0

This is the PYT and SBI project by [**Elba Raimúndez**](https://github.com/elbaraim), [**Clàudia Fontserè**](https://github.com/claudefa) and **Lucas Michel**!

- [What is mirrorTree 1.0?](#description)
- [Requirements before installing](#requirements)
- [Installation](#installation)
- [Organization of the code](#codeorganization)
    * [mirrorTree](#mirrortree)
    * [setup.py](#setup)
    * [modules.py](#modules)
- [How to execute mirrorTree 1.0?](#execute)

##What is mirrorTree 1.0?
**mirrorTree 1.0** is a program to predict protein-protin interaction using similarity between two protein families (Pazos and Valencia, 2001).

It can be initiated from 2 protein sequences (in the same fasta file with fa/fasta extension) and from 2 multiple sequence alignments (two files with aln extension).

The **workflow** of mirrorTree is the following:
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
- Beware that this script is written in **Python 3.4.2**. Check your version before executing mirrorTree.
- Required python modules
    * [**Biopython**](http://biopython.org/)
    * [**Numpy**](http://www.numpy.org/)
    * [**Scipy and PyLab**](http://www.scipy.org/)
    * [**matplotlib**](http://matplotlib.org/)
- [**ClustalW**](http://www.clustal.org/) locally installed is required to run this script. Remember to change the clustalw path to your local path in [modules.py](https://github.com/claudefa/PythonProject/blob/master/Mirror_Tree_v1/MirrorTree/modules.py).  
- Internet connection is needed to perform [**BLAST**](http://blast.ncbi.nlm.nih.gov/Blast.cgi). 

##Installation

The **instructions** for the proper installation of mirrorTree 1.0 are the following:
 - Remember! Make sure your ClustalW path is modified in [modules.py](https://github.com/claudefa/PythonProject/blob/master/Mirror_Tree_v1/MirrorTree/modules.py).
 - Root privileges are needed for installation.
 - The following commands should be called in the command line inside the directory /install/:
 * `python3 setup.py build`
 * `sudo python3 setup.py install`

##Organization of the code
mirrorTree 1.0 package contains:

This program is split in different modules and scripts. Here, you can find the reference of what is contained in each file.

#####[mirrorTree](https://github.com/claudefa/PythonProject/tree/master/Mirror_Tree_v1/MirrorTree)
The main scripts are found in this directory. 
- [**\_\_init\_\_.py**](https://github.com/claudefa/PythonProject/blob/master/Mirror_Tree_v1/MirrorTree/__init__.py): A python file required for a proper installation of the program.
- [**mirrorTree**](https://github.com/claudefa/PythonProject/blob/master/Mirror_Tree_v1/MirrorTree/mirrorTree): The workflow of the program. In this script all functions needed are called. It also contains the argument parser. 
- [**functions.py**](https://github.com/claudefa/PythonProject/blob/master/Mirror_Tree_v1/MirrorTree/functions.py): All the core, helper functions and classes are found here. The documentation of each function is available in this file. 
- [**modules.py**](https://github.com/claudefa/PythonProject/blob/master/Mirror_Tree_v1/MirrorTree/modules.py): In this files you can find the ClustalW path and all modules that the program needs to run. From this file they are imported to the others. 

#####[setup.py](https://github.com/claudefa/PythonProject/blob/master/Mirror_Tree_v1/setup.py)
The script needed to install the package. 


##How to execute mirrorTree

This script is executed in the command line as following:

If starting from the two protein sequences (fasta format with extension .fa/.fasta):
`mirrorTree -i input.fa`

If starting form the two alignments (clustalW alignment with extension .aln)
`mirrorTree -i input1.aln input2.aln`

There are other parameter:
- -v: verbose. To see more detailed stderr. 
- -s: save. To save tmp files.
- -o: output. To save data into a predefined output file. By default they are shown as stdout.
- -e: evalue. Change threshold of Blast filtering. By default it is set in 0.00001.
- -id: identity. Change threshold of Blast filtering. By default it is set in 30. 

If any doubt:
- `mirrorTree -h`
This will show you how to execute the script.



