#!/usr/bin/python3

#############################################################################################
##################### MODULES NEEDED FOR EXECUTING 'mirrortree.py' ##########################
#############################################################################################

path_clustal = "/usr/bin/clustalw"

import sys, re, numpy, os, argparse, math, os.path

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import AlignIO
from Bio import Seq
from Bio.Align.Applications import ClustalwCommandline

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.Consensus import *

from pylab import show, title, savefig
import matplotlib.pyplot as plt

from numpy import corrcoef, arange

