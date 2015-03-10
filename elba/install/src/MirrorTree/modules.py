#!/usr/bin/python3

#############################################################################################
##################### MODULES NEEDED FOR EXECUTING 'mirrortree.py' ##########################
#############################################################################################

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

from numpy import corrcoef, arange
#from pylab import show, title, savefig
#import matplotlib.pyplot as plt


#from Bio import Phylo
# from Bio.Phylo.TreeConstruction import DistanceCalculator
# from Bio import AlignIO
#import numpy
#import math
# import sys
#import os
#import re
#import argparse
# from Bio.Align.Applications import ClustalwCommandline
# import numpy

