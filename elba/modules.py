
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Seq
import sys
import re
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
# from Bio.Phylo.TreeConstruction import DistanceCalculator
# from Bio import AlignIO
import numpy

# import sys
import os
import os.path
import argparse
# from Bio.Align.Applications import ClustalwCommandline
# from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio.Phylo.Consensus import *
# import numpy
from numpy import corrcoef, arange
import math
from pylab import show, title, savefig
import matplotlib.pyplot as plt