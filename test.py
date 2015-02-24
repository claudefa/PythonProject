from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Seq

handle = open("thrombin.fa")
for record in SeqIO.parse(handle, "fasta") :
	result = NCBIWWW.qblast("blastp", "nr", record.seq)
handle.close()

for n in NCBIXML.parse(result):
	for alignment in n.alignments:
		for hsp in alignment.hsps:
#			if hsp.expect < 0.1:
				print ("****Alignment****")
				print ("sequence:", alignment.title)
				print ("length:", alignment.length)
				print ("E-value:", hsp.expect)
				print (hsp.query[:50] + "...")
				print (hsp.match[:50] + "...")
				print (hsp.sbjct[:50] + "...")
