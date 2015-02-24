class Sequence(object):

	mw = {}
	alphabet = ''

	def __init__(self, identifier, sequence):
		self.__identifier = str(identifier)
		for s in sequence:	
			if not s in list(self.alphabet):
				raise ValueError("Impossible to create instance: %s is not possible" %s)
			else:
				self.__sequence = str(sequence)
				
	def get_identifier(self):
		return self.__identifier
	
	def get_sequence(self):
		return self.__sequence
	
	def get_mw(self):
		weight = 0
		for aa in self.__sequence:
			weight += self.mw.get(aa, 0.0)
		return weight
			
	def has_subsequence(self, sequence_obj):
		return sequence_obj.get_sequence() in self.__sequence		
	
	def get_length(self):
		return len(self.__sequence)


class NucleotideSequence(Sequence):	
	table = {}
	start_codons = []
	stop_codons = []
	
	def translate(self):
		seq = self.get_sequence()
		codons = []
		start = ''
		stop = ''
		defi = []
		for i in range(0, len(seq),3):
			codon = seq[i:i+3]
			if codon in self.start_codons:
				start = i
			if codon in self.stop_codons:
				stop = i
		codons = seq[start:stop]	
		for i in range(0, len(codons),3):
			triplet = codons[i:i+3]
			defi.append(triplet)
		aminoacids = ''
		for codon in defi:
			aminoacids += self.table.get(codon.upper())
		return ProteinSequence(self.get_identifier(),aminoacids)


class DNASequence(NucleotideSequence):
	alphabet = 'GATC'
	mw = {'A': 347.0, 'C': 323.0, 'T': 322.0, 'G': 363.0}
	table = {'CTT': 'L', 'ATG': 'M', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'AGC': 'S', 'AGA': 'R', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'ACT': 'T', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'TAC': 'Y', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GAC': 'D', 'GAA': 'E', 'AAG': 'K', 'AAA': 'K', 'AAC': 'N', 'CTC': 'L', 'CAT': 'H', 'AAT': 'N', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'TGT': 'C', 'TCT': 'S', 'GAT': 'D', 'TTT': 'F', 'TGC': 'C', 'TGG': 'W', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TCA': 'S', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A'}
	start_codons = ['TTG', 'CTG', 'ATG']
	stop_codons = ['TAA', 'TAG', 'TGA']
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

	def transcribe(self):
		seq = self.get_sequence()
		compl = ''
		for nt in seq:
			compl += self.complement.get(nt)
		rna = compl.replace('T', 'U')
		return RNASequence(self.get_identifier(),rna)
	
	
class RNASequence(NucleotideSequence):
	alphabet = 'GAUC'
	mw = {'A': 363.0, 'C': 339.0, 'U': 340.0, 'G': 379.0}
	table = {'GUC': 'V', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'AAC': 'N', 'AGG': 'R', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'AGA': 'R', 'AAU': 'N', 'ACU': 'T', 'CAC': 'H', 'GUG': 'V', 'CCG': 'P', 'CCA': 'P', 'AGU': 'S', 'CCC': 'P', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'CGA': 'R', 'CAG': 'Q', 'CGC': 'R', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'CCU': 'P', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GAG': 'E', 'UCC': 'S', 'UAC': 'Y', 'CGU': 'R', 'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'UCA': 'S', 'AUG': 'M', 'CUG': 'L', 'AUU': 'I', 'CAU': 'H', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', 'GUA': 'V', 'UGC': 'C', 'GCU': 'A', 'UGU': 'C', 'CUC': 'L', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', 'UUC': 'F'}
	stop_codons = ['UAA', 'UAG', 'UGA']
	start_codons = ['UUG', 'CUG', 'AUG']
	complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}

	def reverse_transcribe(self):
		seq = self.get_sequence()
		compl = ''
		for nt in seq:
			compl += self.complement.get(nt)
		cdna = compl.replace('U','T')
		return DNASequence(self.get_identifier(),cdna)
			
	
class ProteinSequence(Sequence):
	mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
	alphabet = 'ACDEFGHIKLMNPQRSTVWY'
	pass