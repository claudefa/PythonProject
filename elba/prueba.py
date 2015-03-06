import re

class Protein(object):
	def __init__(self, seq_id, seq, sp, e_val, score, identity):
		self.seq_id=seq_id
		self.seq=seq
		self.sp=sp
		self.e_val=e_val
		self.score=score
		self.identity=identity
    
	def __str__(self):
		return ("%s\n%s\n%s\n%s\n%s\n%s\n" %(self.seq_id, self.seq, self.sp, self.e_val, self.score, self.identity))
	def get_specie (self):
		return str(self.sp)
	def get_id (self):
		return self.seq_id
	def get_seq (self):
		return self.seq

def Protein_creator(fin):
	fd=open(fin, "r")
	info = []
	p = re.compile("\[(.*)\]")
	for field in fd:
		if "#" in field:
			if not info == []:
				yield Protein(info[1], info[6], info[0], info[3], info[4], info[5])
				info = []
		else:
			if "[" in field:
				sp = p.search(field).group(1)
				info.append(sp)
			else:
				line = field.split()
				info.append(line[1])			
	yield Protein(info[1], info[6], info[0], info[3], info[4], info[5])
	fd.close()

for element in Protein_creator("prueba.blast"):
	print(element.get_id())
