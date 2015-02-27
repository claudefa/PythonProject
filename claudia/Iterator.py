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
        

def Protein_creator(fin):
    fd=open(fin, "r")
    i=0
    for line in fd:
        if i==0:
            i+=1
        elif i==1:
            sp=line[14:-4]
            i+=1
        elif i==2:
            seq_id=line.strip()[8:]
            i+=1
        elif i==3:
            i+=1
        elif i==4:
            e_val=line.strip()[9:]
            i+=1
        elif i==5:
            score=line.strip()[7:]
            i+=1
        elif i==6:
            identity=line.strip()[10:]
            i+=1
        elif i==7:
            i+=1
        elif i==8:
            seq=line.strip()
            yield Protein(seq_id, seq, sp, e_val, score, identity)
            i=0
            
if __name__ == "__main__":
    for prot in Protein_creator("1COW.out.blast"):
        print(prot)



        
        
    
    
    