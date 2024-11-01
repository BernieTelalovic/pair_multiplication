# pair.py
import warnings
from .young_diagram import YoungDiagram
from . import utils
import itertools as it

class Pair(YoungDiagram):
    
    def __init__(self, partition_tuple, Nc = None, weight:int = 1, inherited_N0:int=0):
        
        warnings.warn('Class Pair under construction. Proceed with caution!')
        
        permB,permA = partition_tuple
        self.partition  = partition_tuple
        self.barred = False
        
        self.Nc = Nc
        
        diagA = YoungDiagram(permA, Nc = Nc)
        diagB = YoungDiagram(permB,barred=True, Nc = Nc)
        
        self.N0 = max([diagA.N0+diagB.N0,inherited_N0])
        
        if (not Nc is None):
            if Nc < self.N0:
                weight = 0
                warnings.warn('Young diagram pair not admissible under given Nc.'+\
                              ' Weight/multiplicity will be set to 0.')
        
        self.pair = (diagB,diagA)
        self.word = [diagB.word,diagA.word]
        self.weight = weight
        
        
    def get_str(self):
        
        return '1_{'+str(self.N0)+r'}\left('+self.pair[0].get_str()+','+self.pair[1].get_str()+r'\right)'
        
    def __str__(self):
        
        return self.get_str()
    
    def __repr__(self):
        
        return '['+str(self.N0)+']('+self.pair[0].__repr__().replace('_','')+','+self.pair[1].__repr__()+')'
    
    def _repr_latex_(self):
        # Return LaTeX-formatted equation for Jupyter Notebook display

        return f"$$ {self.get_str()} $$"
        
    def __mul__(self, other):
        
        warnings.warn('Pair multiplication not yet implemented! Do not trust me!')
        
        if type(other) is Pair:
            
            mulA = self.pair*other.pair[0]
            
            mulB = self.pair*other.pair[1]
            
            return mulA+mulB
        
        elif type(other) is YoungDiagram:
            
            list_pairs = None
            
            if other.barred:
                list_pairs = self.pair[1]*other
            else:
                list_pairs = self.pair[0]*other
                            
        else:
            try:
                check_number = float(other)
                
                self.weight *= other
            except Exception as e:
                print(e)
                     
                    
    def set_Nc(self, va:int):
        
        self.Nc = val
        
    def multiplicity(self,val):
        return np.heaviside(val-self.N0,1)
        
    def admissible_under_Nc(Nc):
        
        return self.multiplicity(Nc)>0
    
    def conjugate(self,Nc=None):
        
        permA_original = np.array(self.partition[1]).astype(int)
        
        permA = np.array(extend_partition(permA_original,Nc)).astype(int)
        
        Nc = self.check_Nc(Nc)
        
        diagB_conj_partition = np.array(extend_partition(
                                    self.pair[0].conjugate(Nc = Nc).partition,Nc)
                                         ).astype(int)

        
        if np.all(diagB_conj_partition[0:len(permA_original)] == np.max(diagB_conj_partition)):
            
            new_perm = tuple(diagB_conj_partition.astype(int)+permA.astype(int))
            return YoungDiagram(new_perm,weight=self.weight, Nc=Nc, barred=False)
        
        else:
            warnings.warn('Conjugate not admissible under given Nc.')
            return NullDiagram(Nc)
