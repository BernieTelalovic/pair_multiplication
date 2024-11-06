# young_diagram.py
import warnings
import numpy as np
import scipy as sp
from .utils import *
import itertools as it

class YoungDiagram:
    
    def __init__(self, partition, Nc = None, barred:bool = False, weight:int = 1):
        
        #if not type(partition) is np.array:
        #partition=list(np.array([partition]).flatten())
        #partition=tuple(partition)
        
        partition = partition_tuplify(partition)
        
        if len(partition)>1:
            if not all(i >= j for i, j in zip(partition, partition[1:])):
                raise ValueError('Not a young diagram.')
        
        self.barred = barred
        self.Nc = Nc
        N0 = len(partition)
        
        if (not Nc is None):
            if Nc < N0:
                weight = 0
                warnings.warn('Young diagram not admissible under given Nc. Weight/multiplicity will be set to 0.')
            elif Nc==N0:
                new_perm = tuple(np.trim_zeros(np.array(partition).astype(int)-np.array(partition).astype(int)[-1]))
                N0 = len(new_perm)
                partition = new_perm
                
        self.N0 = N0 
        self.weight = weight
        
        self.n = sum(partition)
        
        self.partition = partition
        self.word = self.get_word()
        
        self.hook_length = self.hook_length()
        
    def add_weight(self,added_weight):
        
        try:
            check_number = float(added_weight)
                
            self.weight += added_weight
        except Exception as e:
                print(e)
        
    def __eq__(self,other):
        
        if type(other)==type(self):
            if (self.Nc == other.Nc) and (self.N0==other.N0) and (self.partition==other.partition) \
                and (self.barred==other.barred):
                return True
            else:
                return False
        else:
            return False
        
    def get_str(self):
        
        barred_left = ''
        barred_right = ''
        if self.barred:
            barred_left = r'\overline{'
            barred_right = r'}'
            
        partition_str = str(self.partition).replace(',)', ')')
        #display(Math(r"E = mc^2"))
        
        return barred_left+partition_str+barred_right
        
    def get_cmdline_str(self):
    
        brd = ''
        if self.barred:
            brd='_'
        return str(self.partition).replace(',)', ')')+brd
        
    def __hash__(self):
        return hash((self.partition,self.N0,self.Nc,self.barred,self.weight))
            
    def __repr__(self):
        return self.get_cmdline_str()
        
    def __str__(self):
        return self.get_cmdline_str()
    
    def _repr_latex_(self):
        # Return LaTeX-formatted equation for Jupyter Notebook display

        return f"$$ {self.get_str()} $$"
        
    def __copy__(self):
        
        return self.deepcopy()
    
    def __deepcopy__(self):
        
        return YoungDiagram(self.partition, Nc=self.Nc, barred=self.barred, weight = self.weight)
        
    def __add__(self,other):
        
        if isinstance(other,YoungDiagram):
        
            elements = [self,other]
            multiplicities = [self.weight,other.weight]
            
            return DirectSum(elements,multiplicities)
        else:
            try:
                return other+self
            except Exception as e:
                raise NotImplemented

            
    def __radd__(self,other):
        return self+other
        
    def __mul__(self, other):
        
        Nc = None
        try:
            NcA = self.Nc
            NcB = other.Nc

            if not NcA==NcB:

                raise ValueError('The two diagrams must at least have equal Nc.')

            elif type(NcA)==int and type(NcB)==int:
                Nc = NcA
                warnings.warn('Diagram multiplication performed under specific Nc.')
        except exception as e:
            pass
        
        finally:
            if type(other) is Pair:

                mulA = self*other.pair[0]

                mulB = self*other.pair[1]

                return mulA+mulB

            elif type(other) is YoungDiagram:

                if other.barred==self.barred:
                    return self.multiplyLR(other,Nc)
                else:
                    return self.multiplyQ(other,Nc)
                
            elif type(other) is NullDiagram:
                return other
            
            else:
                try:
                    check_number = float(other)

                    self.weight *= other
                except Exception as e:
                    raise e
                    
    def __rmul__(self, other):
    
        return self*other
                
    def check_Nc(self,Nc):
        
        if Nc is None:
            Nc = self.Nc
            
        if Nc is None:
            raise ValueError('Nc>0 must be given as argument or the diagram must know it.')
            
        return Nc
                
    def conjugate(self, Nc = None, remember_Nc:bool=True, remember_mult:bool=False):
        
        Nc = self.check_Nc(Nc)
        new_diag = NullDiagram(Nc)
        
        if Nc < self.N0:
            warnings.warn('Young diagram pair not admissible under given Nc.'+\
                           'Returning null diagram.')
            return new_diag
        
        perm = extend_partition(self.partition,Nc)
        
        new_perm = tuple(np.trim_zeros(np.flip((np.ones(Nc)*np.max(perm)).astype(int) -\
                                               np.array(perm).astype(int))))
        
        if not remember_Nc:
            Nc = None
        weight = self.weight
        if not remember_mult:
            weight = 1
        
        new_diag = YoungDiagram(new_perm, Nc=Nc, weight=weight, barred = not self.barred)
        
        return new_diag
                
                
    def multiplyQ(self,other,Nc):
    
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            tmax=min([self.n,other.n])
            
            wordA = self.word
            wordB = other.word
            if self.barred:
                wordB = self.word
                wordA = other.word
                
            combinations = []
            weights = []
            for t in range(tmax+1):
                tcombs = superpose(t,wordB,wordA)
                
                combinations.append([Pair(tuple(ky), Nc=Nc) for ky in tcombs.keys()])
                weights.append([tcombs[ky] for ky in tcombs.keys()])
            combinations = list(it.chain(*combinations))
            weights = list(it.chain(*weights))
            return DirectSum(combinations,weights)
    
    
    def multiplyLR(self,other,Nc):
    
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            results_dict = LR_multiple(self.partition,other.partition)
            
            diag_results = [YoungDiagram(tuple(ky), Nc = Nc, barred = self.barred) for ky in results_dict.keys()]
            diag_weights = list(results_dict.values())
            
            return DirectSum(diag_results,diag_weights)
    
    
    def get_word(self,imaginary = True):

        partition = self.partition
        word = 0+0j
        if len(self.partition)>0:

            word = np.array(list(it.chain(*[[ind+1]*partition[ind] 
                                    for ind in range(len(partition))])),dtype=complex)

            if imaginary:
                word+= 1j*np.concatenate([np.arange(1,partition[ind]+1) for ind in range(len(partition))])

        return word
    
    def set_barred(self, val: bool):
    
        self.barred = val
        
    def set_Nc(self,val:int):
        
        self.Nc = val
        
    def bar(self):
        
        self.barred = not self.barred
        
    def multiplicity(self,val):
        return np.heaviside(val-self.N0,1)
    
    def set_N0(self,val):
        
        if val < self.N0:
            warnings.warn('Cannot lower N0 for this diagram.')
        self.N0=max([val,self.N0])
        
    def hook_length(self):
    
        if self.barred:
        
            return None
        else:
            modarr = make_2d_array_ones(self.partition)
            
            hook = np.array([[np.sum(modarr[indx,indy:])+np.sum(modarr[indx:,indy])-1 for 
                            indy in range(len(modarr[0]))]
                            for indx in range(len(modarr))])
            return np.prod(hook[hook>0])
            
    def dimension_Sn(self):
    
        if self.barred:
            warnings.warn('Barred diagrams currently do not support Sn dimentions.')
            return None
        else:
            return sp.special.factorial(np.sum(self.partition))/self.hook_length
            
    def dimension_Nc(self,Nc=None):
    
        Nc = self.check_Nc(Nc)
        if self.barred:
        
            new_diag = self.conjugate(Nc=Nc)
            return new_diag.dimension_Nc()
        else: 
            if len(self.partition)>0:
                
                modarr = make_2d_array_ones(self.partition)
               
                index_arr = ((modarr*np.arange(0,max(self.partition))).T - np.arange(0,len(self.partition))).T
                
                return int(np.prod(Nc + index_arr[modarr>0])/self.hook_length)
            else:
                return 1
        
            
    def evaluate_for_Nc(self,Nc = None):
    
        Nc = self.check_Nc(Nc)
        if self.barred:
            return self.conjugate(Nc = Nc)
        else:
            new_diag = YoungDiagram(self.partition, barred = self.barred,Nc = Nc, weight = self.weight)
            return new_diag
            
    def pair_with(self,other):
            
        if isinstance(other, YoungDiagram):
        
            if self.barred and not other.barred:
                return Pair((self,other),Nc=self.Nc, inherited_N0=max([self.N0,other.N0]))
            elif other.barred and not self.barred:
                return Pair((other,self),Nc=self.Nc, inherited_N0=max([self.N0,other.N0]))
            else:
                raise AttributeError("One of the diagrams must be barred to form a pair")
        else:
            try:
                other_partition = partition_tuplify(other)
                
                other_diag = YoungDiagram(other_partition,barred = not self.barred)
                
                return self.pair_with(other_diag)
                
            except Exception as e:
                raise TypeError('I can only be paired with another diagram or a partition.')
        
            
            
        
        
        
        
#######################################################################################################################
############################################NullDiagram################################################################
#######################################################################################################################
    
class NullDiagram(YoungDiagram):
    
    def __init__(self,Nc):
        
        YoungDiagram.__init__(self,(), Nc = Nc, weight = 0)
        self.partition = None
        self.barred = False
        
    def multiplicity(self,val):
        return 0*val
        
    def dimension_Nc(self,Nc=None):
        return 0
        

#######################################################################################################################
##########################################Pair#########################################################################
#######################################################################################################################

class Pair(YoungDiagram):
    
    def __init__(self, partition_tuple, Nc = None, weight:int = 1, inherited_N0:int=0):
        
        warnings.warn('Class Pair under construction. Proceed with caution!')
        
        permB = partition_tuple[0]
        permA = partition_tuple[1]
        self.barred = False
        
        self.Nc = Nc
        
        diagA = permA
        diagB = permB
        if not isinstance(diagA,YoungDiagram):
            diagA = YoungDiagram(permA, Nc = Nc)
        if not isinstance(diagB,YoungDiagram):
            diagB = YoungDiagram(permB,barred=True, Nc = Nc)
            
        self.partition = (diagB.partition, diagA.partition)
        
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
        
    def get_cmdline_str(self):
    
        return '['+str(self.N0)+']('+self.pair[0].__repr__().replace('_','')+','+self.pair[1].__repr__()+')'
        
    def __str__(self):
        
        return self.get_cmdline_str()
    
    def __repr__(self):
        
        return self.get_cmdline_str()
    
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
            
            if other.barred:
                return barred_diagram_mult(other)
            else:
                return unbarred_diagram_mult(other)
                            
        else:
            try:
                check_number = float(other)
                
                self.weight *= other
            except Exception as e:
                print(e)
                
    def __rmul__(self, other):
        return self*other
        
        
    def __add__(self,other):
        
        if isinstance(other,Pair):
        
            elements = [self,other]
            multiplicities = [self.weight,other.weight]
            
            return DirectSum(elements,multiplicities)
        else:
            try:
                return other+self
            except Exception as e:
                raise NotImplemented

            
    def __radd__(self,other):
        return self+other
                     
                    
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
        
        if Nc < self.N0:
            warnings.warn('Conjugate not admissible under given Nc.')
            return NullDiagram(Nc)
        
        diagB_conj_partition = np.array(extend_partition(
                                    self.pair[0].conjugate(Nc = Nc).partition,Nc)
                                         ).astype(int)

        
        if np.all(diagB_conj_partition[0:len(permA_original)] == np.max(diagB_conj_partition)):
            
            new_perm = tuple(diagB_conj_partition.astype(int)+permA.astype(int))
            return YoungDiagram(new_perm,weight=self.weight, Nc=Nc, barred=False)
        
        else:
            warnings.warn('Conjugate not admissible under given Nc.')
            return NullDiagram(Nc)
            
    def evaluate_for_Nc(self,Nc = None):
    
        return self.conjugate(Nc=Nc)
        
    def dimension_Nc(self,Nc=None):
    
        diag_from_pair = self.conjugate(Nc=Nc)
        
        return diag_from_pair.dimension_Nc()
            
            
#######################################################################################################################
##########################################DirectSum####################################################################
#######################################################################################################################        
            
            
class DirectSum(dict):
    
    def __init__(self, keys, values):
        # Ensure keys and values are of the same length
        if len(keys) != len(values):
            raise ValueError("List of diagrams must have a corresponding list of multiplicities.")
        
        ky_arr = np.array(keys)
        vl_arr = np.array(values)
        
        container = np.array([[ky, np.sum(vl_arr[ky_arr==ky])] for ky in list(set(ky_arr))])
        
        
        # Use the dict constructor to initialize the dictionary with key-value pairs
        super().__init__(zip(container.T[0], container.T[1]))

        
    def conjugate(self,Nc):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            keys = self.keys()
            values = list(self.values())
            
            new_keys = np.array([ky.conjugate(Nc=Nc) for ky in keys])
            new_vals = np.array([values[ind]*new_keys[ind].multiplicity(Nc) 
                                 for ind in range(len(new_keys))])
            
            return DirectSum(new_keys[new_vals>0], new_vals[new_vals>0].astype(type(Nc)))
        
    def evaluate_for_Nc(self,Nc):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            keys = self.keys()
            values = list(self.values())
            
            new_keys = np.array([ky.evaluate_for_Nc(Nc=Nc) for ky in keys])
            new_vals = np.array([values[ind]*new_keys[ind].multiplicity(Nc) 
                                 for ind in range(len(new_keys))])
            
            return DirectSum(new_keys[new_vals>0], new_vals[new_vals>0].astype(type(Nc)))
        
    def dimension_Nc(self,Nc=None):
    
        elements = [ky.dimension_Nc(Nc) for ky in self.keys()]
        multiplicities = list(self.values())
            
        return DimensionDirectSum(elements,multiplicities)
        
    def get_str(self):
    
        strin = ''
        elements = list(self.keys())
        multiplicities = list(self.values())
        oplus = r'\oplus'
        
        for ind in range(len(elements)):
            el_str = elements[ind].get_str()
            mult = str(int(multiplicities[ind]))
            
            if not '1_{' in el_str:
                el_str = mult+r'_{'+str(int(elements[ind].N0))+'}\,'+el_str
            else:
                el_str = el_str.replace('1_{', mult+r'_{')
                
            if ind == len(elements)-1:
                oplus = ''
                
            strin += el_str + oplus
        return strin
        
    def get_cmdline_str(self):
    
        strin = ''
        elements = list(self.keys())
        multiplicities = list(self.values())
        oplus = '+'
        
        for ind in range(len(elements)):
            
            mult = str(int(multiplicities[ind]))
            el_str = mult+elements[ind].get_cmdline_str()
                
            if ind == len(elements)-1:
                oplus = ''
                
            strin += el_str + oplus
        return strin
        
    def __str__(self):
        
        return self.get_cmdline_str()
    
    def _repr_latex_(self):            

        return f"$$ {self.get_str()} $$"
        
        
    def __add__(self,other):
        
        keys = []
        values = []
        
        if type(other)==DirectSum:
            keys = list(self.keys())+list(other.keys())
            values = list(self.values())+list(other.values())
            
        elif type(other) in [YoungDiagram,NullDiagram,Pair]:
            
            keys = list(self.keys())
            keys.append(other)
            
            values = list(self.values())
            values.append(other.weight)
        
        return DirectSum(keys,values)
        
    def __iadd__(self,other):
    
        added = self+other
        self = added
            
    def __radd__(self,other):
        return self+other
    
    def set_N0(self,val):
        
        new_keys = [ky.set_N0(val) for ky in self.keys()]
        
        self = DirectSum(new_keys,list(self.values()))
        
    def __mul__(self, other):
    
        warnings.warn('DirectSum multiplication only works with scalars at the moment! Proceed with caution!')
        
        if type(other) is DirectSum:
        
            return compose_direct_sums()
        
        elif type(other) is Pair:
            
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
                
    def __rmul__(self, other):
        return self*other
                
    def elements(self):
        return list(self.keys())
    def multiplicities(self):
        return list(self.values())
    def lowest_Nc(self):
        return [el.N0 for el in self.keys()]
          




class DimensionDirectSum(DirectSum):

    def __init__(self, keys, values):
        # Ensure keys and values are of the same length
        if len(keys) != len(values):
            raise ValueError("List of dimensions must have a corresponding list of multiplicities.")

        ky_arr = np.array(keys)
        vl_arr = np.array(values)

        super().__init__(ky_arr[ky_arr>0], vl_arr[ky_arr>0])

    def get_cmdline_str(self):
        strin = ''
        elements = list(self.keys())
        multiplicities = list(self.values())
        oplus = r'+'
        
        for ind in range(len(elements)):
            el_str = str(int(elements[ind]))
            mult = str(int(multiplicities[ind]))
            
            el_str = mult+'\u00D7'+el_str
                
            if ind == len(elements)-1:
                oplus = ''
                
            #if int(elements[ind])>0:
            strin += el_str + oplus
        return strin
    
    def get_str(self):
    
        strin = ''
        elements = list(self.keys())
        multiplicities = list(self.values())
        oplus = r'+'
        
        for ind in range(len(elements)):
            el_str = str(int(elements[ind]))
            mult = str(int(multiplicities[ind]))
            
            el_str = mult+r'\cdot'+el_str
                
            if ind == len(elements)-1:
                oplus = ''
                
            #if int(elements[ind])>0:
            strin += el_str + oplus
        return strin
        
    def __str__(self):
        
        return self.get_cmdline_str()
    
    def _repr_latex_(self):            

        return f"$$ {self.get_str()} $$"
        
    def sum(self):
    
        dims = np.array(list(self.keys()))
        mults = np.array(list(self.values()))
        
        return np.sum(dims*mults).astype(int)
        
    

