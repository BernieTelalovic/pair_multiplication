# young_diagram.py
import warnings
import numpy as np
import scipy as sp
from .utils import *
from .cutils_pair_multiplication import *
from .tableau_classes import *
import itertools as it

#######################################################################################################################
############################################NullDiagram################################################################
#######################################################################################################################
    
class NullDiagram:

    """
        Represents a null/inadmissible Young diagram. Returned when trying to evaluate a Young diagram or Pair for an Nc that is too low.
        
        Attributes:
            partition: None
            weight: 0 
            barred: False
            
        Returns 0 for multiplicity and dimension.
            
    """


    
    def __init__(self,Nc=None):
        
        self.partition = None
        self.N0 = 0
        self.Nc = None
        self.barred = False
        self.weight = 0
        self.n = 0
        #self.hook_length = 0
        self.width = 0
        self._hash = 0

        
    def multiplicity(self,val):
        return 0*val
        
    def dimension_Nc(self,Nc=None):
        return 0        
        
    def get_str(self):
        
        barred_left = ''
        barred_right = ''
        if self.barred:
            barred_left = r'\overline{'
            barred_right = r'}'
            
        partition_str = str(clean_numpy(self.partition)).replace(',)', ')')
        
        strn = barred_left+partition_str+barred_right
        
        if self.partition is None:
            strn = '0'
        
        return strn
        
    def get_cmdline_str(self):
    
        brd = ''
        if self.barred:
            brd='_'
            
        strn = str(clean_numpy(self.partition)).replace(',)', ')')+brd
        
        if self.partition is None:
            strn = '0'
        
        return strn

        
    def __hash__(self):
        return self._hash
            
    def __repr__(self):
        return self.get_cmdline_str()
        
    def __str__(self):
        return self.get_cmdline_str()
    
    def _repr_latex_(self):
        # Return LaTeX-formatted equation for Jupyter Notebook display

        return f"$$ {self.get_str()} $$"
        
    def __add__(self,other):
    
        if isinstance(other,NullDiagram):
        
            return other
        
        elif isinstance(other,YoungDiagram) or isinstance(other,Pair):
        
            elements = [other]
            multiplicities = [other.weight]
            
            return DirectSum(elements,multiplicities)
            
        elif isinstance(other,DirectSum):

            return other
        else:
            try:
                return other+self
            except Exception as e:
                raise NotImplemented
                
    def __sub__(self,other):
    
        if self==other:
        
            return NullDiagram()#DirectSum([self],[0])
        
        elif isinstance(other,YoungDiagram) or isinstance(other,Pair) or isinstance(other,DirectSum):
        
            return other

        else:
            try:
                return other+self
            except Exception as e:
                raise NotImplemented     

            
    def __radd__(self,other):
        return self+other
        
        
    def __mul__(self, other):
        
        return NullDiagram()
                    
    def __rmul__(self, other):
    
        return self*other
        
    def __lt__(self,other):
    
        try:
            if self == other:
                return False
            elif self.N0 < other.N0:
                return True
            elif self.N0 == other.N0:
                if self.n < other.n:
                    return True
                elif self.n == other.n:
                
                    if self.width < other.width:
                        return True
                    else:
                        return np.all(self.get_hook_length() < other.get_hook_length())
                else:
                    return False
            else:
                return False

        except Exception as e:
            raise NotImplemented
            
    def __le__(self,other):
    
        try:
        
            lt = self < other
            eq = self == other
            return lt + eq

        except Exception as e:
            raise NotImplemented
            
    def __ge__(self,other):
    
        try:
        
            gt = not (self < other)
            eq = self == other
            return gt + eq

        except Exception as e:
            raise NotImplemented
            
    def __ne__(self,other):
    
        try:
            eq = self == other
            return not eq

        except Exception as e:
            raise NotImplemented
            
    def __eq__(self,other):
        
        if type(other)==type(self):
            return hash(self) == hash(other)
        elif ((other == 0) or (other is None)) and self.partition is None:
            return True
        elif type(other)==NullDiagram:
            return False
            
        elif type(other)==Pair:
            dex = 0 if other.barred else 1
            return hash(other[dex]) == hash(self) and (other[dex-1].n == 0)
            
        elif type(other)==YoungDiagram:
            dex = 0 if other.barred else 1
            return hash(self.pair[dex]) == hash(other) and (self.pair[dex-1].n == 0)
            
        elif type(other)==DirectSum:
            equal = DirectSum([self],[self.weight])==other
            return equal
        else:
            return False
            
    
    def simplify(self):
        return self
        
    def LR_multiply(self, other):
        return self
        
    def evaluate_for_Nc(self,Nc):
        return self

    def lowest_Nc(self):
        return self.N0
        
    def print(self, tex = False):
        strin = ''
        if tex:
            strin = self._repr_latex_()
        else:
            strin = self.get_cmdline_str()
            
        print(strin)
        
    def to_str(self,tex=False):
        strin = ''
        if tex:
            strin = self._repr_latex_()
        else:
            strin = self.get_cmdline_str()
            
        return strin


#######################################################################################################################
############################################YoungDiagram###############################################################
#######################################################################################################################


class YoungDiagram(NullDiagram):

    """
    A subclass of NullDiagram. Represents a young diagram with labelled with a given partition, 
        and provides methods to calculate its dimension given a certain Nc, 
        its representation given a certain Nc and handles multiplication by 
        other Young diagram objects (barred and unbarred).

        Attributes:
            partition (tuple): the partition labelling the diagram. Same as the one given in the constructor, or reduced if the longest column is equal to the given Nc, only if an Nc is given in the constructor.
            N0 (int): the first Nc where this diagram appears. If the given Nc is None, N0 is the number of rows in the diagram.
            Nc (int or None): the Nc under which the diagram appears. 
            n (int): the number of boxes in the Young diagram
            hook_length: the hook length of the Young diagram. Set to None if the diagram is barred.
            [internal] word (list): a list of complex entries with the real parts labelling the rows and the imaginary parts labelling the coulmns of each box in the diagram. Used for multiplication with other Young diagrams.
            [internal] weight (int): the multiplicity of the diagram. Used for removing diagram objects from direct sums when they are not admissible under given Nc. 
            """
            
    def __new__(cls, partition, Nc = None, barred:bool = False, weight:int = 1,inherited_N0:int=0):
    
        if partition is None:
            return NullDiagram()
        elif not (Nc is None):
            if Nc < len(partition_tuplify(partition)):
                return NullDiagram()
            else:
                return super(YoungDiagram, cls).__new__(YoungDiagram)
        else:
            return super(YoungDiagram, cls).__new__(YoungDiagram)
    
    def __init__(self, partition, Nc = None, barred:bool = False, weight:int = 1,inherited_N0:int=0):
        
        """
        Constructs a Young diagram from a given partition
        
        Parameters:
            partition (list,tuple): the partition labelling the diagram.
            Nc (int, default None): the Nc under which the diagram is constructed. If not given, the diagram is constructed independent of Nc.
            barred (bool, default: False): whether the diagram is barred or not.
            [internal] weight (int, default 1): the multiplicity of the diagram. Used for removing diagram objects from direct sums when they are not admissible under given Nc. 

        Raises:
            ValueError: If the given partition is not monotonic decreasing, i.e., does not label a Young diagram.

        """
        if isinstance(self, YoungDiagram):
            super().__init__()

            self.weight = 1
            self.width = 0
            partition = partition_tuplify(partition)
            
            if len(partition)>1:
                self.width = partition[0]
                if not all(i >= j for i, j in zip(partition, partition[1:])):
                    raise ValueError('Not a young diagram.')
            elif len(partition) > 0:
                self.width = partition[0]
            
            
            N0 = len(partition)
            
            if not (Nc is None):
                if (Nc==N0) and (N0>0):
                    new_perm = tuple(np.trim_zeros(np.array(partition).astype(int)-np.array(partition).astype(int)[-1]))
                    N0 = len(new_perm)
                    partition = new_perm
                    
            self.inherited_N0 = inherited_N0
            self.barred = barred
            self.Nc = Nc
            self.N0 = N0
            
            self.n = sum(partition)
            
            self.partition = partition
            self.word = self.get_word()
            self._hash = hash((self.partition,self.N0,self.Nc,self.barred,self.weight))
            
            #self.hook_length = self.get_hook_length()

            
    def _null_me(self):
        self = NullDiagram()
        
    def __eq__(self,other):

        if type(other)==type(self):
            return hash(self) == hash(other)
        elif type(other)==NullDiagram:
            return False
            
        elif type(other)==Pair:
            dex = 0 if other.barred else 1
            return hash(other[dex]) == hash(self) and (other[dex-1].n == 0)
            
        elif type(other)==YoungDiagram:
            dex = 0 if other.barred else 1
            return hash(self.pair[dex]) == hash(other) and (self.pair[dex-1].n == 0)
            
        elif type(other)==DirectSum:
            equal = DirectSum([self],[self.weight])==other
            return equal
        else:
            return False
    
        
    def __copy__(self):
        
        return self.deepcopy()
    
    def __deepcopy__(self):
        
        return YoungDiagram(self.partition, Nc=self.Nc, barred=self.barred, weight = self.weight,inherited_N0=self.N0)
        
    def __hash__(self):
        return self._hash
        
    def __mul__(self, other):
        
        Nc = None
        try:
            NcA = self.Nc
            NcB = other.Nc

            if not NcA==NcB:

                raise ValueError('The two diagrams must at least have equal Nc.')

            elif type(NcA)==int and type(NcB)==int:
                Nc = NcA
                #warnings.warn('Diagram multiplication performed under specific Nc.')
        except exception as e:
            pass
        
        finally:
            if type(other) is NullDiagram:
                return other
            if type(other) is Pair:
            
                sl_pair = self.pair_with(YoungDiagram((), barred = not self.barred))
                mult = other*sl_pair
                
                if not (Nc is None):
                    mult = mul.evaluate_for_Nc(Nc)
                
                return mult

            elif type(other) is YoungDiagram:
            
                sl_pair = self.pair_with(YoungDiagram((), barred = not self.barred))
                oth_pair = other.pair_with(YoungDiagram((), barred = not other.barred))
                
                mul = oth_pair*sl_pair
                
                if (not Nc is None):
                    mul = mul.evaluate_for_Nc(Nc)
                
                return mul
            
            else:
                try:
                    check_number = float(other)

                    return DirectSum([self], [self.weight*other])
                except Exception as e:
                    raise NotImplemented
                    
    def __rmul__(self, other):
    
        return self*other
        
    def __add__(self,other):
        
        if isinstance(other,YoungDiagram) or isinstance(other,Pair):
        
            elements = [self,other]
            multiplicities = [self.weight,other.weight]
            
            return DirectSum(elements,multiplicities)
        elif isinstance(other,NullDiagram):
            return self
        else:
            try:
                return other+self
            except Exception as e:
                raise NotImplemented

            
    def __radd__(self,other):
        return self+other
        
    def __sub__(self,other):
        
        if isinstance(other,DirectSum):
        
            elements = other.elements()+[self]
            multiplicities = other.multiplicities()+[-self.weight]
            
            return DirectSum(elements,multiplicities)
        elif isinstance(other,NullDiagram):
            return self
        else:
            if self==other:
                return NullDiagram()
            else:
                raise NotImplemented

            
    def __rsub__(self,other):
        return self+other
        
                
    def check_Nc(self,ncc, allow_none = False):
        
        Nc = ncc
        if ncc == None:
            Nc = self.Nc
            
        if Nc == None:
            if not allow_none:
                raise ValueError('Nc>0 must be given as argument or the diagram must know it.')
        
        return Nc
                
    def conjugate(self, Nc = None, remember_Nc:bool=True, remember_mult:bool=False):
        
        Nc = self.check_Nc(Nc)
        new_diag = NullDiagram(Nc)
        
        
        if Nc < self.N0:
            warnings.warn('Young diagram pair not admissible under given Nc. Returning null diagram.')
            return new_diag
        
        if len(self.partition) == 0:
            
            return YoungDiagram((), Nc=Nc, weight=self.weight, barred = not self.barred)
            
        else:
            perm = extend_partition(self.partition,Nc)
            
            new_perm = tuple(np.trim_zeros(np.flip((np.ones(Nc)*np.max(perm)).astype(int) -\
                                                   np.array(perm).astype(int))))
            
            if not remember_Nc:
                Nc = None
            weight = self.weight
            if not remember_mult:
                weight = 1
            
            new_diag = YoungDiagram(new_perm, Nc=Nc, weight=weight, barred = not self.barred)
            if new_diag.weight == 0:
                return NullDiagram(Nc=Nc)
            
            return new_diag
                
                
    def multiplyQ(self,other,Nc):
    
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            tmax=min([self.n,other.n])
            n0min = max(self.N0,other.N0)
            
            wordA = self.word
            wordB = other.word
            if self.barred:
                wordB = self.word
                wordA = other.word
                
            combinations = []
            weights = []
            for t in range(tmax+1):
                tcombs,cands, counts = superpose(t,wordB,wordA)
                
                combinations.append([Pair(tuple(ky), Nc=Nc, inherited_N0=n0min) for ky in tcombs.keys()])
                weights.append([tcombs[ky] for ky in tcombs.keys()])
            combinations = list(it.chain(*combinations))
            weights = list(it.chain(*weights))
            return DirectSum(combinations,weights)
    
    
    def multiplyLR(self,other,Nc):
    
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            n0min = max(self.N0,other.N0)
            
            results_dict,combs,cands = LR_multiple(self.partition,other.partition)
            
            diag_results = [YoungDiagram(tuple(ky), Nc = Nc, barred = self.barred, inherited_N0=n0min) 
                            for ky in results_dict.keys() 
                            ]
                            
            acceptable_diag = Nc
            if Nc == None:
                acceptable_diag = np.inf
                            
            diag_weights = list(results_dict.values())
            diag_weights *= np.array([diag.N0 <= acceptable_diag for diag in diag_results])
            
            return DirectSum(diag_results,diag_weights)
    
    
    def get_word(self,imaginary = True):

        partition = self.partition
        word = 0+0j
        if len(self.partition)>0:

            word = np.array(list(it.chain(*[[ind+1]*partition[ind] 
                                            for ind in range(len(partition))])),
                                    dtype=complex)

            if imaginary:
                word+= 1j*np.concatenate([np.arange(1,partition[ind]+1) 
                                          for ind in range(len(partition))])

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
        
    def get_hook_length(self):
    
        modarr = make_2d_array_ones(self.partition)
            
        hook = np.array([[np.sum(modarr[indx,indy:])+np.sum(modarr[indx:,indy])-1 
                            for indy in range(len(modarr[0]))]
                            for indx in range(len(modarr))])
        return np.prod(hook[hook>0])
            
    def dimension_Sn(self):
    
        if self.barred:
            warnings.warn('Barred diagrams currently do not support Sn dimentions.')
            return None
        else:
            return sp.special.factorial(np.sum(self.partition))/self.get_hook_length()
            
    def dimension_Nc(self,Nc=None):
    
        Nc = self.check_Nc(Nc)
        if self.barred:
        
            new_diag = self.conjugate(Nc=Nc)
            return new_diag.dimension_Nc()
        else: 
            if len(self.partition)>0:
                
                modarr = make_2d_array_ones(self.partition)
               
                index_arr = ((modarr*np.arange(0,max(self.partition))).T - np.arange(0,len(self.partition))).T
                
                return int(np.prod(Nc + index_arr[modarr>0])/self.get_hook_length())
            else:
                return 1
        
            
    def evaluate_for_Nc(self,Nc = None):
    
        Nc = self.check_Nc(Nc)
        if self.barred:
            return self.conjugate(Nc = Nc)
        else:
            new_diag = YoungDiagram(self.partition, barred = self.barred, Nc = Nc, weight = self.weight,inherited_N0=self.N0)
            if new_diag.weight == 0:
                return NullDiagram(Nc=Nc)
            return new_diag
            
    def pair_with(self,other, inherited_N0 = 0,Nc = None):
            
        if isinstance(other, YoungDiagram):
        
            if self.barred and not other.barred:
                return Pair((self,other),Nc=None, 
                            inherited_N0=max(self.N0,other.N0, inherited_N0))
            elif other.barred and not self.barred:
                return Pair((other,self),Nc=None, 
                            inherited_N0=max(self.N0,other.N0, inherited_N0))
            else:
                raise AttributeError("One of the diagrams must be barred to form a pair")
        else:
            try:
                other_partition = partition_tuplify(other)
                
                other_diag = YoungDiagram(other_partition,barred = not self.barred)
                
                return self.pair_with(other_diag)
                
            except Exception as e:
                raise TypeError('I can only be paired with another diagram or a partition.')
                
    def LR_multiply(self,other):
    
        Nc = None
        try:
            NcA = self.Nc
            NcB = other.Nc

            if not (NcA==NcB):

                raise ValueError('The two diagrams must at least have equal Nc.')

            elif type(NcA)==int and type(NcB)==int:
                Nc = NcA
                #warnings.warn('Diagram multiplication performed under specific Nc.')
        except Exception as e:
            pass
            
        Nc = self.Nc
    
        if (other == 0) or (self == 0):
            return NullDiagram()
        elif type(other) is YoungDiagram:
            
                if other.barred==self.barred:
                
                    if np.all(other.word == 0j):
                        return DirectSum([self], [self.weight*other.weight])
                    elif np.all(self.word == 0j):
                        return DirectSum([other], [self.weight*other.weight])
                    else:
                        return self.multiplyLR(other,Nc)
                else:
                    if np.all(other.word == 0j) or np.all(self.word == 0j):
                        return DirectSum([self.pair_with(other)], [self.weight*other.weight])
                    return self.multiplyQ(other,Nc)
            
        else:
            try:
                check_number = float(other)

                self.weight *= other
            except Exception as e:
                raise e
        
    
    def as_ydiagram(self):
        pr = self.pair_with(())
        return pr._as_inner_ytab()
        

class BarredDiagram(YoungDiagram):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, barred=True, **kwargs)  # Force `barred=True`




#######################################################################################################################
##########################################Pair#########################################################################
#######################################################################################################################

class Pair(YoungDiagram):

    """
    A subclass of YoungDiagram. Composes a barred and unbarred diagram into a pair.

        Attributes:
            partition (tuple): the partition labelling the pair in the form ((barred_diagram,unbarred_diagram)). 
            N0 (int): the first Nc where this diagram appears. If the given Nc is None, N0 is the sum of the number of rows in the barred and unbarred diagram.
            Nc (int or None): the Nc under which the diagram appears. 
            [Not used] n (int): the number of boxes in the Young diagram
            [Not used] hook_length: the hook length of the Young diagram. Set to None if the diagram is barred.
            [internal] word (list): in the form [word of barred diagram, word of unbarred diagram]
            [internal] weight (int): the multiplicity of the pair. Used for removing diagram objects from direct sums when they are not admissible under given Nc. 
    """
    def __new__(cls, partition_tuple, Nc = None, weight:int = 1, inherited_N0:int=0):
    
        if not Nc is None:
            permA_original = np.array(partition_tuplify(partition_tuple[1])).astype(int)
        
            permA = np.array(extend_partition(permA_original,Nc),dtype=object).astype(int)
            diagB = YoungDiagram(partition_tuple[0],barred=True)
            
            diagB_conj_partition = np.array(extend_partition(
                                            diagB.conjugate(Nc = Nc).partition,Nc),dtype=object
                                             ).astype(int)

            
            if np.all(diagB_conj_partition[0:len(permA_original)] == np.max(diagB_conj_partition,initial=0)):
                
                new_perm = tuple(diagB_conj_partition.astype(int)+permA.astype(int))
                
                return YoungDiagram(new_perm, Nc = Nc)
            else:
                return NullDiagram()

        else:
            return object.__new__(cls)

    
    def __init__(self, partition_tuple, Nc = None, weight:int = 1, inherited_N0:int=0):
        
        #warnings.warn('Class Pair under construction. Proceed with caution!')
        if isinstance(self, Pair):
            super(YoungDiagram).__init__()
            self.partition = None
            self.N0 = 0
            self.Nc = None
            self.barred = False
            self.weight = 0
            self.n = 0
            #self.hook_length = 0
            permB = partition_tuple[0]
            permA = partition_tuple[1]
            
            diagA = permA
            diagB = permB
            if not isinstance(diagA,YoungDiagram):
                diagA = YoungDiagram(permA)
            if not isinstance(diagB,YoungDiagram):
                diagB = YoungDiagram(permB,barred=True)
                
            if (diagA.weight == 0) or (diagB.weight == 0):
                weight = 0
                return
            
            N0 = max([len(diagA.partition)+len(diagB.partition),inherited_N0])
                    
            self.N0 = N0
            self.n = diagA.n+diagB.n
            self.partition = (diagB.partition, diagA.partition)
            
            self.pair = (diagB,diagA)
            self.word = [diagB.word,diagA.word]
            self.width = diagB.width+diagA.width
            self.weight = weight
            self._hash = hash((self.partition,self.N0,self.Nc,self.barred,self.weight))
            #self.hook_length = self.get_hook_length()

        
    def get_str(self):
        strn = ''
        if len(self.partition[0])==0:
            strn = self.pair[1].get_str()
        elif len(self.partition[1])==0:
            strn = self.pair[0].get_str()
        else:
            strn = r'\left('+self.pair[0].get_str()+','+self.pair[1].get_str()+r'\right)'
            
        if len(strn)==0:
            strn = '0'
        else:
            strn = '1_{'+str(self.N0)+r'} ' + strn
        
        return strn
        
    def get_cmdline_str(self):
        strn = ''
        if len(self.partition[0])==0:
            strn = self.pair[1].__repr__()
        elif len(self.partition[1])==0:
            strn = self.pair[0].__repr__()
        else:
            strn = '('+self.pair[0].__repr__().replace('_','')+','+self.pair[1].__repr__()+')'
        if len(strn)==0:
            strn = '0'
        else:
            strn = '['+str(self.N0)+']' + strn
        
        return strn
        
    def __str__(self):
        
        return self.get_cmdline_str()
    
    def __repr__(self):
        
        return self.get_cmdline_str()
    
    def _repr_latex_(self):
        # Return LaTeX-formatted equation for Jupyter Notebook display

        return f"$$ {self.get_str()} $$"
        
    def __mul__(self, other):
        
        Nc = None
        try:
            NcA = self.Nc
            NcB = other.Nc

            if not NcA==NcB:

                raise ValueError('The two diagrams must at least have equal Nc.')

            elif type(NcA)==int and type(NcB)==int:
                Nc = NcA
                #warnings.warn('Diagram multiplication performed under specific Nc.')
        except exception as e:
            pass
        finally:
            if type(other) is Pair:
            
                if self.n == 0:
                    return DirectSum([other], [other.weight])
                if other.n == 0:
                    return DirectSum([self], [self.weight])
                
                candidates = pair_multipy(self.partition,other.partition)
                tuples = candidates.get_partitions()
                n0s = candidates.n0s
                pair_arr = []
                if len(tuples)>0:
                    pair_arr = [Pair(tuples[ind],inherited_N0=n0s[ind], Nc = Nc) for ind in range(len(n0s))]
        
                return DirectSum(pair_arr,np.ones(len(pair_arr)))
            
            elif type(other) is YoungDiagram:
            
                if self.n==0:
                    return DirectSum([other.pair_with(())], [other.weight])
                if other.n ==0:
                    return DirectSum([self], [self.weight])
                
                oth = other.pair_with(())
                
                candidates = pair_multipy(self.partition,oth.partition)
                tuples = candidates.get_partitions()
                n0s = candidates.n0s
                pair_arr = []
                if len(tuples)>0:
                    pair_arr = [Pair(tuples[ind],inherited_N0=n0s[ind], Nc = Nc) for ind in range(len(n0s))]
        
                return DirectSum(pair_arr,np.ones(len(pair_arr)))
                
            elif type(other) is NullDiagram:
                return NullDiagram()
                                
            else:
                try:
                    check_number = float(other)
                    
                    return DirectSum([self], [self.weight*other])
                except Exception as e:
                    raise NotImplemented
                
                
    def __rmul__(self, other):
        return self*other
                     
                    
    def set_Nc(self, va:int):
        
        self.Nc = val
        
    def multiplicity(self,val):
        return np.heaviside(val-self.N0,1)
        
    def admissible_under_Nc(Nc):
        
        return self.multiplicity(Nc)>0
    
    def conjugate(self,Nc=None):
        
        Nc = self.check_Nc(Nc)
        
        if Nc < self.N0:
            #warnings.warn('Conjugate not admissible under given Nc.')
            return NullDiagram(Nc)
        
        permA_original = np.array(partition_tuplify(self.partition[1])).astype(int)
        permA = np.array(extend_partition(permA_original,Nc),dtype=object).astype(int)
        
        diagB_conj_partition = np.array(extend_partition(
                                    self.pair[0].conjugate(Nc = Nc).partition,Nc),dtype=object
                                         ).astype(int)

        
        if np.all(diagB_conj_partition[0:len(permA_original)] == np.max(diagB_conj_partition,initial=0)):
            
            new_perm = tuple(diagB_conj_partition.astype(int)+permA.astype(int))
            diag = YoungDiagram(new_perm,weight=self.weight, Nc=Nc, barred=False)
            return diag
        
        else:
            warnings.warn('Conjugate not admissible under given Nc.')
            return NullDiagram(Nc)
            
    def evaluate_for_Nc(self,Nc = None):
    
        return self.conjugate(Nc=Nc)
        
    def dimension_Nc(self,Nc=None):
    
        diag_from_pair = self.conjugate(Nc=Nc)
        
        return diag_from_pair.dimension_Nc()
        
    def get_hook_length(self):
    
        return np.array([self.pair[0].get_hook_length(), self.pair[1].get_hook_length()])
        
    def n(self,Nc):
    
        return self.evaluate_for_Nc(Nc).n
        
    def simplify(self):
    
        permB = self.partition[0]
        permA = self.partition[1]
        
        if len(permB) == 0:
            return YoungDiagram(permA, Nc = self.Nc, weight = self.weight, inherited_N0=self.N0)
        elif len(permA) == 0:
            return YoungDiagram(permB,barred = True, Nc = self.Nc, weight = self.weight, inherited_N0=self.N0)
        else:
            return self
            
    def _as_inner_ytab(self):
    
        ubr = self.partition[1]
        br = self.partition[0]
        
        if len(ubr)==0 and len(br)==0:
            return r'\bullet'
        
        num_ubr_rows = len(ubr)
        num_br_cols = 0
        br_str = ''
        if len(br) > 0:
            num_br_cols = max(br)
            br_str = '0,'*num_ubr_rows+','.join([f"{num_br_cols-row}+{row}" for row in list(reversed(br))])
        
        ubr_str = ','.join([f"{num_br_cols}+{row}" for row in ubr])
        
        ubr_part = ''
        br_part = ''
        if len(ubr_str)>0:
            ubr_part = "[*(white)]{"+ubr_str+"}"
        if len(br_str)>0:
            br_part = r"[*(white)\bullet]{"+br_str+"}"
            
        if len(ubr_part)>0 and len(br_part)>0:
            return "\ydiagram"+br_part+'*'+ubr_part
        elif len(ubr_part)==0:
            return "\ydiagram"+br_part
        elif len(br_part)==0:
            return "\ydiagram"+ubr_part
            
    def as_ydiagram(self):
        return self._as_inner_ytab()
           
    def _as_ytab(self,name = None):
    
        name = fix_name(self.partition,name)
        
        inner_tab = self._as_inner_ytab()
        
        strin = r"\newcommand{"+"\\"+name+"}{"+inner_tab+"}"
        return strin
            
    def print_as_ytab(self,name = None):
        
        print(self._as_ytab(name = name))
        
        
    def verbosely_multiply_with(self, other, name = None):
    
        name = fix_name(self.partition,name)
        strin_start = r"\newcommand{"+"\\"+name+"}{\n"
        strin_end = "}"
        Nc = None
        try:
            NcA = self.Nc
            NcB = other.Nc

            if not NcA==NcB:

                raise ValueError('The two diagrams must at least have equal Nc.')

            elif type(NcA)==int and type(NcB)==int:
                Nc = NcA
                #warnings.warn('Diagram multiplication performed under specific Nc.')
        except exception as e:
            pass
        finally:
            if type(other) is Pair:
                    
                left = r'\left(\ '
                right = r'\ \right)'
                
                other_array = enumerate_pair(other.partition)
                
                tot_string = self._as_inner_ytab()+ r'\,\otimes\,'+compose_ytab_from_array(other_array)
                
                candidates_in_steps = pair_multipy_verbosely(self.partition,other.partition)
                
                for ind in range(len(candidates_in_steps)):
                
                    cand_str = candidates_in_steps[ind].string_list()
                    leftover_arr = compose_ytab_from_array(other_array[:,1+ind:])
                    if len(leftover_arr)>0:
                        cand_str = left+cand_str+right+r'\,\otimes\,'+leftover_arr+r'\\'
                    tot_string += r'&='+cand_str
                    
                    if ind < len(candidates_in_steps)-1:
                        tot_string +='\n'
                
                return tot_string
            
            elif type(other) is YoungDiagram:
                
                oth = other.pair_with(())
                
                left = r'\left(\ '
                right = r'\ \right)'
                
                other_array = enumerate_pair(oth.partition)
                
                tot_string = self._as_inner_ytab()+ r'\,\otimes\,'+compose_ytab_from_array(other_array)
                
                candidates_in_steps = pair_multipy_verbosely(self.partition,oth.partition)
                
                for ind in range(len(candidates_in_steps)):
                
                    cand_str = candidates_in_steps[ind].string_list()
                    leftover_arr = compose_ytab_from_array(other_array[:,1+ind:])
                    if len(leftover_arr)>0:
                        cand_str = left+cand_str+right+r'\,\otimes\,'+leftover_arr+r'\\'
                    tot_string += r'&='+cand_str
                    
                    if ind < len(candidates_in_steps)-1:
                        tot_string +='\n'
                
                return tot_string
                
            elif type(other) is NullDiagram:
                return('0')
                                
            else:
                try:
                    check_number = float(other)
                    
                    self.weight *= other
                except Exception as e:
                    print(e)
        
    

        
#######################################################################################################################
##########################################DirectSum####################################################################
####################################################################################################################### 
        
class DirectSum:

    def __new__(cls, elements=None, multiplicities=None):

        if elements is None:
            return NullDiagram()
        else:
            if len(elements)>0:
                if isinstance(multiplicities,list) or isinstance(multiplicities,type(np.array([1]))):
                    if len(elements)==len(multiplicities):
                        if np.any(np.array(multiplicities) > 0):
                        
                            ky_arr = np.array(elements)
                            vl_arr = np.array(multiplicities)
                            
                            container = np.array([[ky, np.sum(vl_arr[ky_arr==ky])] 
                                                    for ky in list(set(ky_arr)) 
                                                    if getattr(ky,'weight',1)!=0 and np.sum(vl_arr[ky_arr==ky])!=0]).T
                            
                            if len(container)>0:
                            
                                instance = super(DirectSum, cls).__new__(DirectSum)
                                instance._elements = container[0]
                                instance._multiplicities = container[1]
                                return instance 
                            else:
                                return NullDiagram()
                        else:
                            return NullDiagram()
                    else:
                        instance = super(DirectSum, cls).__new__(DirectSum)
                        instance._elements = []
                        instance._multiplicities = []
                        
                        return instance 
                        
                else:
                    
                    vl_arr = np.ones(len(elements))
                    ky_arr = np.array(elements)
                    
                    container = np.array([[ky, np.sum(vl_arr[ky_arr==ky])] 
                                            for ky in list(set(ky_arr)) 
                                            if getattr(ky,'weight',1)!=0 and np.sum(vl_arr[ky_arr==ky])!=0]).T
                                            
                    if len(container)>0:
                    
                        instance = super(DirectSum, cls).__new__(DirectSum)  # Create an instance of the class
                        instance._elements = container[0]
                        instance._multiplicities = container[1]

                        return instance 
                    else:
                        return NullDiagram()
            else:
                return NullDiagram()

    
    def __init__(self, elements=None, multiplicities=None):

        if hasattr(self, "_elements") and hasattr(self, "_multiplicities"):
        
            if len(self._multiplicities) == 0:
                raise ValueError("List of multiplicities must have equal length to list of diagrams/pairs.")
            else:
                elements = self._elements
                multiplicities = self._multiplicities
                del self._elements  # Cleanup temporary attributes
                del self._multiplicities
                
        
        elements = np.array(elements,dtype=object)
        inds = elements.argsort()
        multiplicities = np.array(multiplicities,dtype = int)[inds]
        elements = elements[inds]
        n0s = np.array([el.N0 if hasattr(el, "N0") else 0 for el in elements])
            
        self.elements = elements
        self.multiplicities = multiplicities
        self.N0 = n0s

    def __getitem__(self, condition):
        """
        Allow NumPy-like filtering to return a new DirectSum object containing only
        the elements that satisfy the condition.
        """
        if isinstance(condition, np.ndarray):
            filtered_elements = self.elements[condition]
            filtered_multiplicities = self.multiplicities[condition]
            return DirectSum(filtered_elements, filtered_multiplicities)
        else:
            raise TypeError("Indexing must be done with a boolean NumPy array.")
            
            
            
    def conjugate(self,Nc):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            values = self.multiplicities
            
            new_keys = np.array([ky.conjugate(Nc=Nc) for ky in self.elements])
            new_vals = np.array([values[ind]*new_keys[ind].multiplicity(Nc) 
                                 for ind in range(len(new_keys))])
            
            return DirectSum(new_keys[new_vals>0], new_vals[new_vals>0].astype(type(Nc)))
        
    def evaluate_for_Nc(self,Nc):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            values = self.multiplicities
            
            new_keys = np.array([ky.evaluate_for_Nc(Nc=Nc) for ky in self.elements])
            new_vals = np.array([values[ind]*new_keys[ind].multiplicity(Nc) 
                                 for ind in range(len(new_keys))])
                                 
            new_keys_full = new_keys[new_vals>0]
            new_values_full = new_vals[new_vals>0].astype(type(Nc))
            
            if len(new_keys_full)>0:
                return DirectSum(new_keys_full, new_values_full)
            else:
                return NullDiagram()
        
    def dimension_Nc(self,Nc=None):
    
        elements = np.array([ky.dimension_Nc(Nc) for ky in self.elements],dtype = float)
        #multiplicities = np.array(self.multiplicities)[elements>0]
        #elements = elements[elements>0]
        
        return elements
            
        #return DimensionDirectSum(elements,multiplicities)
        
    def sum_dimensions(self, Nc=None):
        dims = self.dimension_Nc(Nc)
        
        return (dims*self.multiplicities).sum()
        
    def get_str(self):
    
        strin = ''
        lines = []  # List of lines for LaTeX output
        current_line = []  # Objects in the current line
        elements = self.elements
        multiplicities = self.multiplicities.astype(int)
        
        oplus = r'\oplus'
        N0_prev = 0#
        if len(elements)> 0:
            N0_prev = int(self.N0[0])
        
        for ind in range(len(elements)):
            N0 = int(self.N0[ind])
            if N0 > N0_prev:
                #start_str += "\,$$ \n $$"
                if current_line:
                    lines.append(''.join(current_line))
               
                N0_prev = N0
                current_line = []
                
            el_str = elements[ind].get_str()
            mult = str(int(multiplicities[ind]))
            
            if not '1_{' in el_str:
                el_str = mult+r'_{'+str(N0)+'}\,'+el_str
            else:
                el_str = el_str.replace('1_{', mult+r'_{')
                
            if ind == len(elements)-1:
                oplus = ''
                
            current_line.append(el_str + oplus)
            
            if ind == len(elements)-1:
                lines.append(''.join(current_line))
                
        if len(lines)==0:
            return '0'

        latex_string = r"\[\begin{array}{c}" + "\n\n".join(lines) + r"\end{array}\]"
        return latex_string
        
    def get_cmdline_str(self):
    
        strin = ''
        elements = self.elements
        multiplicities = self.multiplicities
        oplus = '+'
        
        for ind in range(len(elements)):
            
            mult = str(int(multiplicities[ind]))
            el_str = mult+elements[ind].get_cmdline_str()
                
            if ind == len(elements)-1:
                oplus = ''
                
            strin += el_str + oplus
        if len(strin)==0:
            strin = '0'
        
        return strin
        
    def __str__(self):
        
        return self.get_cmdline_str()
        
    def __repr__(self):
        return self.get_cmdline_str()
    
    def _repr_latex_(self):            

        return self.get_str()
        
        
    def __add__(self,other):
        
        elements = []
        multiplicities = []
        
        if type(other)==DirectSum:
            elements = list(self.elements)+list(other.elements)
            multiplicities = list(self.multiplicities)+list(other.multiplicities)
            return DirectSum(elements,multiplicities)
        elif type(other) in [YoungDiagram,NullDiagram,Pair]:
            
            elements = list(self.elements)
            elements.append(other)
            
            multiplicities = list(self.multiplicities)
            multiplicities.append(other.weight)
            return DirectSum(elements,multiplicities)
        else:
            raise NotImplemented

        
    def __iadd__(self,other):
    
        added = self+other
        self = added
            
    def __radd__(self,other):
        return self+other
        
    def __sub__(self,other):
        
        elements = []
        multiplicities = []
        
        if type(other)==DirectSum:
            elements = list(self.elements)+list(other.elements)
            multiplicities = list(self.multiplicities)+list(-1*np.array(list(other.multiplicities)))
            return DirectSum(elements,multiplicities)
        elif type(other) in [YoungDiagram,NullDiagram,Pair]:
            
            elements = list(self.keys())
            elements.append(other)
            
            multiplicities = list(self.multiplicities)
            multiplicities.append(-1*other.weight)
            
            return DirectSum(elements,multiplicities)
        else:
            raise NotImplemented
            
        
        
    def __rsub__(self,other):
        return self-other
        
    def __isub__(self,other):
    
        added = self-other
        self = added
    
    def set_N0(self,val):
        
        new_keys = [ky.set_N0(val) for ky in self.elements]
        
        self = DirectSum(new_keys,self.multiplicities)
        
    def __mul__(self, other):
    
        warnings.warn('DirectSum multiplication only works with scalars at the moment! Proceed with caution!')
        
        if type(other) is DirectSum:
            #TODO
            keys_me = self.elements
            keys_other = other.elements
            
            mults_me = self.multiplicities
            mults_other = other.multiplicities
            
            container = [(keys_me[mind]*keys_other[ond])*mults_me[mind]*mults_other[ond]
                         for mind in range(len(keys_me)) for ond in range(len(keys_other))]
            
            return np.sum(container)
        
        elif (type(other) is Pair) or (type(other) is YoungDiagram):
            keys_me = self.elements
            mults_me = self.multiplicities
            
            container = [(keys_me[mind]*other)*mults_me[mind] for mind in range(len(keys_me))]
            return np.sum(container)
            
        elif other==0:
            
            return NullDiagram
 
        else:
            try:
                check_number = float(other)
                
                values = np.array(self.multiplicities)*other
                keys = self.elements
                
                return DirectSum(keys,values)
                
            except Exception as e:
                raise NotImplemented
                
    def __rmul__(self, other):
        return self*other
                
    def keys(self):
        return self.elements
    def values(self):
        return self.multiplicities
    def N0(self):
        return self.N0
    def lowest_Nc(self):
        return self.N0
        
    def __eq__(self,other):
    
        try:
            if (self-other)==0:
                return True
            else:
                return False
        except Exception as e:
            return False
        
    def as_ydiagram(self):
    
        elements = self.elements
        multiplicities = self.multiplicities.astype(int)
        n0s = self.N0

        lis = [str(multiplicities[ind])+r'_{'+str(n0s[ind])+'}\,'+elements[ind].as_ydiagram()
                for ind in range(len(elements))]
    
        return '\,\oplus\,'.join(lis)

    def print(self, tex = False):
        strin = ''
        if tex:
            strin = '$$'+self._repr_latex_()+'$$'
        else:
            strin = self.get_cmdline_str()
            
        print(strin)
        
    def to_str(self,tex=False):
        strin = ''
        if tex:
            strin = self._repr_latex_()
        else:
            strin = self.get_cmdline_str()
            
        return strin

class DimensionDirectSum:

    def __init__(self, keys, values):
        # Ensure keys and values are of the same length
        if len(keys) != len(values):
            raise ValueError("List of dimensions must have a corresponding list of multiplicities.")

        #ky_arr = np.array(keys)
        #vl_arr = np.array(values)
        
        #container = np.array([[keys[ind],values[ind]] for ind in range(len(keys)) if keys[ind]>0]).T
        
        #container = np.array([[ky, np.sum(values[keys==ky])] 
        #                        for ky in list(set(keys)) 
        #                        if np.sum(values[keys==ky])!=0 and ky>0]).T
        
        elements = np.array(keys)#np.array(container[0])
        inds = elements.argsort()
        multiplicities = np.array(values)[inds]
        elements = elements[inds]
        self.elements = elements
        self.multiplicities = multiplicities

    def get_cmdline_str(self):
        strin = ''
        elements = self.elements
        multiplicities = self.multiplicities
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
        elements = self.elements
        multiplicities = self.multiplicities
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
    
        dims = self.elements
        mults = self.multiplicities
        
        return np.sum(dims*mults).astype(int)
        
    def lowest_Nc(self):
        raise AttributeError("DimensionDirectSum only carries dimension info and cannot recover partitions or lowest Nc.")
        
    def N0(self):
        return self.lowest_Nc()
        
    

