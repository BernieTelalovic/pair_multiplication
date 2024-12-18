# young_diagram.py
import warnings
import numpy as np
import scipy as sp
from .utils import *
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
        self.hook_length = 0

        
    def multiplicity(self,val):
        return 0*val
        
    def dimension_Nc(self,Nc=None):
        return 0
        
    def add_weight(self,added_weight):
        
        self.weight=0
        
        
    def get_str(self):
        
        barred_left = ''
        barred_right = ''
        if self.barred:
            barred_left = r'\overline{'
            barred_right = r'}'
            
        partition_str = str(self.partition).replace(',)', ')')
        
        strn = barred_left+partition_str+barred_right
        
        if len(strn)==0:
            strn = '0'
        
        return strn
        
    def get_cmdline_str(self):
    
        brd = ''
        if self.barred:
            brd='_'
            
        strn = str(self.partition).replace(',)', ')')+brd
        if len(strn)==0:
            strn = '0'
        
        return strn

        
    def __hash__(self):
        return hash((self.partition,self.N0,self.Nc,self.barred,self.weight))
            
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
        
            return DirectSum([self],[0])
        
        elif isinstance(other,YoungDiagram) or isinstance(other,Pair):
        
            elements = [other]
            multiplicities = [-other.weight]
            
            return DirectSum(elements,multiplicities)
            
        elif isinstance(other,DirectSum):

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
                    return np.all(self.hook_length < other.hook_length)
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
            if (self.Nc == other.Nc) and (self.N0==other.N0) and (self.partition==other.partition) \
                and (self.barred==other.barred):
                return True
            else:
                return False
        elif type(other)==YoungDiagram:
            if (self.Nc == other.Nc) and (self.N0==other.N0) and (self.partition==other.partition) \
                and (self.barred==other.barred):
                return True
            else:
                return False
        elif other is None:
            return True
        elif other == 0:
            return True
        else:
            return False
            
    
    def simplify(self):
        return self
        
    def LR_multiply(self):
        return self




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
        super().__init__()
        self.partition = None
        self.N0 = 0
        self.Nc = None
        self.barred = False
        self.weight = 0
        self.n = 0
        self.hook_length = 0
        partition = partition_tuplify(partition)
        
        if len(partition)>1:
            if not all(i >= j for i, j in zip(partition, partition[1:])):
                raise ValueError('Not a young diagram.')
        
        
        N0 = len(partition)
        
        if not (Nc is None):
            if Nc < N0:
                weight = 0
                warnings.warn('Young diagram not admissible under given Nc. Weight/multiplicity will be set to 0.')
                
                return
                
            elif Nc==N0:
                new_perm = tuple(np.trim_zeros(np.array(partition).astype(int)-np.array(partition).astype(int)[-1]))
                N0 = len(new_perm)
                partition = new_perm
                
        self.inherited_N0 = inherited_N0
        self.barred = barred
        self.Nc = Nc
        self.N0 = N0
        self.weight = weight
        
        self.n = sum(partition)
        
        self.partition = partition
        self.word = self.get_word()
        
        self.hook_length = self.get_hook_length()

            
    def _null_me(self):
        self = NullDiagram()
        
    def __eq__(self,other):

        if type(other)==type(self):
            if (self.Nc == other.Nc) and (self.N0==other.N0) and (self.partition==other.partition) \
                and (self.barred==other.barred):
                return True
            else:
                return False
        elif type(other)==NullDiagram:
            if (self.Nc == other.Nc) and (self.N0==other.N0) and (self.partition==other.partition) \
                and (self.barred==other.barred):
                return True
            else:
                return False
                
        elif type(other)==Pair:
            dex = 0 if self.barred else 1
            if (self.partition == other.partition[dex]) and ((self.N0 == other.N0) or self.inherited_N0 == other.N0)\
                and (len(other.partition[dex-1])==0):
                return True
            else:
                return False
        elif type(other)==YoungDiagram:
            dex = 0 if other.barred else 1
            if (other.partition == self.partition[dex]) and ((self.N0 == other.N0) or other.inherited_N0 == self.N0)\
                and (len(self.partition[dex-1])==0):
                return True
            else:
                return False
            
        else:
            return False

        
    def add_weight(self,added_weight):
        
        try:
            check_number = float(added_weight)
                
            self.weight += added_weight
        except Exception as e:
                print(e)
    
        
    def __copy__(self):
        
        return self.deepcopy()
    
    def __deepcopy__(self):
        
        return YoungDiagram(self.partition, Nc=self.Nc, barred=self.barred, weight = self.weight,inherited_N0=self.N0)
        
    def __hash__(self):
        return hash((self.partition,self.N0,self.Nc,self.barred,self.weight))
        
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
                sl_pair = self.pair_with(YoungDiagram((), barred = not self.barred),Nc = Nc)
                mult = other*sl_pair
                if not (Nc is None):
                    mult = mul.evaluate_for_Nc(Nc)
                
                return mult

            elif type(other) is YoungDiagram:
            
                sl_pair = self.pair_with(YoungDiagram((), barred = not self.barred),Nc = Nc)
                oth_pair = other.pair_with(YoungDiagram((), barred = not other.barred),Nc = Nc)
                
                mul = oth_pair*sl_pair
                if (not Nc is None):
                    mul = mul.evaluate_for_Nc(Nc)
                
                return mul
            
            else:
                try:
                    check_number = float(other)

                    self.weight *= other
                except Exception as e:
                    raise e
                    
    def __rmul__(self, other):
    
        return self*other
                
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
            new_diag = YoungDiagram(self.partition, barred = self.barred, Nc = Nc, weight = self.weight,inherited_N0=self.N0)
            if new_diag.weight == 0:
                return NullDiagram(Nc=Nc)
            return new_diag
            
    def pair_with(self,other, inherited_N0 = 0,Nc = None):
            
        if isinstance(other, YoungDiagram):
        
            if self.barred and not other.barred:
                return Pair((self,other),Nc=self.check_Nc(Nc, allow_none = True), 
                            inherited_N0=max([self.N0,other.N0, inherited_N0]))
            elif other.barred and not self.barred:
                return Pair((other,self),Nc=self.check_Nc(Nc, allow_none = True), 
                            inherited_N0=max([self.N0,other.N0, inherited_N0]))
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

            if not NcA==NcB:

                raise ValueError('The two diagrams must at least have equal Nc.')

            elif type(NcA)==int and type(NcB)==int:
                Nc = NcA
                #warnings.warn('Diagram multiplication performed under specific Nc.')
        except Exception as e:
            pass
    
        if type(other) is YoungDiagram:
            
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
                
        elif type(other) is NullDiagram:
                return other
            
        else:
            try:
                check_number = float(other)

                self.weight *= other
            except Exception as e:
                raise e
        
    
        




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

    
    def __init__(self, partition_tuple, Nc = None, weight:int = 1, inherited_N0:int=0):
        
        #warnings.warn('Class Pair under construction. Proceed with caution!')
        super(YoungDiagram).__init__()
        self.partition = None
        self.N0 = 0
        self.Nc = None
        self.barred = False
        self.weight = 0
        self.n = 0
        self.hook_length = 0
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
        
        if not (Nc is None):
            if Nc < N0:
                weight = 0
                return
                
        self.barred = False
        self.Nc = Nc
        self.N0 = N0
        self.n = diagA.n+diagB.n
        self.partition = (diagB.partition, diagA.partition)
        
        self.pair = (diagB,diagA)
        self.word = [diagB.word,diagA.word]
        self.weight = weight
        self.hook_length = self.get_hook_length()

        
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
            strn = '1_{'+str(self.N0)+r'}\,' + strn
        
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
                
                tuples,n0s = pair_multipy(self.partition,other.partition)
                pair_arr = [Pair(tuples[ind],inherited_N0=n0s[ind], Nc = Nc) for ind in range(len(n0s))]
        
                return DirectSum(pair_arr,np.ones(len(pair_arr)))
            
            elif type(other) is YoungDiagram:
                
                if other.barred:
                
                    oth = Pair((other.partition,()), inherited_N0=other.N0, Nc = Nc)
                    return self*oth
                else:
                    oth = Pair(((),other.partition), inherited_N0=other.N0, Nc = Nc)
                    return self*oth
                                
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
        
        Nc = self.check_Nc(Nc)
        
        if Nc < self.N0:
            warnings.warn('Conjugate not admissible under given Nc.')
            return NullDiagram(Nc)
        
        permA_original = np.array(self.partition[1]).astype(int)
        
        permA = np.array(extend_partition(permA_original,Nc)).astype(int)
        
        diagB_conj_partition = np.array(extend_partition(
                                    self.pair[0].conjugate(Nc = Nc).partition,Nc)
                                         ).astype(int)

        
        if np.all(diagB_conj_partition[0:len(permA_original)] == np.max(diagB_conj_partition)):
            
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
    
        return np.array([self.pair[0].hook_length, self.pair[1].hook_length])
        
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
        
        container = np.array([[ky, np.sum(vl_arr[ky_arr==ky])] 
                                for ky in list(set(ky_arr)) 
                                if getattr(ky,'weight',1)>0 and np.sum(vl_arr[ky_arr==ky])!=0])
        
        
        
        # Use the dict constructor to initialize the dictionary with key-value pairs
        if len(container.T)>0:
            super().__init__(zip(container.T[0], container.T[1]))
        else:
            super().__init__()

        
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
                                 
            new_keys_full = new_keys[new_vals>0]
            new_values_full = new_vals[new_vals>0].astype(type(Nc))
            
            if len(new_keys_full)>0:
                return DirectSum(new_keys_full, new_values_full)
            else:
                return NullDiagram()
        
    def dimension_Nc(self,Nc=None):
    
        elements = [ky.dimension_Nc(Nc) for ky in self.keys()]
        multiplicities = list(self.values())
            
        return DimensionDirectSum(elements,multiplicities)
        
    def get_str(self):
    
        strin = ''
        elements = np.array(list(self.keys()))
        inds = elements.argsort()
        multiplicities = np.array(list(self.values()))[inds]
        elements = elements[inds]
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
        if len(strin)==0:
            strin = '0'
        
        return strin
        
    def get_cmdline_str(self):
    
        strin = ''
        elements = np.array(list(self.keys()))
        inds = elements.argsort()
        multiplicities = np.array(list(self.values()))[inds]
        elements = elements[inds]
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
    
    def _repr_latex_(self):            

        return f"$$ {self.get_str()} $$"
        
        
    def __add__(self,other):
        
        keys = []
        values = []
        
        if type(other)==DirectSum:
            keys = list(self.keys())+list(other.keys())
            values = list(self.values())+list(other.values())
            return DirectSum(keys,values)
        elif type(other) in [YoungDiagram,NullDiagram,Pair]:
            
            keys = list(self.keys())
            keys.append(other)
            
            values = list(self.values())
            values.append(other.weight)
            return DirectSum(keys,values)
        else:
            raise NotImplemented

        
    def __iadd__(self,other):
    
        added = self+other
        self = added
            
    def __radd__(self,other):
        return self+other
        
    def __sub__(self,other):
        
        keys = []
        values = []
        
        if type(other)==DirectSum:
            keys = list(self.keys())+list(other.keys())
            values = list(self.values())+list(-1*np.array(list(other.values())))
            return DirectSum(keys,values)
        elif type(other) in [YoungDiagram,NullDiagram,Pair]:
            
            keys = list(self.keys())
            keys.append(other)
            
            values = list(self.values())
            values.append(-1*other.weight)
            
            return DirectSum(keys,values)
        else:
            raise NotImplemented
            
        
        
    def __rsub__(self,other):
        return self-other
        
    def __isub__(self,other):
    
        added = self-other
        self = added
    
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
        
    def __eq__(self,other):
    
        try:
            if len((self-other).elements())==0:
                return True
            else:
                return False
        except Exception as e:
            return False
        
          




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
        
    def lowest_Nc(self):
        raise AttributeError("DimensionDirectSum only carries dimension info and cannot recover partitions or lowest Nc.")
        
    

