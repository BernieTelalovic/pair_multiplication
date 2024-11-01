# direct_sum.py
import numpy as np
from .young_diagram import YoungDiagram
from .young_diagram import NullDiagram
from .young_diagram import Pair
from .utils import *
import itertools as it

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
        
        keys = self.keys()
        values = list(self.values())
        
        new_keys = np.array([ky.conjugate(Nc=Nc) for ky in keys])
        new_vals = np.array([values[ind]*new_keys[ind].multiplicity(Nc) 
                             for ind in range(len(new_keys))])
        
        return DirectSum(new_keys[new_vals>0], new_vals[new_vals>0].astype(type(Nc)))
    
    
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
    
    def set_N0(self,val):
        
        new_keys = [ky.set_N0(val) for ky in self.keys()]
        
        self = DirectSum(new_keys,list(self.values()))
