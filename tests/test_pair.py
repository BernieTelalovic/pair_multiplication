import unittest
from pair_multiplication import *
import numpy as np

class TestPair(unittest.TestCase):

    def setUp(self):
    
        self.test_partitions2 = [((),()),((),(1,)),((1,),(1,)), 
                                 ((),(2,1)), ((1,1),(2,)), ((2,1),(2,1)), 
                                 ((3,2),(3,2,1)), ((2,2),(2,2))
                               ]
        self.test_partitions1 = [((),()),((1,),(1,)),((),(2,1)),
                                      ((2,2),(2,2)), #((3,2),(3,2,1))
                                       ]
        self.pairs1 = [Pair(tup) for tup in self.test_partitions1]                  
        self.pairs2 = [Pair(tup) for tup in self.test_partitions2]
        self.pr1 = None
        self.pr2 = None
        
    def test_Nc_inheritence(self):

        pair = Pair((1,1),inherited_N0 = 5)
        
        self.assertEqual(pair.N0, 5, 
                        "Incorrect N0.")
        
        
    def test_evaluation_for_Nc3(self):

        pair = Pair((1,1))
        
        expected = YoungDiagram((2,1),Nc=3)
        
        self.assertEqual(pair.evaluate_for_Nc(Nc=3), expected, 
                        "Evaluating pair for Nc=3 did not produce correct diagram.")
                        
    def test_LR_compare_with_CA_for_minNc(self):
        yd_barred = Pair(((2,1),(1)))
        yd_unbarred = Pair(((3,2),(2,1)))
        
        barred_tensor_unbarred_CA = yd_barred*yd_unbarred
        minN0 = min(barred_tensor_unbarred_CA.lowest_Nc())
        
        barred_tensor_unbarred_LR = yd_barred.evaluate_for_Nc(minN0).LR_multiply(yd_unbarred.evaluate_for_Nc(minN0))
        
        
        self.assertEqual(barred_tensor_unbarred_CA.evaluate_for_Nc(minN0), barred_tensor_unbarred_LR.evaluate_for_Nc(minN0), 
                        "Multiplying pairs of diags with CA =/= LR for min Nc.")
                        
    def test_pair_multiple_commutative(self):
        pr1 = Pair(((1,1),(2)))
        pr2 = Pair(((1),(1)))
        
        self.assertEqual(pr1*pr2,pr2*pr1, 
                        "Multiplying pairs of diags does not commute.")
                        
                        
    def test_many_multiplications(self):
                
        for pr1 in self.pairs1:
            for pr2 in self.pairs2:
            
                self.pr1 = pr1
                self.pr2 = pr2
            
                with self.subTest():
                    mul = self.pr1*self.pr2
                    minN0 = 0
                    if len(list(mul.lowest_Nc())) > 0:
                        minN0 = min(list(mul.lowest_Nc()))
                                
                    mul_LR = (self.pr1.evaluate_for_Nc(minN0).LR_multiply(self.pr2.evaluate_for_Nc(minN0))).evaluate_for_Nc(minN0)
                                
                    self.assertEqual(mul.evaluate_for_Nc(minN0), mul_LR,  
                                         "Multiplying pairs of diags with CA =/= LR for min Nc for pairs: "+\
                                         str(self.pr1)+' and '+str(self.pr2)+' with Nc = '+str(minN0))
