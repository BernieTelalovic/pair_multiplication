import unittest
from pair_multiplication import *

class TestDirectSum(unittest.TestCase):

    def setUp(self):
        self.instances = [
                          YoungDiagram((2,1))
                         ]
        self.expected_partitions_Nc3 = []
        self.expected_dimensions_Nc3 = []
        
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
        
        
        self.assertEqual(barred_tensor_unbarred_CA.evaluate_for_Nc(minN0), barred_tensor_unbarred_LR, 
                        "Multiplying pairs of diags with CA =/= LR for min Nc.")
