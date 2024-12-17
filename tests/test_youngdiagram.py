import unittest
from pair_multiplication import *
#from .utils_for_testing import *


class TestYoungDiagram(unittest.TestCase):
    
    def setUp(self):
        self.instances = [
                          YoungDiagram((2,1))
                         ]
        self.expected_partitions_Nc3 = []
        self.expected_dimensions_Nc3 = []
        
    def test_yd_correct_construction_no_Nc(self):
        # tests that we can construct an admissible young diagram
        partition = (3,3,2,1)
        yd = YoungDiagram((3,3,2,1))
        correct_info_tuple = (partition, len(partition), None)
        info_tuple = (yd.partition,yd.N0,yd.Nc)
        
        self.assertEqual(info_tuple, correct_info_tuple, 
                        "Construction parameter do not match.")
    
    def test_yd_incorrect_construction(self):
    
        with self.assertRaises(ValueError) as context:
            yd = YoungDiagram((1,1,2))
        self.assertEqual(str(context.exception), 'Not a young diagram.', 
                        "Incorrect error message during inadmissible construction.")
                        
    def test_null_diagram(self):
    
        null = NullDiagram()
        
        self.assertEqual(null.weight, 0, 
                        "Incorrect weight for null diagram.")
        self.assertEqual(null.multiplicity(100), 0, 
                        "Incorrect multiplicity for null diagram.")
        self.assertEqual(null.dimension_Nc(100), 0, 
                        "Incorrect dimension for null diagram.")

    def test_partition_given_Nc3(self):
        """Test that multiple instances give expected results."""
        yd_nc_indep = YoungDiagram((3,2,1))
        yd_nc3 = YoungDiagram((3,2,1),Nc=3)
        
        expected = YoungDiagram((2,1),Nc=3)
        
        self.assertEqual(yd_nc_indep.evaluate_for_Nc(Nc=3), expected, 
                        "Diagram construction failed for evaluated Nc.")
        self.assertEqual(yd_nc3, expected, 
                        "Diagram construction failed for given Nc.")
                        
    def test_reduction_for_too_small_Nc(self):
    
        yd_nc_indep = YoungDiagram((3,2,1))
        yd_nc1 = YoungDiagram((3,2,1),Nc=1)
        
        expected = NullDiagram()
        
        self.assertEqual(yd_nc_indep.evaluate_for_Nc(Nc=1), expected, 
                        "Diagram construction failed for evaluated Nc.")
        self.assertEqual(yd_nc1, expected, 
                        "Diagram construction failed for given Nc.")
                        
                        
    def test_dimension_is_correct(self):
    
        yd = YoungDiagram((2,1))
        
        self.assertEqual(yd.dimension_Nc(3), 8, 
                        "Incorrect dimension for Nc=3.")
                    
                        
    def test_multiply_with_null_diagram(self):
    
        diag = YoungDiagram((2,1))
        null = NullDiagram()
        
        self.assertEqual(diag*null, null, 
                        "Multiplying by null diagram does not return null diagram.")
                        
                        
    def test_multiply_two_unbarred_diags(self):
    
        yd1 = YoungDiagram((2))
        yd2 = YoungDiagram((1,1))
        expected = DirectSum([YoungDiagram((3,1)),YoungDiagram((2,1,1))],[1,1])
        
        self.assertEqual(yd1*yd2, expected, 
                        "Multiplying by unbarred diagrams is incorrect.")
                        
    def test_multiply_two_barred_diags(self):
    
        yd1 = YoungDiagram((2),barred=True)
        yd2 = YoungDiagram((1,1),barred=True)
        expected = DirectSum([YoungDiagram((3,1),barred=True),YoungDiagram((2,1,1),barred=True)],[1,1])
        
        self.assertEqual(yd1*yd2, expected, 
                        "Multiplying by barred diagrams is incorrect.")
                        
    def test_multiply_barred_unbarred_diags(self):
    
        yd1 = YoungDiagram((2))
        yd2 = YoungDiagram((1,1),barred=True)
        expected = DirectSum([Pair(((1),(1)),inherited_N0=2),Pair(((1,1),(2)),inherited_N0=3)],[1,1])
        
        self.assertEqual(yd1*yd2, expected, 
                        "Multiplying by barred and unbarred diagrams is incorrect.")
         
    def test_LR_compare_with_CA_unbarred_unbarred(self):
        yd_barred = YoungDiagram((2,1),barred = False)
        yd_unbarred = YoungDiagram((2,1),barred = False)
        
        barred_tensor_unbarred_LR = yd_barred.LR_multiply(yd_unbarred)
        barred_tensor_unbarred_CA = yd_barred*yd_unbarred
        
        self.assertEqual(barred_tensor_unbarred_CA, barred_tensor_unbarred_LR, 
                        "Multiplying unbarred diags with CA =/= LR.")
         
    def test_LR_compare_with_CA_barred_barred(self):
        yd_barred = YoungDiagram((2,1),barred = True)
        yd_unbarred = YoungDiagram((2,1),barred = True)
        
        barred_tensor_unbarred_LR = yd_barred.LR_multiply(yd_unbarred)
        barred_tensor_unbarred_CA = yd_barred*yd_unbarred
        
        self.assertEqual(barred_tensor_unbarred_CA, barred_tensor_unbarred_LR, 
                        "Multiplying barred diags with CA =/= LR.")
                        
    def test_LR_compare_with_CA_barred_unbarred(self):
        yd_barred = YoungDiagram((2,1),barred = True)
        yd_unbarred = YoungDiagram((2,1))
        
        barred_tensor_unbarred_LR = yd_barred.LR_multiply(yd_unbarred)
        barred_tensor_unbarred_CA = yd_barred*yd_unbarred
        
        self.assertEqual(barred_tensor_unbarred_CA, barred_tensor_unbarred_LR, 
                        "Multiplying barred and unbarred diags with CA =/= LR.")
                        
                        
                        
                        
    
        
                        
                        
                        
