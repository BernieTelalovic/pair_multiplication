import unittest
from pair_multiplication import *

class TestDirectSum(unittest.TestCase):

    def setUp(self):
        self.instances = [
                          YoungDiagram((2,1))
                         ]
        self.expected_partitions_Nc3 = []
        self.expected_dimensions_Nc3 = []
        
    def test_incorrect_constructor(self):

        with self.assertRaises(ValueError) as context:
            to_simplify = DirectSum([YoungDiagram((3,1),barred=True),
                                 YoungDiagram((2,1,1),barred=True),
                                 YoungDiagram((2,1,1),barred=True)],
                                [1])
        
        self.assertEqual(str(context.exception), 'List of multiplicities must have equal length to list of diagrams/pairs.', 
                        "Incorrect error message during inadmissible construction.")
        
    def test_correct_constructor_simplification(self):

        to_simplify = DirectSum([YoungDiagram((3,1),barred=True),
                                 YoungDiagram((2,1,1),barred=True),
                                 YoungDiagram((2,1,1),barred=True)],
                                [1,1,1])
        
        expected = DirectSum([YoungDiagram((3,1),barred=True),YoungDiagram((2,1,1),barred=True)],
                             [1,2])
        
        self.assertEqual(to_simplify, expected, 
                        "Simplification during construction is incorrect.")
                        
                        
    def test_addition_by_null(self):
    
        expected = DirectSum([YoungDiagram((3,1),barred=True),YoungDiagram((2,1,1),barred=True)],
                             [1,2])
                             
        self.assertEqual(expected+NullDiagram(), expected, 
                        "Adding Null diagram should not change a direct sum.")
                        
    def test_addition_by_diagram(self):
    
        initial = DirectSum([YoungDiagram((3,1),barred=True),YoungDiagram((2,1,1),barred=True)],
                             [1,1])
        expected = DirectSum([YoungDiagram((3,1),barred=True),YoungDiagram((2,1,1),barred=True)],
                             [1,2])
                             
        self.assertEqual(initial+YoungDiagram((2,1,1),barred=True), expected, 
                        "Diagram cannot be added.")
                        
    def test_subtraction_by_diagram(self):
    
        initial = DirectSum([YoungDiagram((3,1),barred=True),YoungDiagram((2,1,1),barred=True)],
                             [1,2])
        expected = DirectSum([YoungDiagram((3,1),barred=True),YoungDiagram((2,1,1),barred=True)],
                             [1,1])
                             
        self.assertEqual(initial-YoungDiagram((2,1,1),barred=True), expected, 
                        "Diagram cannot be added.")
                        
                        
    def test_direct_sum_multiplicities(self):

        ds = YoungDiagram((2,1),Nc=3)*YoungDiagram((2,1),Nc=3)
        expected_multiplicities = sorted([2, 1, 1,1, 1])
        expected_elements = sorted([8, 1, 10,10, 27])
        
        self.assertEqual(sorted(ds.multiplicities), expected_multiplicities, 
                        "Incorrect dimensional multiplicities.")
        self.assertEqual(sorted(ds.dimension_Nc()), expected_elements, 
                        "Incorrect dimensional elements.")
                        
    def test_direct_sum_dimension(self):
    
        yd = YoungDiagram((2,1),Nc=3)
        ds = yd*yd
        expected_sum = 64
        
        self.assertEqual(ds.sum_dimensions(), expected_sum, 
                        "Incorrect dimensional sum.")
                        
                        
    def test_evaluating_under_Nc3(self):
    
        yd = YoungDiagram((2,1))
        ds = yd*yd
        
        expected_ds = DirectSum([YoungDiagram((),Nc=3), YoungDiagram((3),Nc=3), YoungDiagram((3,3),Nc=3),
                                 YoungDiagram((4,2),Nc=3),YoungDiagram((2,1),Nc=3)],
                                [1,1,1,1,2])
        self.assertEqual(ds.evaluate_for_Nc(3), expected_ds, 
                        "Incorrect direct sum evaluation.")
