import unittest
from pair_multiplication import *

class TestDirectSum(unittest.TestCase):

    def setUp(self):
        self.instances = [
                          YoungDiagram((2,1))
                         ]
        self.expected_partitions_Nc3 = []
        self.expected_dimensions_Nc3 = []
        
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
