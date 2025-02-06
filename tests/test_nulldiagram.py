import unittest
from pair_multiplication import *


class TestNullDiagram(unittest.TestCase):

    def setUp(self):
        self.null = NullDiagram()
        
    def test_null(self):
    
        self.assertEqual(self.null, 0, 
                        "NullDiagram is not 0.")
                        
    def test_multiplication(self):
    
        yd = YoungDiagram((3,3,2,1))
        
        self.assertEqual(self.null*yd+yd*self.null, 0, 
                        "Multiplying by a NullDiagram does not return 0.")
                        
    def test_addition(self):
    
        yd = YoungDiagram((3,3,2,1))
        
        self.assertEqual(self.null+yd, yd, 
                        "Adding a NullDiagram to a YoungDiagram does not return the YoungDiagram.")
                        
        self.assertEqual(yd-self.null, yd, 
                        "Subtracting a NullDiagram from a YoungDiagram does not return the YoungDiagram.")
