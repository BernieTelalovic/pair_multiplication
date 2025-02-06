import unittest
from pair_multiplication import *


class TestPrinting(unittest.TestCase):

    def setUp(self):
        self.null = NullDiagram()
        self.diag = YoungDiagram((3,2,1))
        self.diag_barred = YoungDiagram((3,2),barred = True)
        self.pair = Pair(((1,1),(2,1)))
        self.ds = self.diag+self.pair+self.diag_barred
        
    def test_null_cmd_print(self):
    
        self.assertEqual(self.null.get_cmdline_str(), '0', 
                        "Command line string incorrect for NullDiagram.")
                        
    def test_null_latex_print(self):
    
        self.assertEqual(self.null._repr_latex_(), '$$ 0 $$', 
                        "LaTeX string incorrect for NullDiagram.")
                        
    def test_yd_cmd_print(self):
    
        self.assertEqual(self.diag.get_cmdline_str(), '(3, 2, 1)', 
                        "Command line string incorrect for YoungDiagram.")
                        
    def test_yd_latex_print(self):
    
        self.assertEqual(self.diag._repr_latex_(), '$$ (3, 2, 1) $$', 
                        "LaTeX string incorrect for YoungDiagram.") 
                        
    def test_ydbr_cmd_print(self):
    
        self.assertEqual(self.diag_barred.get_cmdline_str(), '(3, 2)_', 
                        "Command line string incorrect for barred YoungDiagram.")
                        
    def test_ydbr_latex_print(self):
    
        self.assertEqual(self.diag_barred._repr_latex_(), '$$ \\overline{(3, 2)} $$', 
                        "LaTeX string incorrect for barred YoungDiagram.") 
                        
    def test_pair_cmd_print(self):
    
        self.assertEqual(self.pair.get_cmdline_str(), '[4]((1, 1),(2, 1))', 
                        "Command line string incorrect for Pair.")
                        
    def test_pair_latex_print(self):
    
        self.assertEqual(self.pair._repr_latex_(), '$$ 1_{4} \\left(\\overline{(1, 1)},(2, 1)\\right) $$', 
                        "LaTeX string incorrect for Pair.") 
                        
                        
    def test_directsum_cmd_print(self):
    
        self.assertEqual(self.ds.get_cmdline_str(), '1(3, 2)_+1(3, 2, 1)+1[4]((1, 1),(2, 1))', 
                        "Command line string incorrect for DirectSum.")
                        
    def test_directsum_latex_print(self):
    
        self.assertEqual(self.ds._repr_latex_(), '\\[\\begin{array}{c}1_{2}\\,\\overline{(3, 2)}\\oplus\n\n1_{3}\\,(3, 2, 1)\\oplus\n\n1_{4} \\left(\\overline{(1, 1)},(2, 1)\\right)\\end{array}\\]', 
                        "LaTeX string incorrect for DirectSum.") 
