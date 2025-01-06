import unittest
from pair_multiplication import *
from pair_multiplication.utils import *
import numpy as np

class TestUtils(unittest.TestCase):

    def setUp(self):
        self.info = None
        
    def test_Nc_positive(self):
        
        self.assertEqual(Nc_positive(-1), 0, 
                        "Incorrect Nc.")
        self.assertEqual(Nc_positive(1), 1, 
                        "Incorrect Nc.")
    
    def test_partition_tuplify(self):
    
        self.assertEqual(partition_tuplify(0), (), 
                        "Incorrect tuple.")
                        
        self.assertEqual(partition_tuplify((2)), (2,), 
                        "Incorrect tuple.")
                        
    def test_is_monotonic_increasing(self):
    
        word_increasing = [1,1,2,2]
        word_not_increasing = [1,1,2,1]
        no_word = []
        
        self.assertTrue(is_monotonic_increasing(word_increasing), 
                        "Monotonic increasing word marked as not.")
        self.assertFalse(is_monotonic_increasing(word_not_increasing), 
                        "Not monotonic increasing word marked as is.")
        self.assertTrue(is_monotonic_increasing(no_word), 
                        "Empty word marked as not.")
                        
    def test_is_serial_increasing(self):
    
        word_increasing = [1,2]
        no_word = []
        
        self.assertTrue(is_serial_increasing(word_increasing), 
                        "Serial increasing word marked as not.")
        self.assertTrue(is_serial_increasing(no_word), 
                        "Empty word marked as not serially increasing.")
                        
    def test_is_serial_increasing_from_its_min(self):
    
        word_increasing = [2,3]
        no_word = []
        
        self.assertTrue(is_serial_increasing_from_its_min(word_increasing), 
                        "Serial increasing word marked as not.")

        self.assertTrue(is_serial_increasing_from_its_min(no_word), 
                        "Empty word marked as not serially increasing.")
                        
    def test_is_serial_decreasing(self):
    
        word_not_decreasing = [1,2]
        no_word = []

        self.assertFalse(is_serial_decreasing(word_not_decreasing), 
                        "Not serial decreasing word marked as is.")
        self.assertTrue(is_serial_decreasing(no_word), 
                        "Empty word marked as not serially decreasing.")
                        
                        
    def test_is_serial_decreasing_from_its_min(self):
    
        word_not_decreasing = [2,3]
        no_word = []
        
        self.assertFalse(is_serial_decreasing_from_its_min(word_not_decreasing), 
                        "Not serial decreasing word marked as is.")
        self.assertTrue(is_serial_decreasing_from_its_min(no_word), 
                        "Empty word marked as not serially decreasing.")
                        
    def test_superpose_edge_cases(self):
    
        wordA_0 = np.array([0j])
        wordB_0 = np.array([0j])
        
        dic,lis,cont = superpose(0,wordA_0,wordB_0)
        self.assertEqual(lis[0], ((),()), 
                        "Incorrect pair tuple returned.")
        
        wordA = np.array([1+1j])
        dic,lis,cont = superpose(0,wordA,wordB_0)
        self.assertEqual(cont[0], '[]', 
                        "Incorrect pair tuple returned.")
        dic,lis,cont = superpose(0,wordB_0,wordA)
        self.assertEqual(cont[0], '[]', 
                        "Incorrect pair tuple returned.")
                        
                        
    def test_get_entries_to_add_empty_list(self):
    
        entries = get_entries_to_add(np.array([]),1)
        
        self.assertEqual(len(entries), 0, 
                        "Incorrect cendidates returned.")
                        
                        
    
