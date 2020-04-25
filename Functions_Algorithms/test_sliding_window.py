#!/usr/bin/env python3
# test_basic_functions.py

import unittest
from sliding_window import sliding_window, gc_content

class TestSlidingWindow(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(sliding_window(1, ''), [], "Sliding window of empty string should be an empty list")
    
    def test_basic(self):
        self.assertEqual(sliding_window(1, 'ABC'), ['A', 'B', 'C'], "Sliding window of size 1 should just be a string's characters")
    
    def test_complex(self):
        self.assertEqual(sliding_window(5, "ATCGAGCTAGCGAGAGCTAGC"), ['ATCGA', 'TCGAG', 'CGAGC', 'GAGCT', 
        'AGCTA', 'GCTAG', 'CTAGC', 'TAGCG', 'AGCGA', 'GCGAG', 'CGAGA', 'GAGAG', 'AGAGC', 'GAGCT', 'AGCTA', 'GCTAG', 'CTAGC'],
        "A more complex sliding window test")
    
class TestGC(unittest.TestCase):
    def test_empty(self):
        self.assertRaises(ZeroDivisionError, gc_content, "")
    
    def test_allA(self):
        self.assertEqual(gc_content("AAAA"), 0, "GC content of a string of A's should be zero")
    
    def test_random(self):
        self.assertAlmostEqual(gc_content("ATTCGGCGCGCTAT"), 8/14, "More complex test")

if __name__ == '__main__':
    unittest.main()