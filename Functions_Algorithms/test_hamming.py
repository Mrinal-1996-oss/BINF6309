#!/usr/bin/env python3
# test_basic_functions.py

import unittest
from hamming import hamming

class TestHamming(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(hamming("", ""), 0, "Hamming distance of two empty strings should be zero")
    
    def test_same(self):
        self.assertEqual(hamming("ATCG", "ATCG"), 0, "Hamming distance of two identical strings should be zero")
    
    def test_different(self):
        self.assertEqual(hamming("GCTA", "ATCG"), 4, "Hamming distance of two completely different strings")

    def test_complex(self):
        self.assertEqual(hamming("ATCGATCGGATC", "AGCGATCGGAGC"), 2, "A more complex example")
    
if __name__ == '__main__':
    unittest.main()
