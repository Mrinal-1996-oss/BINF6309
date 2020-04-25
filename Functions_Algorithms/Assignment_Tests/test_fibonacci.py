#!/usr/bin/env python3
# test_basic_functions.py

import unittest
from fibonacci import population

class TestFibonacci(unittest.TestCase):
    def test_1_lots(self):
        self.assertEqual(population(1,1000), 1, "At day 1, population should be 1, regardless of breeding")

    def test_5_3(self):
        self.assertEqual(population(5, 3), 19, "n=5 and k=3 should return 19")
    
    def test_29_4(self):
        self.assertEqual(population(29, 4), 170361678269, "Slightly more complex test")


if __name__ == '__main__':
    unittest.main()