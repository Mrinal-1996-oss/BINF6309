#!/usr/bin/env python3
# test_basic_functions.py

import unittest
from basic_functions import multiply, hello_name, less_than_ten

class TestMultiply(unittest.TestCase):
    def test_zero_zero(self):
        self.assertEqual(multiply(0,0), 0, "0*0 = 0")
    
    def test_n_zero(self):
        self.assertEqual(multiply(37279, 0), 0, "37279*0 = 0")

    def test_n_1(self):
        self.assertEqual(multiply(28372, 1), 28372, "28372*1 = 28372")
    
    def test_swap(self):
        self.assertEqual(multiply(29374, 29), multiply(29, 29374), "xy=yx")

    def test_some_numbers(self):
        self.assertEqual(multiply(250, 137), 34250, "250*137 = 34250")


class TestHelloName(unittest.TestCase):
    def test_no_args(self):
        self.assertEqual(hello_name(), "Hello, you!", "If no args are passed, hello_name() should return 'Hello, you!'")
    
    def test_basic_name(self):
        self.assertEqual(hello_name("Chesley"), "Hello, Chesley!", "Basic hello_name functionality.")
    
class TestLessThanTen(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(less_than_ten([]), [], "Empty list should also return an empty list")
    
    def test_ten(self):
        self.assertEqual(less_than_ten([10]), [], "Should keep numbers less than, but not equal to, ten")

    def test_greater_ten(self):
        self.assertEqual(less_than_ten([20, 30, 40]), [], "List with only numbers greater than ten should return empty list")

    def test_less_ten(self):
        self.assertEqual(less_than_ten([1,2,3]), [1,2,3], "List with only numbers greater than ten should be unchanged")

    def test_mixed_list(self):
        self.assertEqual(less_than_ten([1,10,2, 20, 3, 30]), [1, 2, 3], "Should preserve order, and remove numbers less than ten")

if __name__ == '__main__':
    unittest.main()