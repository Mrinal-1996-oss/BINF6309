#!/usr/bin/env python3
#hamming.py
import sys

def hamming(char_a,char_b):
    '''Returns the hamming distance of two sequence of equal length using list loops and zip'''
    initial= 0
    for sub1,sub2 in zip(char_a,char_b):
        if sub1 != sub2:
           initial +=1
    return initial

if __name__ =="__main__":
    arg_count=len(sys.argv) -1
    if arg_count < 2 :
        raise Exception("This script requires two arguments of equal length")

    a= sys.argv[1]
    b= sys.argv[2]
    if len(a) != len(b):
        raise ValueError("Sequence need to be the same length")
    result =hamming(a,b)
    print('{}\t{}\t{}'.format(a,b,result))
