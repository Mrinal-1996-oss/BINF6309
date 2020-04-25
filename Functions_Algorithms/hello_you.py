#!/usr/bin/env python3
#hello_you.py

import sys
def hello_name(name="you"):
    print("Hello, {}!".format(name))

if __name__=='__main__':
    arg_count = len(sys.argv) -1
    if arg_count > 0:
        hello_name(sys.argv[1])
    else:
        hello_name()
          
