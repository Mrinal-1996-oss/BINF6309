#!/usr/bin/env python3
#basic_functions.py

#multiply.py
def multiply(a,b):
    return a*b

x=multiply(5,10)
print("The value of x is {}".format(x))

#hello_name.py
def hello_name(name="you"):
        print("Hello, {}!".format(name))
hello_name()

#less_than_ten.py
def less_than_ten(list):
    number_list = []
    for i in list:
        if i < 10:
            number_list.append(i)
    return number_list
number_list1 = less_than_ten([1, 5, 81, 10, 8, 2, 102])
print(number_list1)
