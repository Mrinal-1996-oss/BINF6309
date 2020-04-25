#!/usr/bin/env python3
# sliding_window.py


import re 
import sys


def sliding_window(k, dna):
    ''' Returns a list of all k-mers in the given string
    '''
    kmers = []
    end = len(dna) - int(k) + 1
    for start in range(0, end):
        kmers.append(dna[start:start + int(k)])

    return kmers
def gc_content(dna):
    
    dna=dna.lower()
    gc = 0
    for nucleotide in dna:
        if nucleotide in ['g' , 'c']:
            gc += 1
    return float(gc/len(dna))

#def gc_content(string):
#    ''' Returns [0, 1], the fraction of GCs in the given string'''
    # For consistency, make the sequence lowercase
#    string = string.lower()

    # Count the number of g's and c's
#    gc = 0
#    for nucleotide in string:
#        if nucleotide in ['g', 'c']:
#           gc += 1

#    return gc/len(string)

if __name__ == "__main__":
# Check to make sure there are at least two arguments
    arg_count = len(sys.argv) - 1
    if arg_count < 2:
       raise Exception("This script requires 2 arguments: a kmer size and then a string")

    k = int(sys.argv[1])
    dna = sys.argv[2]
    kmers = sliding_window(k,dna)
    kmer_gc={}
   # Calculate the GC content
    #result = gc_content(dna)
    for i in kmers:
        kmer_gc[i] = gc_content(i)

    for key,value in kmer_gc.items():
        # Print the result, rounding GC content to 2 decimal places
       # print(key + '\t' +str(round(value,2)))#".format(kmers[i],gc_content(kmers[i])))
       print("{0}\t{1:.2f}".format(key,value))
       #print ("This site is {0:f}% securely {1}!!". 
                                  # format(100, "encrypted"))

    for i in range(len(kmer_gc)):
        print("{}\t{:.2f}".format(i,gc_content())
