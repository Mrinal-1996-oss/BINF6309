#!/usr/bin/env python3
#sliding_window_fasta.py

import re
import sys

def sliding_window(k, dna):
    
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

if __name__ == "__main__":
    arg_count= len(sys.argv) -1
    if arg_count < 2 :
        raise Exception("This script requires two arguments")
    k= int(sys.argv[1])
    filehandle = sys.argv[2]
    sequences = ''
    with open(filehandle, 'r') as dna:
        sequence= ''
        for line in dna:
            line = line.rstrip()
            if len(line) <1:
                continue
            elif line[0] =='>':
                 print(line)
                 sequence = ''
            else:
                 sequence += line



    kmers= sliding_window(k, sequence)
    for kmer in kmers:
        print("{}\t{:.2f}".format(kmer,gc_content(kmer)))
