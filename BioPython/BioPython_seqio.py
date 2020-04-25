#!/usr/bin/env python3
# BioPython_seqio.py

#Reads the multi-sequence FASTA file with SeqIO

import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
if __name__ == "__main__":
    arg_count = len(sys.argv) -1 
    if arg_count != 2:
        #script should take two arguments from the command line:
        #'the name of the original FASTA file' and
        #' the desired name of the new FASTA file'
        raise Exception("This script requires exactly 2 arguments")
    file_name_1 = sys.argv[1]
    file_name_2 = sys.argv[2]

   #Outputs a new FASTA file with contents are the reverse complements of the sequences from the original FASTA file.
   Records = [rec.reverse_complement(id - "rc_" + rec.id, description ='')\
           for rec in SeqIO.parse(file_name_1, "fasta")]
   SeqIO.write(records, file_name_2, "fasta")
