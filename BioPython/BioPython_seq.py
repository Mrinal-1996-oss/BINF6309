#!/usr/bin/env python3
#BioPython_seq.py

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

#Creates a SeqRecord object with 'seq','id', 'description and alphabet parameters
seq = Seq("aaaatgggggggggggccccgtt")
SeqRecord = SeqRecord(seq, id = "#12345", description = "example1")
#writes the object to a sequence file in GenBank format
SeqIO.write(SeqRecord,"BioPython_seq.gb","genbank")
