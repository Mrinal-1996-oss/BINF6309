#!/usr/bin/env python3
#translate_APOE.py
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
read_sample = "/home/mrinalsubash.f/BINF6309/MSA/apoe_aa.fasta"
for record in SeqIO.parse(read_sample,"fasta"):
    rna = record.seq.transcribe()
    protein = rna.translate()
    my_record = SeqRecord(protein,record.id)
   # print(type(rna))
   # print(record)
   # print(type(record))
   # print(record.id)
   # print(type(record.seq))
   # break
   # print(type(protein))
    #SeqIO.write(protein, "apoe_aa.fasta", "fasta")
    print(my_record)
