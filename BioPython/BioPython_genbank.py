#!/usr/bin/env python3
#BioPython_genbank.py

#Creates a list with 'sequence retrieved from the Genbank and 'a sequence retrieved from Genbank by accession ID

from Bio import Entrez
from Bio import SeqIO

Entrez.email = "mrinalsubash.f@husky.neu.edu"
#create a list for retrieved sequences
Sequences = []
#Entrez.efetch returns the SeqIO object
with Entrez.efetch (
     db = "nucleotide", rettype = "gb",retmode = "text", id = "515056"
     ) as handle:
          seq_record_1 = SeqIO.read(handle,"gb")
with Entrez.efetch(
     db = "nucleotide",rettype = "gb", retmode = "text", id = "J01673.1"
     ) as handle:
          seq_record_2= SeqIO.read(handle,"gb")
Sequences = [seq_record_1, seq_record_2]
#Prints out the sequences from the list
print(Sequences)

# the .features attributes contains a list of objects that have .type, .location and .strand attributes.
Seq_features = [seq_record_1.features, seq_record_2.features]
#Prints out the type, location and strand of each feature.
print(Seq_features)

