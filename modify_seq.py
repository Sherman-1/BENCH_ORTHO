#!/bin/python3


import os
import glob
from Bio import SeqIO
import re


# Scer >NC_001133.9_-_14193-14435_2_nc_intergenic HCA score 0.31 for complete iORF
iORF = "IVSGKGKVIMQWNLRYSHVQISATKSVVGTFRYQNFLSTTREHHYAKQIVAHTSDRKKTVPKTMTYEETSIRIFIINSMAIVSGKGKVIMQWNLRYSHVQISATKSVVGTFRYQNFLSTTREHHYAKQIVAHTSDRKKTVPKTMTYEETSIRIFIINSMA"


files = glob.glob("rna*")


records = (r for r in SeqIO.parse("orthoData/Scer_NCBI_CDS.pep", "fasta"))
sequences = {}
for f in files:

    with open(f,"r") as id_list:

        ids = id_list.readlines()
        ids = [sub.replace('\n', '') for sub in ids]
        ids = [sub.replace('>', '') for sub in ids]
        length = re.findall(r"\d+", f)[0]
        sequences[length] = {}
        for id in ids:

            sequences[length][id] = str()
    
print(sequences)

for record in records:

    for length in sequences: 

        for id in sequences[length]:

            if id in record.id:

                # Found the sequence with the id "my_super_id"
                temp_seq = record.seq.replace("*","")
                sequences[length][id] = temp_seq
        

            
print(sequences)
       
        

        

