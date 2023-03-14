#!/bin/python3

'''
This program takes care of generating data for benchmarking OrthoFINDER on elongated or depleted 
sequences from already known orthougroups. The main goal is to assess wether or not OF is able 
to successfully regroup sequences even in the presence of an elongation by an iORF on one of the sequences

We also test for depletion on one of the sequences.

'''
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
import re
from utils import change_seq


# Scer >NC_001133.9_-_14193-14435_2_nc_intergenic HCA score 0.31 for complete iORF
iORF = "IVSGKGKVIMQWNLRYSHVQISATKSVVGTFRYQNFLSTTREHHYAKQIVAHTSDRKKTVPKTMTYEETSIRIFIINSMAIVSGKGKVIMQWNLRYSHVQISATKSVVGTFRYQNFLSTTREHHYAKQIVAHTSDRKKTVPKTMTYEETSIRIFIINSMA"


# Store all the multiCDS for each species
Scer_records = SeqIO.parse("orthoData/Scer_NCBI_CDS.pep", "fasta")
Sarb_records = SeqIO.parse("orthoData/Sarb_CDS.pep", "fasta")
Sbay_records = SeqIO.parse("orthoData/Sbay_CDS.pep", "fasta")
Skud_records = SeqIO.parse("orthoData/Skud_CDS.pep", "fasta")
Spar_records = SeqIO.parse("orthoData/Spar_NCBI_CDS.pep", "fasta")
Smik_records = SeqIO.parse("orthoData/Smik_CDS.pep","fasta")




files = glob.glob("rna*")
sequences = {}

for f in files:

    with open(f,"r") as id_list:

        ids = id_list.readlines()
        ids = [sub.replace('\n', '') for sub in ids]
        ids = [sub.replace('>', '') for sub in ids]
        length = re.findall(r"\d+", f)[0]
        sequences[int(length)] = []

        for id in ids:

            sequences[int(length)].append(id)
    
print(sequences)





# keys = 100, 200, 400
for condition in sequences.keys():

    os.mkdir(f'cond_{condition}')
    os.chdir(f'cond_{condition}')

    # We iterate over the CDS randomly selected by get_seq_from_OG.sh
    for id in sequences[condition]:

        # For each ID of CDSs of interest, we create a folder to store their modified versions
        os.mkdir(f'seq_{id}')
        os.chdir(f'seq_{id}')

    
        # We isolate the CDS of interest so it can be modified 
        temp_CDS = [rec for rec in temp_scer_records if rec.id == id][0]

        # Each CDS is going to be diminished ... 
        os.mkdir('remove')

        # ... by 10% of its length, then 20%, then 30 ...
        os.mkdir('remove/len_10')

        # We remove the original CDS from a temporary Scer multiCDS fasta
        temp_scer_records = Scer_records
        temp_scer_records = [rec for rec in temp_scer_records if rec.id != id]

        changed_CDS = change_seq(temp_CDS)

        



        os.mkdir('remove/len_20')
        os.mkdir('remove/len_30')
        os.mkdir('remove/len_40')
        os.mkdir('remove/len_50')

        os.mkdir('append')
        os.mkdir('append/len_10')
        os.mkdir('append/len_20')
        os.mkdir('append/len_30')
        os.mkdir('append/len_40')
        os.mkdir('append/len_50')
