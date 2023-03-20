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
import shutil
from Bio.SeqRecord import SeqRecord



from utils import change_seq


# Store all the multiCDS for each species
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

os.mkdir("Modified_data")
os.chdir("Modified_data")

# keys = 100, 200, 400
for condition in sequences.keys():

    os.mkdir(f'cond_{condition}')

    print(f'Condition : {condition}')

    Scer_records = list(SeqIO.parse("../orthoData/Scer_NCBI_CDS.pep", "fasta"))

    for id in sequences[condition]:

        print(f'Seq ID : {id}')

        records_without_target_seq = []

        for record in Scer_records:

            ''' Create tmp records object without sequence to be modified '''

            if record.id != id:

                records_without_target_seq.append(record)

            elif record.id == id:

                temp_CDS = SeqRecord(seq = record.seq,
                                     id = record.id, 
                                     description = "")


        buffer = records_without_target_seq


        os.mkdir(f'cond_{condition}/seq_{id}')
        
        for action, bool in zip(("remove","append"),(False,True)):

            os.mkdir(f'cond_{condition}/seq_{id}/{action}')

            for length in (10,20,30,40,50):

                os.mkdir(f'cond_{condition}/seq_{id}/{action}/{length}')
                
                tmp_records = records_without_target_seq.copy() # Make tmp copy of object


                changed_tmp_record = change_seq(temp_CDS, percent = length, lengthen = bool)
                tmp_records.append(changed_tmp_record)

                with open(f"cond_{condition}/seq_{id}/{action}/{length}/Scer_NCBI_CDS.pep","w") as f:
                    SeqIO.write(tmp_records, f, "fasta")
                f.close()

                shutil.copyfile("../orthoData/Sarb_CDS.pep", f"cond_{condition}/seq_{id}/{action}/{length}/Sarb_CDS.pep")
                shutil.copyfile("../orthoData/Sbay_CDS.pep", f"cond_{condition}/seq_{id}/{action}/{length}/Sbay_CDS.pep")
                shutil.copyfile("../orthoData/Skud_CDS.pep", f"cond_{condition}/seq_{id}/{action}/{length}/Skud_CDS.pep")
                shutil.copyfile("../orthoData/Smik_CDS.pep", f"cond_{condition}/seq_{id}/{action}/{length}/Smik_CDS.pep")
                shutil.copyfile("../orthoData/Spar_NCBI_CDS.pep", f"cond_{condition}/seq_{id}/{action}/{length}/Spar_NCBI_CDS.pep")
                        


                

                



