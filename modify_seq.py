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

    
        

        # Each CDS is going to be diminished ... 
        os.mkdir('remove')

         ######################

        # 10

        ######################


        # ... by 10% of its length, then 20%, then 30 ...
        os.mkdir('remove/len_10')

        # For that, we store a temporary version of the Scer multiCDS
        temp_scer_records = Scer_records

        # We isolate the CDS of interest so it can be modified 
        temp_CDS = [rec for rec in temp_scer_records if rec.id == id][0]

        # We remove the original CDS of interest from the temporary Scer multiCDS
        temp_scer_records = [rec for rec in temp_scer_records if rec.id != id]

        # We create a novel version of the CDS of interest ( here : trunkated by 10% )
        changed_CDS = change_seq(temp_CDS, percent = 10, lengthen = False)  

        temp_scer_records = list(temp_scer_records)
        temp_scer_records.append(changed_CDS)

        with open("remove/len_10/Scer_NCBI_CDS.pep","w") as f:
            SeqIO.write(temp_scer_records, f, "fasta")
        f.close()

        shutil.copyfile("../../orthoData/Sarb_CDS.pep", "remove/len_10/Sarb_CDS.pep")
        shutil.copyfile("orthoData/Sbay_CDS.pep", "remove/len_10/Sbay_CDS.pep")
        shutil.copyfile("orthoData/Skud_CDS.pep", "remove/len_10/Skud_CDS.pep")
        shutil.copyfile("orthoData/Smik_CDS.pep", "remove/len_10/Smik_CDS.pep")
        shutil.copyfile("orthoData/Spar_NCBI_CDS.pep", "remove/len_10/Spar_NCBI_CDS.pep")
        
        ######################

        # 20

        ######################

        os.mkdir('remove/len_20')

        temp_scer_records = Scer_records
        temp_scer_records = [rec for rec in temp_scer_records if rec.id != id]
        changed_CDS = change_seq(temp_CDS, percent = 20, lengthen = False)  
        temp_scer_records = list(temp_scer_records)
        temp_scer_records.append(changed_CDS)

        with open("remove/len_20/Scer_NCBI_CDS.pep","w") as f:
            SeqIO.write(temp_scer_records, f, "fasta")
        f.close()

        shutil.copyfile("orthoData/Sarb_CDS.pep", "remove/len_20/Sarb_CDS.pep")
        shutil.copyfile("orthoData/Sbay_CDS.pep", "remove/len_20/Sbay_CDS.pep")
        shutil.copyfile("orthoData/Skud_CDS.pep", "remove/len_20/Skud_CDS.pep")
        shutil.copyfile("orthoData/Smik_CDS.pep", "remove/len_20/Smik_CDS.pep")
        shutil.copyfile("orthoData/Spar_NCBI_CDS.pep", "remove/len_20/Spar_NCBI_CDS.pep")

             
        ######################

        # 30 

        ######################

        os.mkdir('remove/len_30')



        temp_scer_records = Scer_records

        # We remove the original CDS of interest from the temporary Scer multiCDS
        temp_scer_records = [rec for rec in temp_scer_records if rec.id != id]

        # We create a novel version of the CDS of interest ( here : trunkated by 30% )
        changed_CDS = change_seq(temp_CDS, percent = 30, lengthen = False)  

        temp_scer_records = list(temp_scer_records)
        temp_scer_records.append(changed_CDS)

        with open("remove/len_30/Scer_NCBI_CDS.pep","w") as f:
            SeqIO.write(temp_scer_records, f, "fasta")
        f.close()

        shutil.copyfile("orthoData/Sarb_CDS.pep", "remove/len_30/Sarb_CDS.pep")
        shutil.copyfile("orthoData/Sbay_CDS.pep", "remove/len_30/Sbay_CDS.pep")
        shutil.copyfile("orthoData/Skud_CDS.pep", "remove/len_30/Skud_CDS.pep")
        shutil.copyfile("orthoData/Smik_CDS.pep", "remove/len_30/Smik_CDS.pep")
        shutil.copyfile("orthoData/Spar_NCBI_CDS.pep", "remove/len_30/Spar_NCBI_CDS.pep")

        ######################

        # 40

        ######################

        os.mkdir('remove/len_40')

        temp_scer_records = Scer_records

        temp_scer_records = [rec for rec in temp_scer_records if rec.id != id]

        changed_CDS = change_seq(temp_CDS, percent = 40, lengthen = False)  

        temp_scer_records = list(temp_scer_records)
        temp_scer_records.append(changed_CDS)

        with open("remove/len_40/Scer_NCBI_CDS.pep","w") as f:
            SeqIO.write(temp_scer_records, f, "fasta")
        f.close()

        shutil.copyfile("orthoData/Sarb_CDS.pep", "remove/len_40/Sarb_CDS.pep")
        shutil.copyfile("orthoData/Sbay_CDS.pep", "remove/len_40/Sbay_CDS.pep")
        shutil.copyfile("orthoData/Skud_CDS.pep", "remove/len_40/Skud_CDS.pep")
        shutil.copyfile("orthoData/Smik_CDS.pep", "remove/len_40/Smik_CDS.pep")
        shutil.copyfile("orthoData/Spar_NCBI_CDS.pep", "remove/len_40/Spar_NCBI_CDS.pep")


        ######################

        # 50

        ######################


        os.mkdir('remove/len_50')

        temp_scer_records = Scer_records

        temp_scer_records = [rec for rec in temp_scer_records if rec.id != id]

        changed_CDS = change_seq(temp_CDS, percent = 50, lengthen = False)  

        temp_scer_records = list(temp_scer_records)
        temp_scer_records.append(changed_CDS)

        with open("remove/len_50/Scer_NCBI_CDS.pep","w") as f:
            SeqIO.write(temp_scer_records, f, "fasta")
        f.close()

        shutil.copyfile("orthoData/Sarb_CDS.pep", "remove/len_50/Sarb_CDS.pep")
        shutil.copyfile("orthoData/Sbay_CDS.pep", "remove/len_50/Sbay_CDS.pep")
        shutil.copyfile("orthoData/Skud_CDS.pep", "remove/len_50/Skud_CDS.pep")
        shutil.copyfile("orthoData/Smik_CDS.pep", "remove/len_50/Smik_CDS.pep")
        shutil.copyfile("orthoData/Spar_NCBI_CDS.pep", "remove/len_50/Spar_NCBI_CDS.pep")

        os.mkdir('append')
        os.mkdir('append/len_10')
        os.mkdir('append/len_20')
        os.mkdir('append/len_30')
        os.mkdir('append/len_40')
        os.mkdir('append/len_50')
