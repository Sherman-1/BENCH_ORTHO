#!/bin/python3

import math
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def change_seq(record, percent = 10, lengthen = True):
    
    # Scer >NC_001133.9_-_14193-14435_2_nc_intergenic HCA score 0.31 for complete iORF
    iORF = "IVSGKGKVIMQWNLRYSHVQISATKSVVGTFRYQNFLSTTREHHYAKQIVAHTSDRKKTVPKTMTYEETSIRIFIINSMAIVSGKGKVIMQWNLRYSHVQISATKSVVGTFRYQNFLSTTREHHYAKQIVAHTSDRKKTVPKTMTYEETSIRIFIINSMAIVSGKGKVIMQWNLRYSHVQISATKSVVGTFRYQNFLSTTREHHYAKQIVAHTSDRKKTVPKTMTYEETSIRIFIINSMAIVSGKGKVIMQWNLRYSHVQISATKSVVGTFRYQNFLSTTREHHYAKQIVAHTSDRKKTVPKTMTYEETSIRIFIINSMA"

    temp_seq = str(record.seq)
    
    temp_seq = temp_seq[:-1] if temp_seq.endswith('*') else temp_seq

    number_of_char = math.ceil(percent * len(temp_seq) / 100)
    
    
    if lengthen == True:
        
        new_seq = Seq(str(temp_seq + iORF[0:number_of_char])+ "*")
    
    else:
        
        new_seq = Seq(str(temp_seq[:-number_of_char]) + "*")
    
    
    record = SeqRecord(seq = new_seq,
                       id = record.id)



    return record