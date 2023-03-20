#!/bin/python3

import math
from Bio.Seq import Seq

def change_seq(record, percent = 10, lengthen = True):
    
    # Scer >NC_001133.9_-_14193-14435_2_nc_intergenic HCA score 0.31 for complete iORF
    iORF = "IVSGKGKVIMQWNLRYSHVQISATKSVVGTFRYQNFLSTTREHHYAKQIVAHTSDRKKTVPKTMTYEETSIRIFIINSMAIVSGKGKVIMQWNLRYSHVQISATKSVVGTFRYQNFLSTTREHHYAKQIVAHTSDRKKTVPKTMTYEETSIRIFIINSMA"

    temp_seq = record.seq
    
    
    number_of_char = math.ceil(percent * len(temp_seq) / 100)
    
    
    if lengthen == True:
        
        new_seq = Seq(temp_seq + iORF[0:number_of_char])
    
    else:
        
        new_seq = Seq(temp_seq[:-number_of_char])
    
    
    record.seq = new_seq


    return record