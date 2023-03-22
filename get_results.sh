#!/bin/sh


# Store which are the Orthologs of the 9 sequences used for benchmark

for file in rna*.txt; do

    cond=$(grep -Po "\\d+" $file)

    for id in $(cat $file); do

        grep -r $id Single_Hog

    done

done

