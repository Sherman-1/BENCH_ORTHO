#!/bin/sh


wd="Single_Copy_Orthologue_Sequences/"

# Change to folder with selected sequences to change
current_dir=$(pwd)


# Randomly select 3 Orthogroups for each of the different lengths
s100=$(shuf -n 3 OGroup_100.txt)
s200=$(shuf -n 3 OGroup_200.txt)
s400=$(shuf -n 3 OGroup_400.txt)


# Keep touch of selected RNAs to extend 
touch rna100.txt
touch rna200.txt
touch rna400.txt


# Search for S.cere rna ID from the randomly selected orthoFINDER output of given length
cd $wd
for s in $s100; do

	echo $(grep "rna-NM*" $s) >> rna100.txt

done


for s in $s200; do

        echo $(grep "rna-NM*" $s) >> rna200.txt

done


for s in $s400; do

        echo $(grep "rna-NM*" $s) >> rna400.txt

done






mv rna*.txt $current_dir
cd $current_dir

sed -i -e 's/\[translate([0-9]*)\]//g' orthoData/Scer_NCBI_CDS.pep







