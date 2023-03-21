#!/bin/sh 


python3 ~/utils/OrthoFinder/tools/orthogroup_gene_count.py ${cond}_${seq}_${action}_${length}/* 

python3 OrthoFinder/tools/create_files_for_hogs.py Run_2/Results_Mar10/ N0_fasta N0

awk 'NR==FNR{a[$1]++; next} $1 in a' test.txt Run_2/Results_Mar10/Phylogenetic_Hierarchical_Orthogroups/N0.ts



    # Use OrthoFinder tools to count how many genes per species for each Phylogenetic Orthogroup of node N0 ( enclude all species )
    python3 ~/utils/OrthoFinder/tools/orthogroup_gene_count.py ${cond}_${seq}_${action}_${length}/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv 

    # Keep name of HOG that have only 1 gene per specie ( = Single Orthologues groups )
    SOGs=$(awk -F'[\t]+' '{if (($8 == 1) && ($4 == 1) && ($5 == 1) && ($6 == 1) && ($7 == 1) && ($9 == 1)) print $1}' ${cond}_${seq}_${action}_${length}/*/Phylogenetic_Hierarchical_Orthogroups/N0.GeneCount.csv)

    python3 OrthoFinder/tools/create_files_for_hogs.py Run_2/Results_Mar10/ . N0




    for SOG in $SOGs;
    do

        mv N0_fasta
