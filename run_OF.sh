#!/bin/sh 


if [ ! -d "res_bench" ]
then 
    mkdir res_bench 
fi

if [ ! -d "temp" ]
then
    mkdir temp
fi

res_bench

i=1


for dir in $(find Modified_data/ -maxdepth 4 -mindepth 4); do

    
    cond=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $2 }') 
    seq=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $3 }')
    action=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $4 }')
    length=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $5 }')

    # Compute orthologs with OrthoFinder
    orthofinder -t 16 -f Modified_data/$cond/$seq/$action/$length -o ${cond}_${seq}_${action}_${length}

    mkdir res_bench/${cond}_${seq}_${action}_${length}
    mkdir res_bench/${cond}_${seq}_${action}_${length}/GeneCount

    # Use OrthoFinder tools to count how many genes per species for each Phylogenetic Orthogroup of node N0 ( enclude all species )
    python3 ~/utils/OrthoFinder/tools/orthogroup_gene_count.py ${cond}_${seq}_${action}_${length}/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv 

    mv ${cond}_${seq}_${action}_${length}/*/Phylogenetic_Hierarchical_Orthogroups/N0.GeneCount.csv res_bench/${cond}_${seq}_${action}_${length}/GeneCount

    python3 ~/utils/OrthoFinder/tools/create_files_for_hogs.py ${cond}_${seq}_${action}_${length}/* res_bench/${cond}_${seq}_${action}_${length}/HOGS_fa N0

    rm -rf ${cond}_${seq}_${action}_${length}

    

done