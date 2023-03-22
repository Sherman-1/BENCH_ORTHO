#!/bin/sh 

PATH=$PBS_O_PATH
cd $PBS_O_WORKDIR
 

if [ ! -d "res_bench" ]
then 
    mkdir res_bench 
fi

time=$(date +"%T"| sed s/://g)
date=$(date | awk '{print $2 $3}')

touch run_${time}_${date}.log


echo "-------------------------------------"
echo "Current working directory : $(pwd)"
echo "-------------------------------------"


for dir in $(find Modified_data/ -maxdepth 4 -mindepth 4); do

   
    
   
    cond=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $2 }') 
    seq=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $3 }')
    action=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $4 }')
    length=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $5 }')



    # Compute orthologs with OrthoFinder

    echo "---------------------------------" >> run_${time}_${date}.log
    echo "Running OrthoFinder for ${seq} ${action} ${length}" >> run_${time}_${date}.log
    echo "---------------------------------" >> run_${time}_${date}.log

	
    orthofinder -t 20 -f Modified_data/$cond/$seq/$action/$length -o ${cond}_${seq}_${action}_${length}

    echo "\nOrthofinder done ..." >> run_${time}_${date}.log

    echo "\n Creating resbench directories" >> run_${time}_${date}.log

    mkdir res_bench/${cond}_${seq}_${action}_${length}
    mkdir res_bench/${cond}_${seq}_${action}_${length}/GeneCount 



    echo "\n Computing GeneCounts from Hierarchical orthogroups" >> run_${time}_${date}.log

    # Use OrthoFinder tools to count how many genes per species for each Phylogenetic Orthogroup of node N0 ( enclude all species )
    python3 /home/simon.herman/OrthoFinder/tools/orthogroup_gene_count.py ${cond}_${seq}_${action}_${length}/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv 

    echo "\n Moving GeneCount file to res_bench" >> run_${time}_${date}.log

    mv ${cond}_${seq}_${action}_${length}/*/Phylogenetic_Hierarchical_Orthogroups/N0.GeneCount.csv res_bench/${cond}_${seq}_${action}_${length}/GeneCount

    echo "\nWriting mFASTA files from Hierarchical Orthogroups ... " >> run_${time}_${date}.log

    python3 /home/simon.herman/OrthoFinder/tools/create_files_for_hogs.py ${cond}_${seq}_${action}_${length}/* res_bench/${cond}_${seq}_${action}_${length}/HOGS_fa N0

    rm -rf ${cond}_${seq}_${action}_${length}
	
    echo " ${cond}_${seq}_${action}_${length} done " >> run_${time}_${date}.log
    echo " -------------------------------------------" >> run_${time}_${date}.log
    echo " -------------------------------------------" >> run_${time}_${date}.log

done