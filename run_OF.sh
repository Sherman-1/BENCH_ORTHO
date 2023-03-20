#!/bin/sh 




abs_path=$(pwd)/res_bench



for dir in $(find Modified_data/ -maxdepth 4 -mindepth 4); do

    
    cond=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $2 }') 
    seq=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $3 }')
    action=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $4 }')
    length=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $5 }')

    
    orthofinder -f Modified_data/$cond/$seq/$action/$length -og -o temp

    

    # grep -Eo "[0-9]+"

done