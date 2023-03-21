#!/bin/sh 


if [ ! -d "res_bench" ]
then 
    mkdir res_bench 
fi

if [ -d "temp" ]
then
    rm -rf temp
fi

abs_path=$(pwd)/res_bench

for dir in $(find Modified_data/ -maxdepth 4 -mindepth 4); do

    
    cond=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $2 }') 
    seq=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $3 }')
    action=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $4 }')
    length=$(echo  ${dir} | awk 'BEGIN { FS = "/" }; { print $5 }')

    
    orthofinder -f Modified_data/$cond/$seq/$action/$length -op -o temp

    mv temp/Results_Mar21/Log.txt $abs_path/${cond}_${seq}_${action}_${length}

    rm -rf temp


    # grep -Eo "[0-9]+"

done