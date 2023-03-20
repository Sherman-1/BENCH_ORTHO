#!/bin/sh 


#Select sequences of length greater than given threshold in OrthoFINDER output files ( OG0000XXX.fa )



wd=$(pwd)


if [ -z "$1" ]
then

	SingleDir="Single_Copy_Orthologue_Sequences/"

else

	SingleDir=$1

fi

cd $SingleDir


# Keep touch of which OGroup is used for benchmarking
touch OGroup_100.txt
touch OGroup_200.txt
touch OGroup_400.txt


for f in *.fa; do

	temp=$(faSize $f)
	length=$(echo $temp | grep -o "mean \([0-9]*\)" | head -n 1 | grep -o "\([0-9]*\)")
	sd=$(echo $temp | grep -o "sd \([0-9]*\)" | head -n 1 | grep -o "\([0-9]*\)")
		
	if [ $length -ge 190 ] && [ $length -le 200 ] && [ $sd -le 10 ]; then
		
		echo $f >> OGroup_200.txt


	elif [ $length -ge 390 ] && [ $length -le 400 ] && [ $sd -le 20 ]; then

		echo $f >> OGroup_400.txt 

	elif [ $length -ge 90 ] && [ $length -le 100 ] && [ $sd -le 5 ]; then

		echo $f >> OGroup_100.txt

	fi

done 

mv OGroup*.txt $wd
cd $wd









