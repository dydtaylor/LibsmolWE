#!/bin/bash

cd ..

sed -i "1c\corralsBit 1
" "corralsParams.txt"

sed -i "8c\nPart	20
" "dynamicsParams.txt"

sed -i "2c\repsPerBin	100
" "WEParams.txt"

sed -i "3c\initialReps 	5000
" "WEParams.txt"


make

#Date variable. 
if [ 1 -gt $3 ];
then
d=`date +%m.%d.%Y`
fi

if [ $3 -gt 0 ];
then
d=$4
fi

cd scripts

for runNo in $2
do

for corralWidth in 2 2.25 2.5 2.75 3
do
for corralRate in 0 0.00001 0.00002 0.00005 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1
do

SAVEPATH=/pub/robertbt/WELibsmolData

FULLNAME=$1.$d.Corrals.width${corralWidth}.corralRate${corralRate}

./makePubs.pl $FULLNAME ${runNo}

cd ..

mkdir $SAVEPATH/$FULLNAME

sed -i "2c\corralsWidth 	${corralWidth}
" "corralsParams.txt"

sed -i "3c\corralsRate 		${corralRate}
" "corralsParams.txt"

./seedchange ${runNo} 0

cp weSmoldyn seedchange WEParams.txt dynamicsParams.txt binDefinitions.txt binParams.txt corralsParams.txt ISEED scripts/pubs/$FULLNAME.pub $SAVEPATH/$FULLNAME



cd scripts
done
done
done

cd ..
sed -i "1c\corralsBit 0
" "corralsParams.txt"

cd /dfs3/pub/robertbt/WELibsmolData

for runNo in $2
do

for corralWidth in 2 2.25 2.5 2.75 3
do

for corralRate in 0 0.00001 0.00002 0.00005 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1
do

FULLNAME=$1.$d.Corrals.width${corralWidth}.corralRate${corralRate}

cd $FULLNAME
qsub ./$FULLNAME.pub
cd ..

done

done

done

cd /data/users/robertbt/libsmolWEDimerCode/LibsmolWE

sed -i "2c\repsPerBin	50
" "WEParams.txt"

sed -i "3c\initialReps 	500
" "WEParams.txt"