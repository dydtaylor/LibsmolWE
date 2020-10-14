#!/bin/bash

cd ..


sed -i "10c\entryBit 1
" "dynamicsParams.txt"

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

for runNo in $2;
do

for entryRate in 0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.0

#0 0.025 0.05 0.075 0.10 0.125 0.15 0.175 0.20 0.225 0.25 0.275 0.30 0.325 0.35 0.375 0.40 0.425 0.45 0.475 0.50 0.525 0.55 0.575 0.60 0.625 0.65 0.675 0.70 0.725 0.75 0.775 0.80 0.825 0.85 0.875 0.90 0.925 0.95 0.975 1.0
do

SAVEPATH=/pub/robertbt/WELibsmolData
FULLNAME=$1.$d.Reentry.EntryRate${entryRate}

./makePubs.pl $FULLNAME ${runNo}
cd ..

#make

mkdir $SAVEPATH/$FULLNAME

sed -i "11c\entryRate	${entryRate}
" "dynamicsParams.txt"


./seedchange ${runNo} 0

cp weSmoldyn seedchange WEParams.txt dynamicsParams.txt binDefinitions.txt binParams.txt corralsParams.txt ISEED scripts/pubs/$FULLNAME.pub $SAVEPATH/$FULLNAME


cd scripts
done
done


cd ..
sed -i "10c\entryBit 0
" "dynamicsParams.txt"

cd /dfs3/pub/robertbt/WELibsmolData

for runNo in $2;
do

for entryRate in 0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.0
do
FULLNAME=$1.$d.Reentry.EntryRate${entryRate}

cd $FULLNAME
qsub ./$FULLNAME.pub
cd ..
done
done