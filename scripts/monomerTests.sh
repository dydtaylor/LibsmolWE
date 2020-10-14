#!/bin/bash

cd ..

sed -i "28c\#define STOPCOMMAND 1
" "weSmoldyn.h"

sed -i "9c\reactBit	0
" "dynamicsParams.txt"

sed -i "11c\entryRate	1.0
" "dynamicsParams.txt"

sed -i "12c\densityBit 0
" "dynamicsParams.txt"

sed -i "1c\tau		50
" "WEParams.txt"

sed -i "2c\repsPerBin	300
" "WEParams.txt"

sed -i "4c\tauMax		12000
" "WEParams.txt"

sed -i "5c\nBins 20
" "WEParams.txt"

sed -i "6c\fluxBin	0
" "WEParams.txt"

sed -i "7c\ksNT		300
" "WEParams.txt"

sed -i "32c\#define FIXEDTIME 0
" "weSmoldyn.h"

sed -i "14c\monomerStart 1
" "dynamicsParams.txt"

make



cd scripts

for runNo in {1..20};
do

for nPart in {2..19}
do


SAVEPATH=/pub/robertbt/WELibsmolData

./makePubs.pl $1N${nPart}runNo${runNo}
cd ..

#make

mkdir $SAVEPATH/$1N${nPart}runNo${runNo}

sed -i "8c\nPart	${nPart}
" "dynamicsParams.txt"

./seedchange ${runNo} 0

cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED scripts/pubs/$1N${nPart}runNo${runNo}.pub $SAVEPATH/$1N${nPart}runNo${runNo}
./seedchange ${runNo} 1
cd scripts
done
done


cd ..

cd /dfs3/pub/robertbt/WELibsmolData

for runNo in {1..20};
do
for nPart in {2..19}
do
cd $1N${nPart}runNo${runNo}
qsub ./$1N${nPart}runNo${runNo}.pub
cd ..
done
done