#!/bin/bash

cd ..

sed -i "28c\#define STOPCOMMAND 1
" "weSmoldyn.h"

sed -i "1c\corralsBit 1
" "corralsParams.txt"

sed -i "8c\nPart	20
" "dynamicsParams.txt"

sed -i "9c\reactBit	0
" "dynamicsParams.txt"

sed -i "10c\entryBit 0
" "dynamicsParams.txt"

sed -i "12c\densityBit 0
" "dynamicsParams.txt"

sed -i "1c\tau		50
" "WEParams.txt"

sed -i "2c\repsPerBin	50
" "WEParams.txt"

sed -i "4c\tauMax		10
" "WEParams.txt"

sed -i "6c\fluxBin	0
" "WEParams.txt"

sed -i "32c\#define FIXEDTIME 1
" "weSmoldyn.h"

make

cd scripts

for runNo in $2
do

for corralWidth in 2.1
do
for corralRate in 1 
do

SAVEPATH=/pub/robertbt/WELibsmolData

./makePubs2.pl $1corralWidth${corralWidth}corralRate${corralRate} ${runNo}

cd ..

mkdir $SAVEPATH/$1corralWidth${corralWidth}corralRate${corralRate}

sed -i "2c\corralsWidth 	${corralWidth}
" "corralsParams.txt"

sed -i "3c\corralsRate 		${corralRate}
" "corralsParams.txt"

./seedchange ${runNo} 0

cp weSmoldyn seedchange WEParams.txt dynamicsParams.txt binDefinitions.txt binParams.txt corralsParams.txt ISEED \
scripts/pubs/$1corralWidth${corralWidth}corralRate${corralRate}.pub \
$SAVEPATH/$1corralWidth${corralWidth}corralRate${corralRate}

./seedchange ${runNo} 1

cd scripts
done
done
done

cd /dfs3/pub/robertbt/WELibsmolData

for runNo in $2;
do

for corralWidth in 2.1
do
for corralRate in  1
do
cd $1corralWidth${corralWidth}corralRate${corralRate}
qsub ./$1corralWidth${corralWidth}corralRate${corralRate}.pub
cd ..
done
done
done