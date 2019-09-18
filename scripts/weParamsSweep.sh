#!/bin/bash

##Adjust dynamics parameters + wesmoldyn + WE parameters
## Sets number of particles to be 20, sets "if all particles clear, stop sim" command in libsmol, enables reactions, sets flux bin
cd ..
sed -i "28c\#define STOPCOMMAND 1
" "weSmoldyn.h"

sed -i "25c\#define KSCRITICAL .02
" "weSmoldyn.h"

sed -i "26c\#define DKSCRITICAL 0.3
" "weSmoldyn.h"

sed -i "32c\#define FIXEDTIME 0
" "weSmoldyn.h"

sed -i "8c\nPart	20
" "dynamicsParams.txt"

sed -i "9c\reactBit	0
" "dynamicsParams.txt"

sed -i "6c\fluxBin	0
" "WEParams.txt"

sed -i "7c\ksNT		300
" "WEParams.txt"

#dtMax=100000000

##Make the executable
make
cd scripts
SAVEPATH=/pub/robertbt/WELibsmolData

##Copy executable for each run, set save paths
for tau in 1 2 5 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000
do
for mtarg in 10 45 80 115 150 185 220 255 290 325 360 395 430 465 500 
do
for runNo in {1..6};
do
cd ..

#tauMax=$(expr $dtMax / $tau / $mtarg)

#sed -i "4c\tauMax		${tauMax}
#" "WEParams.txt"

sed -i "1c\tau			${tau}
" "WEParams.txt"

sed -i "2c\repsPerBin	${mtarg}
" "WEParams.txt"
cd scripts



./makePubs.pl $1Tau${tau}mTarg${mtarg}runNo${runNo}
cd ..

#make

mkdir $SAVEPATH/$1Tau${tau}mTarg${mtarg}runNo${runNo}

./seedchange ${runNo} 0

cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED scripts/pubs/$1Tau${tau}mTarg${mtarg}runNo${runNo}.pub $SAVEPATH/$1Tau${tau}mTarg${mtarg}runNo${runNo}
./seedchange ${runNo} 0
cd scripts
done
done
done

#Submit jobs to hpc queue


for tau in 1 2 5 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000
do
for mtarg in 10 45 80 115 150 185 220 255 290 325 360 395 430 465 500 
do
for runNo in {1..6};
do
cd ..

cd $SAVEPATH/$1Tau${tau}mTarg${mtarg}runNo${runNo}
qsub ./$1Tau${tau}mTarg${mtarg}runNo${runNo}.pub
done
done
done
