#!/bin/bash

cd ..

sed -i "28c\#define STOPCOMMAND 1
" "weSmoldyn.h"

sed -i "8c\nPart	20
" "dynamicsParams.txt"

sed -i "9c\reactBit	1
" "dynamicsParams.txt"

sed -i "1c\tau		100
" "WEParams.txt"

sed -i "2c\repsPerBin	300
" "WEParams.txt"

sed -i "4c\tauMax		12000
" "WEParams.txt"

sed -i "6c\fluxBin	0
" "WEParams.txt"

sed -i "32c\#define FIXEDTIME 1
" "weSmoldyn.h"
make

cd scripts


for runNo in {1..3};
do

for bindR in 0.0001 0.001 0.01 #$(seq .1 .1 1.0);
do

for unbindK in .001 .002 .005 .01 .02 .05 .1 .2 .5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000
do


SAVEPATH=/pub/robertbt/WELibsmolData

./makePubs.pl $1bind${bindR}unbind${unbindK}runNo${runNo}
cd ..

#make

mkdir $SAVEPATH/$1bind${bindR}unbind${unbindK}runNo${runNo}

sed -i "6c\bindR	${bindR}
" "dynamicsParams.txt"

sed -i "7c\unbindK	${unbindK}
" "dynamicsParams.txt"

./seedchange ${runNo} 0

cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED scripts/pubs/$1bind${bindR}unbind${unbindK}runNo${runNo}.pub $SAVEPATH/$1bind${bindR}unbind${unbindK}runNo${runNo}
./seedchange ${runNo} 1
cd scripts
done
done
done

cd /dfs3/pub/robertbt/WELibsmolData

for runNo in {1..3};
do

for bindR in 0.0001 0.001 0.01 #$(seq .1 .1 1.0);
do

for unbindK in .001 .002 .005 .01 .02 .05 .1 .2 .5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000
do
cd $1bind${bindR}unbind${unbindK}runNo${runNo}
qsub ./$1bind${bindR}unbind${unbindK}runNo${runNo}.pub
cd ..
done
done
done