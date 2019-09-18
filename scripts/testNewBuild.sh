#!/bin/bash

cd ..
sed -i "28c\#define STOPCOMMAND 1
" "weSmoldyn.h"

sed -i "25c\#define KSCRITICAL .02
" "weSmoldyn.h"

sed -i "26c\#define DKSCRITICAL 0.3
" "weSmoldyn.h"

sed -i "32c\#define FIXEDTIME 1
" "weSmoldyn.h"

sed -i "8c\nPart	20
" "dynamicsParams.txt"

sed -i "9c\reactBit	0
" "dynamicsParams.txt"

sed -i "6c\fluxBin	0
" "WEParams.txt"

sed -i "7c\ksNT		300
" "WEParams.txt"

sed -i "4c\tauMax 2000
" "WEParams.txt"

make
cd scripts
SAVEPATH=/pub/robertbt/WELibsmolData/TestDirectory

for tau in 30
do
for mtarg in 100
do
for runNo in {1..3};
do
cd ..
sed -i "1c\tau			${tau}
" "WEParams.txt"

sed -i "2c\repsPerBin	${mtarg}
" "WEParams.txt"
cd scripts
./makePubs.pl $1runNo${runNo}
cd ..

#make

mkdir $SAVEPATH/$1runNo${runNo}

./seedchange ${runNo} 0

cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED scripts/pubs/$1runNo${runNo}.pub $SAVEPATH/$1runNo${runNo}
./seedchange ${runNo} 0
cd scripts
done
done
done

for tau in 30
do
for mtarg in 100
do
for runNo in {1..3};
do
cd $SAVEPATH/$1runNo${runNo}
qsub ./$1runNo${runNo}.pub
done
done
done
