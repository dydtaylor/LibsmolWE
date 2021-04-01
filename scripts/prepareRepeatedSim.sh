#!/bin/bash
for runNo in {1..20};
do

./makePubs.pl $1runNo${runNo}

cd ..

SAVEPATH=/pub/robertbt/WELibsmolData

mkdir $SAVEPATH/$1runNo${runNo}
./seedchange ${runNo} 0

cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED scripts/pubs/$1runNo${runNo}.pub $SAVEPATH/$1runNo${runNo}
./seedchange ${runNo} 1
cd scripts
done