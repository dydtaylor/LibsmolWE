#!/bin/bash
for n in {2..25};
do

./makePubs.pl $1${n}
cd ..

SAVEPATH=/pub/robertbt/WELibsmolData

#make
mkdir $SAVEPATH/$1${n}

sed -i "8c\nPart	${n}
" "dynamicsParams.txt"

for i in {1..10}
do
./seedchange $2
done

cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED scripts/pubs/$1${n}.pub $SAVEPATH/$1${n}
cd scripts
done