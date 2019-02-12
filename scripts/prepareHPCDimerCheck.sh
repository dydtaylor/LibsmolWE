#!/bin/bash
for bindR in .0000001 .0000002 .0000005 .000001 .000002 .000005 .00001 .00002 .00005 .0001 .0002 .0005 .001 .002 .005 .01 .02
do

./makePubs.pl $1${n}
cd ..

RUNPATH=/data/users/robertbt/WELibsmolData

#make
mkdir $RUNPATH/$1${n}

sed -i "6c\bindR	${bindR}
" "dynamicsParams.txt"

cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED scripts/pubs/$1${n}.pub $RUNPATH/$1${n}
cd scripts

done