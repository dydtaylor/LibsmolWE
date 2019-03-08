#!/bin/bash
for n in {2..20};
do

./makePubs.pl $1${n}
cd ..

RUNPATH=/data/users/robertbt/WELibsmolData

#make
mkdir $RUNPATH/$1${n}

sed -i "" "5c\
nPart	${n}
" "dynamicsParams.txt"

cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED scripts/pubs/$1${n}.pub $RUNPATH/$1${n}
done