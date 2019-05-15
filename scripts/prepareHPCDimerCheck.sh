#!/bin/bash
for bindR in .0000001 .0000002 .0000005 .000001 .000002 .000005 .00001 .00002 .00005 .0001 .0002 .0005 .001 .002 .005 .01 .02
for unbindK in 1 2 5
do

./makePubs.pl $1bind${bindR}unbind${unbindK}
cd ..

RUNPATH=/data/users/robertbt/WELibsmolData

#make
mkdir $RUNPATH/$1bind${bindR}unbind${unbindK}

sed -i "6c\bindR	${bindR}
" "dynamicsParams.txt"

sed -i "7c\unbindK    ${unbindK}
" "dynamicsParams.txt"

cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED scripts/pubs/1bind${bindR}unbind${unbindK}.pub $RUNPATH/1bind${bindR}unbind${unbindK}
cd scripts

done