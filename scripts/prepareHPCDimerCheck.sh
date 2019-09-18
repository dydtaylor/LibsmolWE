#!/bin/bash
for bindR in 0.1 #$(seq .1 .1 1.0);
do

for unbindK in .1 .2 .5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000
do


SAVEPATH=/pub/robertbt/WELibsmolData

./makePubs.pl $1bind${bindR}unbind${unbindK}
cd ..

#make

mkdir $SAVEPATH/$1bind${bindR}unbind${unbindK}

sed -i "6c\bindR	${bindR}
" "dynamicsParams.txt"

sed -i "7c\unbindK	${unbindK}
" "dynamicsParams.txt"

./seedchange 1 0

cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED scripts/pubs/$1bind${bindR}unbind${unbindK}.pub $SAVEPATH/$1bind${bindR}unbind${unbindK}
cd scripts
done
done