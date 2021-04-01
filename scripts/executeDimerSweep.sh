#!/bin/bash
cd /pub/robertbt/WELibsmolData
for bindR in 0.1 #$(seq .1 .1 1.0);
do

for unbindK in .1 .2 .5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 #20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000 100000000 200000000 500000000
do 

cd $1bind${bindR}unbind${unbindK}
qsub ./$1bind${bindR}unbind${unbindK}.pub
cd ..
done
done