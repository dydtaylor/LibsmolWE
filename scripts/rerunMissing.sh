#!/bin/bash
cd /data/users/robertbt/WELibsmolData

for bindR in .001 .002 .003 .01 .02 .05 .1 .2 .5 1
for unbindk in .1 .2 .5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000

do
cd $1bind${bindR}unbind${unbindk}
if [! -e mCounts.txt]
qsub ./$1bind${bindR}unbind{unbindk}
cd ..
done