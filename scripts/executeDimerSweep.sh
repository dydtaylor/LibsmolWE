#!/bin/bash
cd /data/users/robertbt/WELibsmolData
for bindR in .0000001 .0000002 .0000005 .000001 .000002 .000005 .00001 .00002 .00005 .0001 .0002 .0005 .001 .002 .005 .01 .02;
for unbindK in 1 2 5
do

cd $1bind${bindR}unbind
qsub ./$1bind${bindR}unbind${unbindK}.pub
cd ..
done