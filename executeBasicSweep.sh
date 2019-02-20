#!/bin/bash
cd /data/users/robertbt/WELibsmolData
for n in {2..20};
do

cd $1${n}
qsub ./$1${n}.pub
cd ..
done