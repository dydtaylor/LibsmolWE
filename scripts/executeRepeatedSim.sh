#!/bin/bash
cd /pub/robertbt/WELibsmolData
for n in {1..20};
do

cd $1runNo${n}
qsub ./$1runNo${n}.pub
cd ..
done