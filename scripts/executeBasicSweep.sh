#!/bin/bash
cd /pub/robertbt/WELibsmolData
for n in {105..280..5}
do
for runNo in {1..10}
do

cd $1${n}runNo${runNo}
qsub ./$1${n}runNo${runNo}.pub
cd ..
done
done