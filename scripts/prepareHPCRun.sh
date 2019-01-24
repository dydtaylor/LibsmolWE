#!/bin/bash

cd ..

RUNPATH=/data/users/robertbt/WELibsmolData

make
mkdir $RUNPATH/$1
./makePubs.pl $1

cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED scripts/singleHPCRun.pub $RUNPATH/$1