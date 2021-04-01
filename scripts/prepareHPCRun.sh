#!/bin/bash

./makePubs.pl $1
cd ..

RUNPATH=/pub/robertbt/WELibsmolData

#make
mkdir $RUNPATH/$1

./seedchange
cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED scripts/pubs/$1.pub $RUNPATH/$1