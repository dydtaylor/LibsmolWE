#!/bin/bash

cd ..

RUNPATH=/Users/allardlab/LibsmolWE/Runs/

make
mkdir $RUNPATH/$1
cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED $RUNPATH/$1