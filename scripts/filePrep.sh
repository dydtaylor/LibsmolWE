#!/bin/bash

cd ..

sed -i "28c\#define STOPCOMMAND 1
" "weSmoldyn.h"

sed -i "32c\#define FIXEDTIME 0
" "weSmoldyn.h"

sed -i "1c\dt	.000001
" "dynamicsParams.txt"

sed -i "1c\tau		50
" "WEParams.txt"

sed -i "2c\repsPerBin	200
" "WEParams.txt"

sed -i "3c\initialReps 	1000
" "WEParams.txt"

sed -i "4c\tauMax		12000
" "WEParams.txt"

sed -i "6c\fluxBin	0
" "WEParams.txt"

sed -i "9c\reactBit	0
" "dynamicsParams.txt"

sed -i "10c\entryBit	0
" "dynamicsParams.txt"

sed -i "12c\densityBit	0
" "dynamicsParams.txt"

sed -i "1c\corralsBit 0
" "corralsParams.txt"

sed -i "2c\worldL	3
" "dynamicsParams.txt"

cd scripts