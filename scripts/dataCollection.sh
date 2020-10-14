#!/bin/bash

#This script should not be run without changing it. The purpose of this script is to combine multiple data collection scripts into one to save time on waiting for the terminal to execute job submission. EVERYTHING IS SUBJECT TO CHANGE

# Argv: $1 Series name $2 Number of runs to perform per parameter set. For recollection (and inside the scripts) 3 and 4 are
# 3 = custom date (if 0 then don't use custom date, if > 0 then use custom date) 
#4 = custom date string
#5 (only simple diffusion) world L

cd ..

sed -i "28c\#define STOPCOMMAND 1
" "weSmoldyn.h"

sed -i "32c\#define FIXEDTIME 0
" "weSmoldyn.h"

sed -i "1c\dt	.000001
" "dynamicsParams.txt"

sed -i "1c\tau		50
" "WEParams.txt"

sed -i "2c\repsPerBin	50
" "WEParams.txt"

sed -i "3c\initialReps 	500
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

sed -i "2c\worldL	4
" "dynamicsParams.txt"
make

cd scripts
#Execute beginning scripts. constantNSweepL might need to be changed depending on box size.
./simpleDiffusionSweep.sh $1 $2 0 null 4
./simpleDiffusionConstantNSweepL.sh $1 $2 0 null
./constantDensityScan.sh $1 $2 0 null


##Change parameters + binning that is different for each sim. Use binDefinitions2.txt for N>100
cd ..
sed -i "2c\worldL	4
" "dynamicsParams.txt"

sed -i "8c\nPart	128
" "dynamicsParams.txt"

cp binDefinitions2.txt binDefinitions.txt

cd scripts

./reentryRateScan.sh $1 $2 0 null
./repeatedDimerCheck.sh $1 $2 0 null

##Return binning to low N binning.
cd ..
cp binDefinitions1.txt binDefinitions.txt
cd scripts

./corralsSweep.sh  $1 $2 0 null

#Simple diffusion scripts will go here and follow the same convention. They are separated right now because the scripts need to be changed.

#./simpleDiffusion.sh seriesname

#Other misc. things that need to be changed in individual scripts:

#Dimers Would like to replace the hard coding of n63 in the FULLNAME variable. (and also in sed commands)
#Constant Density: Would like to get rid of hard coding of rho7 in FULLNAME variable (and also in sed commands).

#All the bin definitions need to be more fluid. This can definitely be handled with the binparams.txt, but for runs that have various bindefinitions it would be nice to have a smoother way to modify binDefinitions.txt

#Desired features to include:
#Makepubs.pl includes an rsync command to push towards the imac / desired server. Could potentially take that server as an input argument.

#Move the edits of dynamics params / wesmoldyn / corrals params /etc. that are going to be universal for all the runs in the data collection example to this file and have only the parameters that are changed be edited in the sub scripts.