#!/bin/bash

#This script should not be run without changing it. The purpose of this script is to combine multiple data collection scripts into one to save time on waiting for the terminal to execute job submission. EVERYTHING IS SUBJECT TO CHANGE. This script runs assuming you are attempting to rerun simulations that did not finish / execute.

# Argv: $1 Series name $2 Number of runs to perform per parameter set.
# 3 = custom date (if 0 then don't use custom date, if > 0 then use custom date) 
#4 = custom date string

./simpleDiffusionSweep.sh $1 $2 $3 $4
./simpleDiffusionConstantNSweepL.sh $1 $2 $3 $4
./constantDensityScan.sh $1 $2 $3 $4
./reentryRateScan.sh $1 $2 $3 $4
./repeatedDimerCheck.sh $1 $2 $3 $4
./corralsSweep.sh  $1 $2 $3 $4

#Simple diffusion scripts will go here and follow the same convention. They are separated right now because the scripts need to be changed.

#./simpleDiffusion.sh seriesname

#Other misc. things that need to be changed in individual scripts:

#Dimers Would like to replace the hard coding of n63 in the FULLNAME variable. (and also in sed commands)
#Constant Density: Would like to get rid of hard coding of rho7 in FULLNAME variable (and also in sed commands).

#All the bin definitions need to be more fluid. This can definitely be handled with the binparams.txt, but for runs that have various bindefinitions it would be nice to have a smoother way to modify binDefinitions.txt

#Desired features to include:
#Makepubs.pl includes an rsync command to push towards the imac / desired server. Could potentially take that server as an input argument.