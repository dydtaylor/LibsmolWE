#!/bin/bash
#To make many easily accessed folders, no hash/date included
#allows for easier creation of multiple copies of similar runs


mkdir -p ../wesmoldynresults
for numParticles in 15 25
do
sed -i "5c${numParticles}" "mfiles/input/parameters.txt"
for tauVals in .8 .5 .3 .1 .08 .05 .03 .01 .008 .005 .003 .001 .0008 .0005 .0003 .0001 .00008 .00005 .00003 .00001 .000008 .000005 .000003 .000001
do
sed -i "10ctau=${tauVals}" "mfiles/runTr.sh"
for j in 0 1 2 3 4 5 6 7 8 9
do 
PREFIX="Sweep${numParticles}Tau${tauVals}r${j}.e1fix"
cp -r mfiles "../wesmoldynresults/Sweep/${PREFIX}"
done
done
done 