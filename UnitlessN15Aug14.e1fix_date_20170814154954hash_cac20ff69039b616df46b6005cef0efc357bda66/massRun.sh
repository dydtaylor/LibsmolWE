#!/bin/bash
#Run a massive amount of jobs simultaneously. Jobs SHOULD NOT MAKE EXTENSIVE OUTPUTS
for numParticles in 15 25
do
for tauVals in .003 .5 .3 .1 .08 .05 .03 .01 .008 .005 .003 .001 .0008 .0005 .0003 .0001 .00008 .00005 .00003 .00001 .000008 .000005 .000003 .000001
do
for j in 0 1 2 3 4 5 6 7 8 9
do
PREFIX="Sweep${numParticles}Tau${tauVals}r${j}.e1fix"
cd "${PREFIX}"
qsub runTr.sh
cd ..
done
done
done
