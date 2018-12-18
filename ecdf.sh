#!/bin/bash

make

sed -i "" "1c\ 
10 500 200 10000 25" "WEParams.txt"
sed -i "" "1c\ 
.000001 3 1 1 15" "dynamicsParams.txt"

for run in {1..10}
do
caffeinate ./weSmoldyn ecdfSim${run}.txt ecdfFlux${run}.txt ecdfErr${run}.txt 0 >ecdfOut${run}.txt
mv ecdfSim${run}.txt ECDFs/NR500/ecdfSim${run}.txt
mv ecdfFlux${run}.txt ECDFs/NR500/ecdfFlux${run}.txt
mv ecdfErr${run}.txt ECDFs/NR500/ecdfErr${run}.txt
rm ecdfOut${run}.txt
done

for tau in 1 2 5 10 20 50 100
do
sed -i "" "1c\ 
${tau} 500 200 10000 25" "WEParams.txt"
caffeinate ./weSmoldyn simOutTau${tau}N15.txt fluxTau${tau}N15.txt errTau${tau}N15.txt 0 >consoleOutTau${tau}N15.txt
rm consoleOutTau${tau}N15.txt
mv simOutTau${tau}N15.txt tauSweep/simOutTau${tau}N15.txt
mv fluxTau${tau}N15.txt tauSweep/fluxTau${tau}N15.txt
mv errTau${tau}N15.txt tauSweep/errTau${tau}N15.txt
done