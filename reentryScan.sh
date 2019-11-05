#!/bin/bash

sed -i "" "28c\ 
#define STOPCOMMAND 1
" "weSmoldyn.h"

sed -i "" "8c\ 
nPart	20
" "dynamicsParams.txt"

sed -i "" "9c\ 
reactBit	0
" "dynamicsParams.txt"

sed -i "" "10c\ 
entryBit 1
" "dynamicsParams.txt"

sed -i "" "1c\ 
tau		50
" "WEParams.txt"

sed -i "" "2c\ 
repsPerBin	100
" "WEParams.txt"

sed -i "" "4c\ 
tauMax		12000
" "WEParams.txt"

sed -i "" "6c\ 
fluxBin	0
" "WEParams.txt"

sed -i "" "32c\ 
#define FIXEDTIME 0
" "weSmoldyn.h"
make


for runNo in 1;
do

for entryRate in 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.0
do

sed -i "" "11c\ 
entryRate ${entryRate}
" "dynamicsParams.txt"

./seedchange ${runNo} 0

./weSmoldyn reentryRate${entryRate}run${runNo}.out reentryRate${entryRate}run${runNo}.Flux reentryRate${entryRate}run${runNo}.seed 0 reentryRate${entryRate}run${runNo}.Time &>/dev/null

mv bin1.txt ReentryData/bin1EntryRate${entryRate}runNo${runNo}.txt
mv monomerLocs.txt ReentryData/monomerLocs${entryRate}runNo${runNo}.txt
mv dimerLocs.txt ReentryData/dimerLocs${entryRate}runNo${runNo}.txt
mv ksOut.txt ReentryData/dkOut${entryRate}runNo${runNo}.txt
mv dualKS.txt ReentryData/dualKS.txt${entryRate}runNo${runNo}.txt
mv reentryRate${entryRate}run${runNo}.out ReentryData/reentryRate${entryRate}run${runNo}.out
mv reentryRate${entryRate}run${runNo}.Flux ReentryData/reentryRate${entryRate}run${runNo}.Flux
mv reentryRate${entryRate}run${runNo}.seed ReentryData/reentryRate${entryRate}run${runNo}.seed
mv reentryRate${entryRate}run${runNo}.Time ReentryData/reentryRate${entryRate}run${runNo}.Time

done
done