#!/bin/bash

make

sed -i "" "1c\ 
.000001 3 1 1 15" "dynamicsParams.txt"

for tau in 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100
do
sed -i "" "1c\ 
${tau} 200 10000 25" "WEParams.txt"
caffeinate ./weSmoldyn simOutTau${tau}N15.txt fluxTau${tau}N15.txt errTau${tau}N15.txt 0 >consoleOutTau${tau}N15.txt
rm consoleOutTau${tau}N15.txt
done