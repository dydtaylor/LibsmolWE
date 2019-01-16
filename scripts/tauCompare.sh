#!/bin/bash

make

sed -i "" "1c\ 
.000001 3 1 1 15" "dynamicsParams.txt"

for tau in 50 100
do
sed -i "" "1c\ 
${tau} 500 200 10000 25" "WEParams.txt"
caffeinate ./weSmoldyn simOutTau${tau}N15.txt fluxTau${tau}N15.txt errTau${tau}N15.txt 0 &>/dev/null
rm consoleOutTau${tau}N15.txt
done