#!/bin/bash

make

sed -i "" "1c\ 
10 200 10000 25" "WEParams.txt"
sed -i "" "1c\ 
.000001 3 1 1 15" "dynamicsParams.txt"

for run in {1..10}
do
caffeinate ./weSmoldyn ecdfSim${run}.txt ecdfFlux${run}.txt ecdfErr${run}.txt 0 >ecdfOut${run}.txt
rm ecdfOut${run}.txt
done
