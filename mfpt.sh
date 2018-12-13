#!/bin/bash

make
sed -i "" "1c\ 
10 200 10000 25" "WEParams.txt"

for n in 13 19
do
sed -i "" "1c\ 
.000001 3 1 1 ${n}" "dynamicsParams.txt"
caffeinate ./weSmoldyn simOut${n}.txt flux${n}.txt err${n}.txt 0 >consoleOut${n}.txt
done