#!/bin/bash

make
sed -i "" "1c\ 
10 500 200 100000 25" "WEParams.txt"

for n in {1..4}
do
caffeinate ./weSmoldyn 100ksimOut${n}.txt 100kflux${n}.txt 100kerr${n}.txt 0 100kExcTime${n}.txt >&/dev/null
done