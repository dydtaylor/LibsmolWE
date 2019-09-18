#!/bin/bash
cd /dfs3/pub/robertbt/WELibsmolData/
for tau in 1 2 5 10 20 50 100 200 500
do
for mtarg in 10 80 150 220 290 360 430 500
do
for run in {1..3}
do

cd /dfs3/pub/robertbt/WELibsmolData/
cd ${1}Tau${tau}mTarg${mtarg}runNo${run}
rm mCountsWeighted.txt
cd ..

done
done
done