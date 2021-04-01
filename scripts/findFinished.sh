#!/bin/bash

cd /dfs3/pub/robertbt/WELibsmolData/
printf "\n\n\n$1 \n" >> finishedSims.txt

for tau in 1 2 5 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000
do
for mtarg in 10 45 80 115 150 185 220 255 290 325 360 395 430 465 500 
do
for runNo in {1..3};
do

cd /dfs3/pub/robertbt/WELibsmolData/
cd $1Tau${tau}mTarg${mtarg}runNo${runNo}

if [ -e $1Tau${tau}mTarg${mtarg}runNo${runNo}.Time ]
#if [ -e WEParams.txt ]
then
	cd /dfs3/pub/robertbt/WELibsmolData/
	printf "$1Tau${tau}mTarg${mtarg}runNo${runNo}\n" >> finishedSims.txt

fi

done
done
done