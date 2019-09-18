#!/bin/bash

cd /dfs3/pub/robertbt/WELibsmolData/
printf "\n\n\n$1 \n" >> missingData.txt

for tau in 1 2 5 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000
do
for mtarg in 10 45 80 115 150 185 220 255 290 325 360 395 430 465 500 
do
for runNo in {1..6};
do

cd /dfs3/pub/robertbt/WELibsmolData/
cd $1Tau${tau}mTarg${mtarg}runNo${runNo}

if [ ! -e $1Tau${tau}mTarg${mtarg}runNo${runNo}.Time ]
#if [ -e WEParams.txt ]
then
if [ -e $1Tau${tau}mTarg${mtarg}runNo${runNo}.Out ] #
then

rm $1Tau${tau}mTarg${mtarg}runNo${runNo}.Flux
rm $1Tau${tau}mTarg${mtarg}runNo${runNo}.Out
rm $1Tau${tau}mTarg${mtarg}runNo${runNo}.pub
rm ksOut.txt
rm dualKS.txt
cd /data/users/robertbt/libsmolWEDimerCode/LibsmolWE/scripts
./makePubs.pl $1Tau${tau}mTarg${mtarg}runNo${runNo}
cp pubs/$1Tau${tau}mTarg${mtarg}runNo${runNo}.pub /dfs3/pub/robertbt/WELibsmolData/$1Tau${tau}mTarg${mtarg}runNo${runNo}
cd /dfs3/pub/robertbt/WELibsmolData/$1Tau${tau}mTarg${mtarg}runNo${runNo}
qsub ./$1Tau${tau}mTarg${mtarg}runNo${runNo}.pub

	cd /dfs3/pub/robertbt/WELibsmolData/
	printf "$1Tau${tau}mTarg${mtarg}runNo${runNo}\n" >> missingData.txt
fi
fi

done
done
done