#!/bin/bash

cd /dfs3/pub/robertbt/WELibsmolData/
printf "\n\n\n$1 \n" >> missingData.txt

#for tau in 1 2 5 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000
#do
#for mtarg in 10 45 80 115 150 185 220 255 290 325 360 395 430 465 500 
#do

for bindR in 0.01 #$(seq .1 .1 1.0);
do

for unbindK in .001 .002 .005 .01 .02 .05 .1 .2 .5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000
do
for runNo in {1..20};
do

cd /dfs3/pub/robertbt/WELibsmolData/
#cd $1Tau${tau}mTarg${mtarg}runNo${runNo}
cd $1bind${bindR}unbind${unbindK}runNo${runNo}
if [ ! -e $1bind${bindR}unbind${unbindK}runNo${runNo}.Time ]
#if [ -e WEParams.txt ]
then
if [ -e $1bind${bindR}unbind${unbindK}runNo${runNo}.Out ] #
then

rm $1bind${bindR}unbind${unbindK}runNo${runNo}.Flux
rm $1bind${bindR}unbind${unbindK}runNo${runNo}.Out
#rm $1bind${bindR}unbind${unbindK}runNo${runNo}.pub
rm ksOut.txt
rm dualKS.txt
rm monomerLocs.txt
rm dimerLocs.txt
rm bin1.txt
#cd /data/users/robertbt/libsmolWEDimerCode/LibsmolWE/scripts
#./makePubs.pl $1bind${bindR}unbind${unbindK}runNo${runNo}
#cp pubs/$1bind${bindR}unbind${unbindK}runNo${runNo}.pub /dfs3/pub/robertbt/WELibsmolData/$1bind${bindR}unbind${unbindK}runNo${runNo}
#cd /dfs3/pub/robertbt/WELibsmolData/$1bind${bindR}unbind${unbindK}runNo${runNo}
qsub ./$1bind${bindR}unbind${unbindK}runNo${runNo}.pub

	cd /dfs3/pub/robertbt/WELibsmolData/
	printf "$1bind${bindR}unbind${unbindK}runNo${runNo}\n" >> missingData.txt
fi
fi

if [ -e $1bind${bindR}unbind${unbindK}runNo${runNo}.Time ]
then
if [ ! -e $1bind${bindR}unbind${unbindK}runNo${runNo}.Flux ]
then
rm $1bind${bindR}unbind${unbindK}runNo${runNo}.Time
#rm $1bind${bindR}unbind${unbindK}runNo${runNo}.Out
#rm $1bind${bindR}unbind${unbindK}runNo${runNo}.pub
rm ksOut.txt
rm dualKS.txt
rm monomerLocs.txt
rm dimerLocs.txt
rm bin1.txt
#cd /data/users/robertbt/libsmolWEDimerCode/LibsmolWE/scripts
#./makePubs.pl $1bind${bindR}unbind${unbindK}runNo${runNo}
#cp pubs/$1bind${bindR}unbind${unbindK}runNo${runNo}.pub /dfs3/pub/robertbt/WELibsmolData/$1bind${bindR}unbind${unbindK}runNo${runNo}
#cd /dfs3/pub/robertbt/WELibsmolData/$1bind${bindR}unbind${unbindK}runNo${runNo}
qsub ./$1bind${bindR}unbind${unbindK}runNo${runNo}.pub

	cd /dfs3/pub/robertbt/WELibsmolData/
	printf "$1bind${bindR}unbind${unbindK}runNo${runNo}\n" >> missingData.txt
fi
fi

done
done
done