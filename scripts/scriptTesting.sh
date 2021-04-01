#!/bin/bash

cd ..

make

#Date variable. 
d=`date +%m.%d.%Y`

runNo=$2

cd scripts

for n in 20
do

SAVEPATH=/pub/robertbt/WELibsmolData
FULLNAME=$1.$d.SimpleDiffusion.L3.n${n}

./makePubs.pl $FULLNAME $runNo
cd ..

mkdir $SAVEPATH/$FULLNAME

sed -i "8c\nPart	${n}
" "dynamicsParams.txt"


./seedchange $runNo 0


cp weSmoldyn WEParams.txt dynamicsParams.txt binParams.txt binDefinitions1.txt corralsParams.txt ISEED scripts/pubs/$FULLNAME.pub $SAVEPATH/$FULLNAME

./seedchange $runNo 1

cd scripts
done


cd /dfs3/pub/robertbt/WELibsmolData

for n in 20
do
FULLNAME=$1.$d.SimpleDiffusion.L3.n${n}
cd $FULLNAME
mv binDefinitions1.txt binDefinitions.txt
qsub ./$FULLNAME.pub
cd ..

done
