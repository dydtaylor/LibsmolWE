#!/bin/bash

#Date variable. 
if [ 1 -gt $3 ];
then
d=`date +%m.%d.%Y`
fi

if [ $3 -gt 0 ];
then
d=$4
fi

runNo=$2

for L in 2 2.25 2.5 2.75 3 3.25 3.5 3.75
do
for n in {1..30}
do

SAVEPATH=/pub/robertbt/WELibsmolData
FULLNAME=$1.$d.SimpleDiffusion.L${L}.n${n}

./makePubs.pl $FULLNAME $runNo
cd ..

mkdir $SAVEPATH/$FULLNAME

sed -i "2c\worldL	${L}
" "dynamicsParams.txt"

sed -i "8c\nPart	${n}
" "dynamicsParams.txt"


./seedchange $runNo 0


cp weSmoldyn WEParams.txt dynamicsParams.txt binParams.txt binDefinitions1.txt corralsParams.txt ISEED scripts/pubs/$FULLNAME.pub $SAVEPATH/$FULLNAME

./seedchange $runNo 1

cd scripts
done
done

cd /dfs3/pub/robertbt/WELibsmolData


for L in 2 2.25 2.5 2.75 3 3.25 3.5 3.75
do
for n in {1..30}
do
FULLNAME=$1.$d.SimpleDiffusion.L${L}.n${n}
cd $FULLNAME
mv binDefinitions1.txt binDefinitions.txt
qsub ./$FULLNAME.pub
cd ..

done
done
