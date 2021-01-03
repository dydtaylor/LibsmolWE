#!/bin/bash

## Argv: 1 = runname, 2 = number of runs, 3 = custom date (if 0 then don't use custom date, if > 0 then use custom date) 4 = custom date string, 5 = worldL

#Date variable. 
cd ..
make
cd scripts
if [ 1 -gt $3 ];
then
d=`date +%m.%d.%Y`
fi

if [ $3 -gt 0 ];
then
d=$4
fi

runNo=$2
L=$5

for n in {2..36}
do

if [ $n -gt 25 ];
then
cd ..
./binMaker $n 1 $L
cd scripts
fi

SAVEPATH=/pub/robertbt/WELibsmolData
FULLNAME=$1.$d.SimpleDiffusion.L${L}.n${n}

./makePubs.pl $FULLNAME $runNo
cd ..

mkdir -p $SAVEPATH/$FULLNAME

sed -i "8c\nPart	${n}
" "dynamicsParams.txt"


./seedchange $runNo 0


cp seedchange weSmoldyn WEParams.txt dynamicsParams.txt binParams.txt binDefinitions.txt corralsParams.txt ISEED scripts/pubs/$FULLNAME.pub $SAVEPATH/$FULLNAME

./seedchange $runNo 1

cd scripts
done

for n in {38..54..2}
do
cd ..
./binMaker $n 1 $L
cd scripts

SAVEPATH=/pub/robertbt/WELibsmolData
FULLNAME=$1.$d.SimpleDiffusion.L${L}.n${n}

./makePubs.pl $FULLNAME $runNo
cd ..

mkdir -p $SAVEPATH/$FULLNAME

sed -i "8c\nPart	${n}
" "dynamicsParams.txt"


./seedchange $runNo 0


cp weSmoldyn seedchange WEParams.txt dynamicsParams.txt binParams.txt binDefinitions.txt corralsParams.txt ISEED scripts/pubs/$FULLNAME.pub $SAVEPATH/$FULLNAME

./seedchange $runNo 1

cd scripts
done

for n in {55..200..5} 128
do
cd ..
./binMaker $n 1 $L
cd scripts
SAVEPATH=/pub/robertbt/WELibsmolData
FULLNAME=$1.$d.SimpleDiffusion.L${L}.n${n}

./makePubs.pl $FULLNAME $runNo
cd ..

mkdir -p $SAVEPATH/$FULLNAME

sed -i "8c\nPart	${n}
" "dynamicsParams.txt"


./seedchange $runNo 0


cp weSmoldyn seedchange WEParams.txt dynamicsParams.txt binParams.txt binDefinitions.txt corralsParams.txt ISEED scripts/pubs/$FULLNAME.pub $SAVEPATH/$FULLNAME

./seedchange $runNo 1

cd scripts
done

cd /dfs6/pub/robertbt/WELibsmolData

for n in {2..36}
do
FULLNAME=$1.$d.SimpleDiffusion.L${L}.n${n}
cd $FULLNAME
qsub ./$FULLNAME.pub
cd ..

done

for n in {38..54..2}
do
FULLNAME=$1.$d.SimpleDiffusion.L${L}.n${n}
cd $FULLNAME
qsub ./$FULLNAME.pub
cd ..

done

for n in {55..200..5} 128
do
FULLNAME=$1.$d.SimpleDiffusion.L${L}.n${n}
cd $FULLNAME
qsub ./$FULLNAME.pub
cd ..

done