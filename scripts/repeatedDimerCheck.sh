#!/bin/bash
N=128

cd ..

sed -i "9c\reactBit	1
" "dynamicsParams.txt"

sed -i "5c\nBins 45
" "WEParams.txt"

sed -i "14c\monomerStart 0
" "dynamicsParams.txt"

make

#Date variable. 
if [ 1 -gt $3 ];
then
d=`date +%m.%d.%Y`
fi

if [ $3 -gt 0 ];
then
d=$4
fi

cd scripts


	sed -i "8c\nPart	${N}
	" "dynamicsParams.txt"

	for runNo in $2
	do

		for bindR in 0.0003 0.001 0.003 0.03  #$(seq .1 .1 1.0);
		do

			for unbindK in .001 .002 .005 .01 .02 .05 .1 .2 .5 1 2 5
			do


				SAVEPATH=/pub/robertbt/WELibsmolData
				FULLNAME=$1.$d.DimerMFPT.n${N}.bind$bindR.unbind$unbindK

				./makePubs.pl $FULLNAME ${runNo}
				cd ..

				#make

				mkdir $SAVEPATH/$FULLNAME

sed -i "6c\bindR	${bindR}
" "dynamicsParams.txt"

sed -i "7c\unbindK	${unbindK}
" "dynamicsParams.txt"

				./seedchange ${runNo} 0

				cp weSmoldyn WEParams.txt dynamicsParams.txt binDefinitions.txt corralsParams.txt binParams.txt ISEED seedchange scripts/pubs/$FULLNAME.pub $SAVEPATH/$FULLNAME

				cd scripts
			done
		done
	done


cd ..

sed -i "14c\monomerStart 1
" "dynamicsParams.txt"

cd scripts

for runNo in $2;
do

for bindR in 0.0003 0.001 0.003 0.03 #$(seq .1 .1 1.0);
do

for unbindK in 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000
do


SAVEPATH=/pub/robertbt/WELibsmolData
FULLNAME=$1.$d.DimerMFPT.n${N}.bind$bindR.unbind$unbindK

./makePubs.pl $FULLNAME ${runNo}
cd ..

#make

mkdir $SAVEPATH/$FULLNAME

sed -i "6c\bindR	${bindR}
" "dynamicsParams.txt"

sed -i "7c\unbindK	${unbindK}
" "dynamicsParams.txt"

./seedchange ${runNo} 0

cp weSmoldyn WEParams.txt dynamicsParams.txt binDefinitions.txt corralsParams.txt seedchange binParams.txt ISEED scripts/pubs/$FULLNAME.pub $SAVEPATH/$FULLNAME

cd scripts

done
done
done

cd ..
sed -i "9c\reactBit	0
" "dynamicsParams.txt"

cd /dfs3/pub/robertbt/WELibsmolData


for runNo in $2;
do

for bindR in 0.0003 0.001 0.003 0.03 #$(seq .1 .1 1.0);
do

for unbindK in .001 .002 .005 .01 .02 .05 .1 .2 .5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000
do

FULLNAME=$1.$d.DimerMFPT.n${N}.bind$bindR.unbind$unbindK

cd $SAVEPATH/$FULLNAME
qsub ./$FULLNAME.pub
cd ..
done
done
done
