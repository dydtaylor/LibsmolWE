#!/bin/bash

rho=8

cd ..

sed -i "12c\densityBit 1
" "dynamicsParams.txt"

sed -i "13c\density ${rho}
" "dynamicsParams.txt"

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

	for runNo in  $2
	do

		for nPart in  {30..100..5};
		do

			SAVEPATH=/pub/robertbt/WELibsmolData
			FULLNAME=$1.$d.ConstantDensity.rho${rho}.n${nPart}

			./makePubs.pl $FULLNAME ${runNo}
			cd ..

			mkdir $SAVEPATH/$FULLNAME

sed -i "8c\nPart	${nPart}
" "dynamicsParams.txt"
sed -i "5c\nBins ${nPart}
" "WEParams.txt"

			./seedchange ${runNo} 0
			cp weSmoldyn seedchange WEParams.txt dynamicsParams.txt binDefinitions1.txt binParams.txt corralsParams.txt ISEED scripts/pubs/$FULLNAME.pub $SAVEPATH/$FULLNAME

			cd scripts

		done
	done

	for runNo in  $2
	do

		for nPart in  {105..200..5};
		do

			SAVEPATH=/pub/robertbt/WELibsmolData
			FULLNAME=$1.$d.ConstantDensity.rho${rho}.n${nPart}

			./makePubs.pl $FULLNAME ${runNo}
			cd ..

			mkdir $SAVEPATH/$FULLNAME

sed -i "8c\nPart	${nPart}
" "dynamicsParams.txt"
sed -i "5c\nBins ${nPart}
" "WEParams.txt"

			./seedchange ${runNo} 0
			cp weSmoldyn WEParams.txt seedchange dynamicsParams.txt binDefinitions2.txt binParams.txt corralsParams.txt ISEED scripts/pubs/$FULLNAME.pub $SAVEPATH/$FULLNAME

		cd scripts

		done
	done

cd ..

sed -i "12c\densitybit 0
" "dynamicsParams.txt"

cd /dfs3/pub/robertbt/WELibsmolData


	for runNo in  $2;
	do
		for nPart in  {30..100..5};
		do
			FULLNAME=$1.$d.ConstantDensity.rho${rho}.n${nPart}
			cd $FULLNAME
			mv binDefinitions1.txt binDefinitions.txt
			qsub ./$FULLNAME.pub
			cd ..
		done
	done

	for runNo in  $2;
	do
		for nPart in  {105..200..5};
		do
			FULLNAME=$1.$d.ConstantDensity.rho${rho}.n${nPart}
			cd $FULLNAME
			mv binDefinitions2.txt binDefinitions.txt
			qsub ./$FULLNAME.pub
			cd ..

		done
	done