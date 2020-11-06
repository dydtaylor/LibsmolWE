#!/usr/bin/perl

# create pub files to submit jobs to pbs

# make pubs ###########################################################

#loop over instance names, generated in makeparams2.pl


    #pub file name
    my $run_name_Prefix=$ARGV[0];
	my $file_name=$run_name_Prefix . ".pub";
	my $run_numbers=$ARGV[1];
	my $tempDirPrefix="/tmp/robertbt/${run_name_Prefix}";
	my $saveDirPrefix="/dfs3/pub/robertbt/WELibsmolData/${run_name_Prefix}";
    #open file for writing and print the following
    open(FOOD, ">pubs/$file_name" );
    print FOOD << "EOF";
#!/bin/bash
#\$ -N $run_name_Prefix
#\$ -q rxn,free64,pub64
#\$ -ckpt restart
#\$ -r y
#\$ -e $run_name_Prefix.\$SGE_TASK_ID.err
#\$ -o $run_name_Prefix.\$SGE_TASK_ID.log
#\$ -t 1-$run_numbers

run_name=$run_name_Prefix.runNo\$SGE_TASK_ID
tempDir=$tempDirPrefix.runNo\$SGE_TASK_ID
saveDir=$saveDirPrefix
##First argv being > 0 means we are simply finding simulations that did not run and ONLY rerunning them. 

	if [ ! -e \$run_name.Time ]
	then
	echo \$run_name.Time

	##Make temp folders for execution
		mkdir /tmp/robertbt
		mkdir \$tempDir 

##Copy files from the save directory to the temp folder
		cd \$saveDir
		cp weSmoldyn seedchange WEParams.txt corralsParams.txt dynamicsParams.txt binDefinitions.txt binParams.txt ISEED $run_name_Prefix.pub \$tempDir
		if [ -e savestate\$SGE_TASK_ID.txt ];
		then
			cp \$run_name.Out \$run_name.Flux \$run_name.seed \$run_name.Time savestate\$SGE_TASK_ID.txt ISEED mCountsWeighted.txt bin1.txt monomerLocs.txt dimerLocs.txt seedchange \$tempDir
		fi


##Execute code
		cd \$tempDir
		LD_LIBRARY_PATH=/data/users/robertbt/lib
		export LD_LIBRARY_PATH
		echo Running on host `hostname`
		echo Time is `date`
		echo Directory is `pwd`

##ISEED Changing
		./seedchange \$SGE_TASK_ID 0

## Run executable
		if [ -e savestate\$SGE_TASK_ID.txt ];
		then
			echo Savestate Found
			mv savestate\$SGE_TASK_ID.txt savestate.txt
			./weSmoldyn \$run_name.Out \$run_name.Flux \$run_name.seed 0 \$run_name.Time 1 &>/dev/null
		fi

		if [ ! -e savestate\$SGE_TASK_ID.txt ]
		then
			echo Savestate Not Found
			./weSmoldyn \$run_name.Out \$run_name.Flux \$run_name.seed 0 \$run_name.Time 0 &>/dev/null
		fi


##Move Output to the save directory
		mv savestate.txt savestate\$SGE_TASK_ID.txt
		mv mCountsWeighted.txt mCountsWeighted\$SGE_TASK_ID.txt
		cp \$run_name.Out \$run_name.Flux \$run_name.seed \$run_name.Time savestate\$SGE_TASK_ID.txt ISEED mCountsWeighted\$SGE_TASK_ID.txt bin1.txt monomerLocs.txt dimerLocs.txt \$saveDir

		if [ -e ksOut.txt ];
		then
			cp ksOut.txt \$saveDir
		fi

		if [ -e dualKS.txt ];
		then
			cp dualKS.txt \$saveDir
		fi

		cd ..
	##Cleanup temp directory
		rm -rf /tmp/robertbt/\$run_name
		echo Finished at `date`
	fi

EOF
    close FOOD;
    #done writing