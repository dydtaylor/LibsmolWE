#!/usr/bin/perl

# create pub files to submit jobs to pbs
#December 10, 2020, updating for use on HPC3 rather than HPC.
#This involves changing from SGE to SLURM

# make pubs ###########################################################

#loop over instance names, generated in makeparams2.pl


    #pub file name
    my $run_name_Prefix=$ARGV[0];
	my $file_name=$run_name_Prefix . ".pub";
	my $run_numbers=$ARGV[1];
	my $tempDirPrefix="/tmp/robertbt/${run_name_Prefix}";
	my $saveDirPrefix="/dfs6/pub/robertbt/WELibsmolData/${run_name_Prefix}";
    #open file for writing and print the following
    open(FOOD, ">pubs/$file_name" );
    print FOOD << "EOF";
#!/bin/bash
#SBATCH --job-name=$run_name_Prefix
#SBATCH -A elread_lab
#SBATCH -p standard
#SBATCH --output=$run_name_Prefix.%a.log
#SBATCH --error=$run_name_Prefix.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --array=1-$run_numbers

run_name=$run_name_Prefix.runNo\$SLURM_ARRAY_TASK_ID
tempDir=$tempDirPrefix.runNo\$SLURM_ARRAY_TASK_ID
saveDir=$saveDirPrefix
##If statement checks to make sure that the simulation has not previously finished and has not previously started.

	if [ ! -e \$run_name.Time ]
	then
	echo \$run_name.Time

	##Make temp folders for execution
		mkdir /tmp/robertbt
		mkdir \$tempDir 

##Copy files from the save directory to the temp folder
		cd \$saveDir
		cp weSmoldyn seedchange WEParams.txt corralsParams.txt dynamicsParams.txt binDefinitions.txt binParams.txt ISEED $run_name_Prefix.pub \$tempDir
		if [ -e savestate\$SLURM_ARRAY_TASK_ID.txt ];
		then
			cp \$run_name.Out \$run_name.Flux \$run_name.seed \$run_name.Time savestate\$SLURM_ARRAY_TASK_ID.txt ISEED mCountsWeighted.txt bin1.txt monomerLocs.txt dimerLocs.txt seedchange \$tempDir
		fi


##Execute code
		cd \$tempDir
		LD_LIBRARY_PATH=/data/homezvol2/robertbt/lib
		export LD_LIBRARY_PATH
		echo Running on host `hostname`
		echo Time is `date`
		echo Directory is `pwd`

##ISEED Changing
		./seedchange \$SLURM_ARRAY_TASK_ID 0

## Run executable
		if [ -e savestate\$SLURM_ARRAY_TASK_ID.txt ];
		then
			echo Savestate Found
			mv savestate\$SLURM_ARRAY_TASK_ID.txt savestate.txt
			./weSmoldyn \$run_name.Out \$run_name.Flux \$run_name.seed 0 \$run_name.Time 1 &>/dev/null
		fi

		if [ ! -e savestate\$SLURM_ARRAY_TASK_ID.txt ]
		then
			echo Savestate Not Found
			./weSmoldyn \$run_name.Out \$run_name.Flux \$run_name.seed 0 \$run_name.Time 0 &>/dev/null
		fi

##Move Output to the save directory
		mv savestate.txt savestate\$SLURM_ARRAY_TASK_ID.txt
		mv mCountsWeighted.txt mCountsWeighted\$SLURM_ARRAY_TASK_ID.txt
		cp \$run_name.Out \$run_name.Flux \$run_name.seed \$run_name.Time savestate\$SLURM_ARRAY_TASK_ID.txt ISEED mCountsWeighted\$SLURM_ARRAY_TASK_ID.txt bin1\$SLURM_ARRAY_TASK_ID.txt monomerLocs\$SLURM_ARRAY_TASK_ID.txt dimerLocs\$SLURM_ARRAY_TASK_ID.txt \$saveDir

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
		cd \$saveDir
		echo Finished at `date`
	fi

EOF
    close FOOD;
    #done writing