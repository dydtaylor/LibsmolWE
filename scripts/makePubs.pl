#!/usr/bin/perl

# create pub files to submit jobs to pbs

# make pubs ###########################################################

#loop over instance names, generated in makeparams2.pl


    #pub file name
    my $run_name=$ARGV[0];
	my $file_name=$run_name . ".pub";
	my $tempDir="/tmp/robertbt/${run_name}";
	my $saveDir="/dfs3/pub/robertbt/WELibsmolData/${run_name}";
    #open file for writing and print the following
    open(FOOD, ">pubs/$file_name" );
    print FOOD << "EOF";
#!/bin/bash
#\$ -N $run_name
#\$ -q rxn,free64,pub64
#\$ -ckpt restart
#\$ -r y
#\$ -e $run_name.err
#\$ -o $run_name.log
##Make temp folders for execution
mkdir /tmp/robertbt
mkdir $tempDir 

##Copy files from the save directory to the temp folder
cd $saveDir
cp weSmoldyn WEParams.txt dynamicsParams.txt ISEED $run_name.pub $tempDir

##Execute code
cd $tempDir
LD_LIBRARY_PATH=/data/users/robertbt/lib
export LD_LIBRARY_PATH
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
## Run executable
if [ -e savestate.txt ];
./weSmoldyn $run_name.Out $run_name.Flux $run_name.seed 0 $run_name.Time 1 &>/dev/null
else
./weSmoldyn $run_name.Out $run_name.Flux $run_name.seed 0 $run_name.Time 0 &>/dev/null
fi

##Move Output to the save directory
cp $run_name.Out $run_name.Flux $run_name.seed $run_name.Time savestate.txt ISEED mCountsWeighted.txt bin1.txt monomerLocs.txt dimerLocs.txt $saveDir
if [ -e ksOut.txt ];
then
cp ksOut.txt $saveDir
fi

if [ -e dualKS.txt ];
then
cp dualKS.txt $saveDir
fi

cd ..
##Cleanup temp directory
rm -rf tmp/robertbt/$run_name
echo Finished at `date`
EOF
    close FOOD;
    #done writing