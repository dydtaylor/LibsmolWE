#!/usr/bin/perl

# create pub files to submit jobs to pbs

# make pubs ###########################################################

#loop over instance names, generated in makeparams2.pl


    #pub file name
    my $run_name=$ARGV[0];
	my $file_name=$run_name . ".pub";
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
cd /data/users/robertbt/WELibsmolData/$run_name
LD_LIBRARY_PATH=/data/users/robertbt/lib
export LD_LIBRARY_PATH
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
# Run executable

./weSmoldyn $run_name.Out $run_name.Flux $run_name.seed 0 $run_name.Time &>/dev/null
echo Finished at `date`
EOF
    close FOOD;
    #done writing