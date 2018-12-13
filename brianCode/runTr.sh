#!/bin/bash
#$ -N tm0 #job name
#$ -q rxn
#$ -pe openmp 1


#home=/dfs2/elread/rxn-share/chubk/cme #location of script files
folder=WE_main #script file names
#data=/fast-scratch/chubk/tm0
tau=.00003 #weighted ensemble timestep
NR=500 # number of replicas
typeflag=flux
inputflag=false
stopflag=false
WENitersMax=10000
module load MATLAB/r2014b

./WE_main $tau $NR $typeflag $inputflag $stopflag $WENitersMax
