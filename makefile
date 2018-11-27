all:
	gcc -Wall -O0 -g -c weightedEnsemble.c
	gcc weightedEnsemble.o -o WESim -lsmoldyn_shared