all:
	gcc -Wall -O0 -g -c weSmoldyn.c
	gcc weSmoldyn.o -o weSmoldyn -lsmoldyn_shared