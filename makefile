all:
	gcc -Wall -O0 -g -c MFPT.c
	gcc MFPT.o -o MFPT -lsmoldyn_shared