all:
	gcc -Wall -O0 -g -c weSmoldyn.c -I/data/users/robertbt/include
	gcc weSmoldyn.o -o weSmoldyn -L/data/users/robertbt/lib -lsmoldyn_shared -lm
	gcc -Wall -O0 -g -c seedchange.c -I/data/users/robertbt/include
	gcc seedchange.o -o seedchange -L/data/users/robertbt/lib
	gcc -Wall -O0 -g -c binMaker.c -I/data/users/robertbt/include -lm
	gcc binMaker.o -o binMaker -L/data/users/robertbt/lib -lm