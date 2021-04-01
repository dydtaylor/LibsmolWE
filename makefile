all:
	gcc -Wall -O0 -g -c weSmoldyn.c -I/data/users/robertbt/include
	gcc weSmoldyn.o -o weSmoldyn -L/data/users/robertbt/lib -lsmoldyn_shared -lm
	gcc -Wall -O0 -g -c seedchange.c -I/data/users/robertbt/include
	gcc seedchange.o -o seedchange -L/data/users/robertbt/lib
	gcc -Wall -O0 -g -c binMaker.c -I/data/users/robertbt/include -lm
	gcc binMaker.o -o binMaker -L/data/users/robertbt/lib -lm
jun:
	gcc -Wall -O0 -g -c weSmoldyn.c -I/usr/local/include
	gcc weSmoldyn.o -o weSmoldyn -L/Users/jun/Dropbox/science/library/smoldyn-2.54-mac/lib/ -lsmoldyn_shared
