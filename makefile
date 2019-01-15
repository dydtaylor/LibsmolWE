all:
	gcc -Wall -O0 -g -c weSmoldyn.c
	gcc weSmoldyn.o -o weSmoldyn -lsmoldyn_shared
jun:
	gcc -Wall -O0 -g -c weSmoldyn.c -I/usr/local/include
	gcc weSmoldyn.o -o weSmoldyn -L/Users/jun/Dropbox/science/library/smoldyn-2.54-mac/lib/ -lsmoldyn_shared
