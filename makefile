all:
	gcc -Wall -O0 -g -c weSmoldyn.c
	gcc weSmoldyn.o -L/usr/local/lib -I/System/Library/Frameworks/OpenGL.framework/Headers -I/System/Library/Frameworks/GLUT.framework/Headers -framework GLUT -framework OpenGL -framework Cocoa -L/System/Library/Frameworks/OpenGL.framework/Libraries -o weSmoldyn -lsmoldyn_shared -ltiff
jun:
	gcc -Wall -O0 -g -c weSmoldyn.c -I/usr/local/include
	gcc weSmoldyn.o -o weSmoldyn -L/Users/jun/Dropbox/science/library/smoldyn-2.54-mac/lib/ -lsmoldyn_shared
