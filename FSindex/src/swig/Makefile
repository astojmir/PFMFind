SHELL=/bin/sh

all: wrap

wrap: 	FS.i
	swig -python FS.i; \
	gcc -g -c FS_wrap.c -I ../include -I /usr/include/python2.2; \
	gcc -shared -L ../lib FS_wrap.o -lFSindex -lmd -o _FS.so

clean:
	rm -f *.o; \
	rm -f *.a; 

srcclean: clean
	rm -f *~; 


