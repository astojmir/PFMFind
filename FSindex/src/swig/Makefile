SHELL=/bin/sh

all: wrap

wrap: 	FS.i
	swig -python FS.i; \
	gcc -g -c FS_wrap.c -I ../include -I /usr/include/python2.2; \
	gcc -shared -L ../lib -L ../lib/expat-1.95.7 FS_wrap.o -lFSindex -lmd -lexpat -o _FS.so

clean:
	rm -f *.o; \
	rm -f *.a; \
	rm -f *.so;\
	rm -f *.pyc;\
	rm -f *.py;\
	rm -f *_wrap.c;

srcclean: clean
	rm -f *~; 


