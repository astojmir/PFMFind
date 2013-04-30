SHELL=/bin/sh
VPATH = ../lib
all: _FS.so

FS_wrap.c: FS.i
	swig -python -modern -new_repr -noproxy -v FS.i;\
#	cp FS.py ../ShortFrags;

FS_wrap.o: FS_wrap.c
	$(CC) $(CFLAGS) -fPIC -c FS_wrap.c -I$(PYDIR)

_FS.so: FS_wrap.o libFSindex.a
	$(CC) -shared -g FS_wrap.o $(LDFLAGS) $(LDLIBS) -o FS.so

clean_wrap:
	rm -f *.pyc;\
	rm -f *.py;\
	rm -f *_wrap.c;\
	rm -f ../ShortFrags/FS.py

clean_lib:
	rm -f *.o; \
	rm -f *.a; \
	rm -f *.so;

clean: clean_wrap clean_lib

srcclean: clean
	rm -f *~;
