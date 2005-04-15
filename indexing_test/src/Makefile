SHELL=/bin/sh

BINDIR = ../bin/
INCS = -I../../FSindex/src/include

MAKE = gmake
DEBUG=-g -D DEBUG=2
LDFLAGS = -L../../FSindex/src/lib -L../../FSindex/src/sarray

GPROF = -pg
CC = gcc
CWARNINGS = -pipe -Wall -Wno-char-subscripts -Wshadow \
          -Wwrite-strings -Wstrict-prototypes \
          -Wformat -Wmissing-prototypes -funsigned-char #-Werror
COPTIM=-O3 -fomit-frame-pointer -finline-functions -funroll-loops -D GCC_INLINE
LDLIBS = -lFSindex -lsarray -lm # -ldmalloc"

CFLAGS = $(COPTIM) $(INCS) $(CWARNINGS) -fPIC # -D USE_DMALLOC

all: progs 

progs: MVPT2_create MVPT_test MVPT_perf
	if [ ! -d  $(BINDIR) ]; then \
	   mkdir $(BINDIR); \
	fi
	for x in *; do if [ -x $$x -a -f $$x.c ]; \
	then mv $$x $(BINDIR); fi done; \

MVPT2_create: ufragdb.h mvptree.h avl.h ufragdb.o mvptree.o avl.o MVPT2_create.o
	$(CC) $(LDFLAGS) MVPT2_create.o ufragdb.o mvptree.o avl.o $(LDLIBS) -o MVPT2_create
	rm -rf MVPT2_create.o

MVPT_test: ufragdb.h mvptree.h avl.h ufragdb.o mvptree.o avl.o MVPT_test.o
	$(CC) $(LDFLAGS) MVPT_test.o ufragdb.o mvptree.o avl.o $(LDLIBS) -o MVPT_test
	rm -rf MVPT_test.o

MVPT_perf: ufragdb.h mvptree.h avl.h ufragdb.o mvptree.o avl.o MVPT_perf.o
	$(CC) $(LDFLAGS) MVPT_perf.o ufragdb.o mvptree.o avl.o $(LDLIBS) -o MVPT_perf
	rm -rf MVPT_perf.o

clean:
	rm -f *.o; \
	cd $(BINDIR); \
	rm -f core; \
	for x in *; do if [ -x $$x -a -f ../src/$$x.c ]; \
	then rm -f $$x; fi done; \
	cd ..; \

srcclean: clean
	rm -f *~; 

tarball: 
	tar -czvf index_test.tar.gz Makefile *.c *.h
