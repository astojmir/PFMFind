SHELL=/bin/sh

BINDIR = ../bin/
INCS = -I../../FSindex/src/include

SADEPS= lcp.o sarray.o sarray.h scode.o ssarray.o

MAKE = gmake
DEBUG=-g -D DEBUG=2
LDFLAGS = -L../../FSindex/src/lib 

GPROF = -pg
CC = gcc
CWARNINGS = -pipe -Wall -Wno-char-subscripts -Wshadow \
          -Wwrite-strings -Wstrict-prototypes \
          -Wformat -Wmissing-prototypes -funsigned-char #-Werror
COPTIM=-O3 -fomit-frame-pointer -finline-functions -funroll-loops -D GCC_INLINE
LDLIBS = -lFSindex -lpthread -lm # -ldmalloc"

CFLAGS = $(DEBUG) $(INCS) $(CWARNINGS) -fPIC # -D USE_DMALLOC

all: progs 

progs: fragdist MVPT2_create MVPT_test MVPT_perf SAindex_create
	if [ ! -d  $(BINDIR) ]; then \
	   mkdir $(BINDIR); \
	fi
	for x in *; do if [ -x $$x -a -f $$x.c ]; \
	then mv $$x $(BINDIR); fi done; \

fragdist: avl.h avl.o fragdist.o
	$(CC) $(LDFLAGS) fragdist.o avl.o $(LDLIBS) -o fragdist
	rm -rf fragdist.o

MVPT2_create: ufragdb.h mvptree.h avl.h ufragdb.o mvptree.o avl.o MVPT2_create.o
	$(CC) $(LDFLAGS) MVPT2_create.o ufragdb.o mvptree.o avl.o $(LDLIBS) -o MVPT2_create
	rm -rf MVPT2_create.o

MVPT_test: ufragdb.h mvptree.h avl.h ufragdb.o mvptree.o avl.o MVPT_test.o
	$(CC) $(LDFLAGS) MVPT_test.o ufragdb.o mvptree.o avl.o $(LDLIBS) -o MVPT_test
	rm -rf MVPT_test.o

MVPT_perf: SAindex.h SAindex.o sarray.h lcp.o sarray.o ufragdb.h mvptree.h avl.h ufragdb.o mvptree.o avl.o MVPT_perf.o 
	$(CC) $(LDFLAGS) MVPT_perf.o SAindex.o lcp.o sarray.o mvptree.o ufragdb.o avl.o $(LDLIBS) -o MVPT_perf
	rm -rf MVPT_perf.o

clean:
	rm -f *.o; \
	cd $(BINDIR); \
	rm -f *; \
	cd ..;

SAindex_create: SAindex.h SAindex.o sarray.h lcp.o sarray.o SAindex_create.o
	$(CC) $(LDFLAGS) SAindex_create.o SAindex.o lcp.o sarray.o $(LDLIBS) -o SAindex_create

srcclean: clean
	rm -f *~

tarball: 
	tar -czvf index_test.tar.gz Makefile *.c *.h
