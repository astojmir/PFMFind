SHELL=/bin/sh



SRCDIR = src/
PROGSDIR = src/progs/
LIBSDIR = src/lib/
SWIGDIR = src/swig/
BINDIR = bin/
LIBDIR = ../lib
INCLUDEDIR = src/include/
GSLINCLUDE = -I /u/students/aleksand/exec/include/
GSLLIB = -L /u/students/aleksand/exec/lib

DEBUG=-g -D DEBUG=2
MPATROL = -D USE_MPATROL
MPATROLLIBS = -lmpatrol -lbfd -liberty
MAXOPT=-O3 -fomit-frame-pointer -finline-functions -funroll-loops -D GCC_INLINE
GPROF = -pg
WARNINGS = -pipe -Wall -Wno-char-subscripts -Wshadow \
          -Wwrite-strings -Wstrict-prototypes \
          -Wformat -Wmissing-prototypes -funsigned-char \

CC = gcc

CFLAGS = "-I ../include -Werror $(WARNINGS) $(DEBUG)"
LDFLAGS = "-L $(LIBDIR)"
LDLIBS = "-lm -lFSindex -lmd"


all: progs

progs: libs
	if [ ! -d  $(BINDIR) ]; then \
	   mkdir $(BINDIR); \
	fi 
#	cd $(PROGSDIR); \
#	make CC=$(CC) CFLAGS=$(CFLAGS) LDFLAGS=$(LDFLAGS) LDLIBS=$(LDLIBS); \
#	for x in *; do if [ -x $$x -a -f $$x.c ]; \
#	then mv $$x "../../$(BINDIR)"; fi done; \
#	cd ..; \

swig:
	cd $(SWIGDIR); \
	make; \
	cd ../..;\

libs:
	cd $(LIBSDIR); \
	make CC=$(CC) CFLAGS=$(CFLAGS) \

clean:
	rm -f core; \
	cd $(PROGSDIR); \
	make clean; \
	cd ../..; \
	cd $(LIBSDIR); \
	make clean; \
	cd ../..; \
	cd $(BINDIR); \
	rm -f core; \
	for x in *; do if [ -x $$x -a -f ../$(PROGSDIR)$$x.c ]; \
	then rm -f $$x; fi done; \
	cd ..; \

srcclean:
	rm -f *~; \
	cd $(PROGSDIR); \
	make srcclean; \
	cd ../..; \
	cd $(LIBSDIR); \
	make srcclean; \
	cd ../..; \
	cd $(INCLUDEDIR); \
	rm -f *~; \
	cd ../..; \

tarball: 
	tar -cvf FSindex.tar $(PROGSDIR)*.c; \
	tar -rvf FSindex.tar $(PROGSDIR)Makefile; \
	tar -rvf FSindex.tar $(LIBSDIR)*.c; \
	tar -rvf FSindex.tar $(INCLUDEDIR)*.h; \
	tar -rvf FSindex.tar $(LIBSDIR)Makefile; \
	tar -rvf FSindex.tar Makefile; \
	gzip  FSindex.tar; \



