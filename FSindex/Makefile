SHELL=/bin/sh



SRCDIR = src/
PROGSDIR = src/progs/
LIBSDIR = src/lib/
BINDIR = bin/
LIBDIR = ../lib
INCLUDEDIR = src/include/

DEBUG=-g -D DEBUG=2
MPATROL = -D USE_MPATROL
MPATROLLIBS = -lmpatrol -lbfd -liberty
MAXOPT=-O4
GPROF = -pg
WARNINGS = -pipe -Wall -Wshadow -Wpointer-arith \
          -Wwrite-strings -Wstrict-prototypes \
          -Wformat -Wmissing-prototypes -funsigned-char \

CC = gcc

CFLAGS = "-I ../include $(WARNINGS) $(MAXOPT)"
LDFLAGS = "-L $(LIBDIR)"
LDLIBS = "-lm -lFSindex"


all: progs

progs: libs
	if [ ! -d  $(BINDIR) ]; then \
	   mkdir $(BINDIR); \
	fi 
	cd $(PROGSDIR); \
	gmake CC=$(CC) CFLAGS=$(CFLAGS) LDFLAGS=$(LDFLAGS) LDLIBS=$(LDLIBS); \
	for x in *; do if [ -x $$x -a -f $$x.c ]; \
	then mv $$x "../../$(BINDIR)"; fi done; \
	cd ..; \

libs:
	cd $(LIBSDIR); \
	gmake CC=$(CC) CFLAGS=$(CFLAGS) \

clean:
	rm -f core; \
	cd $(PROGSDIR); \
	gmake clean; \
	cd ../..; \
	cd $(LIBSDIR); \
	gmake clean; \
	cd ../..; \
	cd $(BINDIR); \
	rm -f core; \
	for x in *; do if [ -x $$x -a -f ../$(SRCDIR)$$x.c ]; \
	then rm -f $$x; fi done; \
	cd ..; \

srcclean:
	rm -f *~; \
	cd $(PROGSDIR); \
	gmake srcclean; \
	cd ../..; \
	cd $(LIBSDIR); \
	gmake srcclean; \
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



