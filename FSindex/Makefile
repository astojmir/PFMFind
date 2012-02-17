# This makefile is very old and provisional
# a configure script would be very useful
SHELL=/bin/sh
MAKE = make

SRCDIR = src/
BINDIR = bin/
# All subdirs of $(SRCDIR) are relative to $(SRCDIR)
PROGSDIR = src/progs/
LIBSDIR = src/lib/
SWIGDIR = src/swig
SCRIPTDIR = src/scripts
LIBDIR = ../lib
INCLUDEDIR = src/include/
INCS = -I../include

# # mpatrol is a useful malloc debugger - commented out for now
# # http://mpatrol.sourceforge.net/
# MPATROL = -D USE_MPATROL
# MPATROLLIBS = -lmpatrol -lbfd -liberty

# # dmalloc is another such library
# # http://dmalloc.com/
# DMALLOC =  -D USE_DMALLOC
# DMALLOCLIBS = -ldmalloc

DEBUG=-g -D DEBUG=2

# This should be set to point to Python headers
PYDIR = /home/aleksand/exec/VirtualEnv/pfmfind-2.7/include/python2.7

LDFLAGS = "-L ../lib"

# gcc flags
GPROF = -pg
CC = gcc
CWARNINGS = -pipe -Wall -Wno-char-subscripts -Wshadow \
          -Wwrite-strings -Wstrict-prototypes \
          -Wformat -Wmissing-prototypes -funsigned-char #-Werror
COPTIM=-O3 -fomit-frame-pointer -finline-functions -funroll-loops -D GCC_INLINE
LDLIBS = "-lFSindex -lgcc -lm -lpthread"


# # Solaris cc flags
# GPROF = -xpg
# CC = cc
# CWARNINGS =  -errwarn -v
# COPTIM = -xO5 -fast -xipo -s
# LDLIBS = "-lFSindex"




VPATH = src/lib:src/swig

CFLAGS = $(DEBUG) $(INCS) -fPIC # -D THREADS=4

all: progs

progs: libFSindex swig
	if [ ! -d  $(BINDIR) ]; then \
	   mkdir $(BINDIR); \
	fi
	cd $(PROGSDIR); \
	$(MAKE) CC=$(CC) CFLAGS="$(CFLAGS) $(CWARNINGS)" LDFLAGS=$(LDFLAGS) LDLIBS=$(LDLIBS); \
	for x in *; do if [ -x $$x -a -f $$x.c ]; \
	then mv $$x "../../$(BINDIR)"; fi done; \
	cd ..; \

swig:	libFSindex
	cd $(SWIGDIR); \
	$(MAKE) CC=$(CC) CFLAGS="$(CFLAGS)" LDFLAGS=$(LDFLAGS) \
	LDLIBS=$(LDLIBS) PYDIR=$(PYDIR); \

libFSindex:
	cd $(LIBSDIR); \
	$(MAKE) CC=$(CC) CFLAGS="$(CFLAGS) $(CWARNINGS)"\

clean:
	rm -f core; \
	cd $(PROGSDIR); \
	$(MAKE) clean; \
	cd ../..; \
	cd $(LIBSDIR); \
	$(MAKE) clean; \
	cd ../..; \
	cd $(SWIGDIR); \
	$(MAKE) clean_lib; \
	cd ../..;\
	cd $(BINDIR); \
	rm -f core; \
	for x in *; do if [ -x $$x -a -f ../$(PROGSDIR)$$x.c ]; \
	then rm -f $$x; fi done; \
	rm -f *.py*; \
	cd ..; \

srcclean:
	rm -f *~; \
	cd $(PROGSDIR); \
	$(MAKE) srcclean; \
	cd ../..; \
	cd $(LIBSDIR); \
	$(MAKE) srcclean; \
	cd ../..; \
	cd $(INCLUDEDIR); \
	rm -f *~; \
	cd ../..; \
	cd $(SWIGDIR); \
	$(MAKE) srcclean; \
	cd ../..; \

tarball:
	tar -cvf FSindex.tar Makefile; \
	tar -rvf FSindex.tar data/*.*; \
	tar -rvf FSindex.tar --exclude=CVS matrix/*; \
	tar -rvf FSindex.tar $(PROGSDIR)*.c; \
	tar -rvf FSindex.tar $(PROGSDIR)Makefile; \
	tar -rvf FSindex.tar $(LIBSDIR)*.c; \
	tar -rvf FSindex.tar $(LIBSDIR)Makefile; \
	tar -rvf FSindex.tar $(INCLUDEDIR)*.h; \
	tar -rvf FSindex.tar Makefile; \
	tar -rvf FSindex.tar $(SWIGDIR)/*.i; \
	tar -rvf FSindex.tar $(SWIGDIR)/Makefile; \
	tar -rvf FSindex.tar $(SCRIPTDIR)/ShortFrags; \
	gzip  FSindex.tar; \



