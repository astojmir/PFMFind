SHELL=/bin/sh

# Uncomment if Sun C/C++ compiler is to be used.
#COMPILER= sun
#MACHINE = st-james
#MACHINE = morava

SRCDIR = src/
BINDIR = bin/
# All subdirs of $(SRCDIR) are relative to $(SRCDIR)
PROGSDIR = src/progs/
LIBSDIR = src/lib/
SARRAYDIR = src/sarray/
SWIGDIR = src/swig
SCRIPTDIR = src/scripts
LIBDIR = ../lib
INCLUDEDIR = src/include/
INCS = -I../include

GSLINCLUDE = -I /u/students/aleksand/exec/include/
GSLLIB = -L /u/students/aleksand/exec/lib
MPATROL = -D USE_MPATROL
MPATROLLIBS = -lmpatrol -lbfd -liberty
MAKE = gmake

DEBUG=-g -D DEBUG=2

ifeq ($(MACHINE), st-james)
PYDIR = /usr/pkg/include/python2p3
else
  ifeq ($(MACHINE), morava)
  PYDIR = /usr/include/python2.2
  else
  PYDIR = /usr/include/python2.3
  endif
endif

LDFLAGS = "-L ../lib -L ../sarray"

ifneq ($(COMPILER), sun)
# gcc flags 
GPROF = -pg
CC = gcc
CWARNINGS = -pipe -Wall -Wno-char-subscripts -Wshadow \
          -Wwrite-strings -Wstrict-prototypes \
          -Wformat -Wmissing-prototypes -funsigned-char #-Werror
COPTIM=-O3 -fomit-frame-pointer -finline-functions -funroll-loops -D GCC_INLINE
LDLIBS = "-lFSindex -lsarray -lgcc -lpthread"# -ldmalloc"
else
# Solaris cc flags 
GPROF = -xpg

CC = cc
CWARNINGS =  -errwarn -v
COPTIM = -xO5 -fast -xipo -s
LDLIBS = "-lFSindex -lsarray"
endif



VPATH = src/lib:src/swig

CFLAGS = $(DEBUG) $(INCS) -fPIC # -D USE_DMALLOC -D THREADS=4

all: progs

progs: libFSindex libsarray swig
	if [ ! -d  $(BINDIR) ]; then \
	   mkdir $(BINDIR); \
	fi 
	cd $(PROGSDIR); \
	$(MAKE) CC=$(CC) CFLAGS="$(CFLAGS) $(CWARNINGS)" LDFLAGS=$(LDFLAGS) LDLIBS=$(LDLIBS); \
	for x in *; do if [ -x $$x -a -f $$x.c ]; \
	then mv $$x "../../$(BINDIR)"; fi done; \
	cd ..; \

swig:	libFSindex libsarray
	cd $(SWIGDIR); \
	$(MAKE) CC=$(CC) CFLAGS="$(CFLAGS)" LDFLAGS=$(LDFLAGS) \
	LDLIBS=$(LDLIBS) PYDIR=$(PYDIR); \

libFSindex:
	cd $(LIBSDIR); \
	$(MAKE) CC=$(CC) CFLAGS="$(CFLAGS) $(CWARNINGS)"\

libsarray:
	cd $(SARRAYDIR); \
	$(MAKE) CC=$(CC) CFLAGS="$(CFLAGS)" \

clean:
	rm -f core; \
	cd $(PROGSDIR); \
	$(MAKE) clean; \
	cd ../..; \
	cd $(LIBSDIR); \
	$(MAKE) clean; \
	cd ../..; \
	cd $(SARRAYDIR); \
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
	tar -rvf FSindex.tar $(SARRAYDIR)*.c; \
	tar -rvf FSindex.tar $(SARRAYDIR)Makefile; \
	tar -rvf FSindex.tar $(SARRAYDIR)*.h; \
	tar -rvf FSindex.tar Makefile; \
	tar -rvf FSindex.tar $(SWIGDIR)/*.i; \
	tar -rvf FSindex.tar $(SWIGDIR)/Makefile; \
	tar -rvf FSindex.tar $(SCRIPTDIR)/ShortFrags; \
	gzip  FSindex.tar; \



