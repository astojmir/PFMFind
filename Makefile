# Makefile to compile C libraries and programs

SHELL=/bin/sh
MAKE = make

SRCDIR = ./
BINDIR = ./bin

PROGSDIR = $(SRCDIR)/FSindex-testing/progs
LIBSDIR = $(SRCDIR)/lib
SWIGDIR = $(SRCDIR)/swig
INCLUDEDIR = $(SRCDIR)/include
INCS = -I $(abspath $(INCLUDEDIR))
LDFLAGS = -L $(abspath $(LIBSDIR))
BINSDIR = $(abspath $(BINDIR))

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
PYDIR = /usr/include/python2.7


# gcc flags
GPROF = -pg
CC = gcc
CWARNINGS = -pipe -Wall -Wunused-result -Wno-char-subscripts -Wshadow \
            -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables \
            -fasynchronous-unwind-tables \
          -Wwrite-strings -Wstrict-prototypes \
          -Wformat -Wmissing-prototypes -funsigned-char #-Werror
COPTIM=-O3 -fomit-frame-pointer -finline-functions -funroll-loops -D GCC_INLINE
LDLIBS = "-lFSindex -lgcc -lm -lpthread"


# VPATH = src/lib:src/swig

CFLAGS = $(COPTIM) $(INCS) -fPIC # -D THREADS=4

all: progs $(SWIGDIR)/FS.so

progs: $(LIBSDIR)/libFSindex.a
	if [ ! -d  $(BINDIR) ]; then \
	   mkdir $(BINDIR); \
	fi
	$(MAKE) -C $(PROGSDIR) CC=$(CC) CFLAGS="$(CFLAGS) $(CWARNINGS)" LDFLAGS=$(LDFLAGS) LDLIBS=$(LDLIBS); \
	for x in *; do if [ -x $$x -a -f $$x.c ]; \
	then mv $$x "../../$(BINDIR)"; fi done; \
	cd ..; \

$(SWIGDIR)/FS.so: $(LIBSDIR)/libFSindex.a
	$(MAKE) -C $(SWIGDIR) CC=$(CC) CFLAGS="$(CFLAGS)" LDFLAGS=$(LDFLAGS) LDLIBS=$(LDLIBS) PYDIR=$(PYDIR)

$(LIBSDIR)/libFSindex.a:
	$(MAKE) -C $(LIBSDIR) CC=$(CC) CFLAGS="$(CFLAGS) $(CWARNINGS)"

clean:
	$(MAKE) -C $(PROGSDIR) clean; \
	$(MAKE) -C $(LIBSDIR) clean; \
	$(MAKE) -C $(SWIGDIR) clean_lib; \
	rm -f $(BINSDIR)/*; \
