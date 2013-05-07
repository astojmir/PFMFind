SHELL=/bin/sh

INCLUDES= ../include/misclib.h ../include/fastadb.h \
          ../include/bioseq.h \
	  ../include/partition.h \
	  ../include/FSindex.h ../include/hit_list.h \
          ../include/smatrix.h



OBJECTS= misclib.o fastadb.o  bioseq.o partition.o FSindex.o hit_list.o \
         smatrix.o 

all: libFSindex.a

libFSindex.a: libFSindex.a($(OBJECTS)) $(INCLUDES)
	ranlib libFSindex.a; 

clean:
	rm -f *.o; \
	rm -f *.a; \
	rm -f *.la; \
	cd ..; \

srcclean: clean
	rm -f *~; 
