SHELL=/bin/sh

INCLUDES= ../include/misclib.h ../include/fastadb.h ../include/partition.h \
	  ../include/FSindex.h ../include/hit_list.h ../include/blastlib.h \
          ../include/randseq.h ../include/smatrix.h ../include/quadtree.h \
          ../include/pmatrix.h ../include/clusters.h 

OBJECTS= misclib.o fastadb.o partition.o FSindex.o hit_list.o blastlib.o \
         randseq.o smatrix.o quadtree.o pmatrix.o clusters.o \
	 sequential.o

all: lib

lib: $(OBJECTS) $(INCLUDES)
	ar crv libFSindex.a $(OBJECTS); \
	ranlib libFSindex.a; 

clean:
	rm -f *.o; \
	rm -f *.a; 

srcclean: clean
	rm -f *~; 


