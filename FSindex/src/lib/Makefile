SHELL=/bin/sh

INCLUDES= ../include/misclib.h ../include/fastadb.h \
          ../include/bioseq.h \
	  ../include/partition.h \
	  ../include/FSindex.h ../include/hit_list.h \
          ../include/smatrix.h ../include/avl.h 

# ../include/blastlib.h ../include/clusters.h ../include/quadtree.h

OBJECTS= misclib.o fastadb.o  bioseq.o partition.o FSindex.o hit_list.o \
         smatrix.o avl.o 

#sequential.o randseq.o 

# blastlib.o quadtree.o clusters.o

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
