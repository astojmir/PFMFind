SHELL=/bin/sh

INCLUDES= ../include/misclib.h ../include/fastadb.h ../include/partition.h \
	  ../include/FSindex.h ../include/hit_list.h  \
          ../include/randseq.h ../include/smatrix.h ../include/quadtree.h \
          ../include/pmatrix.h ../include/clusters.h ../include/avl.h \
          ../include/keyword.h ../include/bioseq.h
# ../include/blastlib.h
OBJECTS= misclib.o fastadb.o partition.o FSindex.o hit_list.o \
         randseq.o smatrix.o quadtree.o pmatrix.o clusters.o \
	 sequential.o avl.o keyword.o bioseq.o
# blastlib.o

all: libFSindex.a libmd.a

libFSindex.a: $(OBJECTS) $(INCLUDES)
	ar crv libFSindex.a $(OBJECTS); \
	ranlib libFSindex.a; 

libmd.a:
	cd cephes-double; \
	make -f makefile; \
	cp libmd.a ..;\
	cd ..; \

clean:
	rm -f *.o; \
	rm -f *.a; 

srcclean: clean
	rm -f *~; 


