SHELL=/bin/sh

INCLUDES= ../include/misclib.h ../include/fastadb.h \
          ../include/bioseq.h \
	  ../include/partition.h \
	  ../include/FSindex.h ../include/hit_list.h \
          ../include/randseq.h ../include/smatrix.h \
          ../include/pmatrix.h ../include/avl.h \
          ../include/keyword.h 

# ../include/blastlib.h ../include/clusters.h ../include/quadtree.h

OBJECTS= misclib.o fastadb.o  bioseq.o partition.o FSindex.o hit_list.o \
         randseq.o smatrix.o pmatrix.o \
	 sequential.o avl.o keyword.o

# blastlib.o quadtree.o clusters.o


VPATH = cephes-double

all: libFSindex.a libmd.a libexpat.la

libFSindex.a: $(OBJECTS) $(INCLUDES)
	ar crv libFSindex.a $(OBJECTS); \
	ranlib libFSindex.a; 

libmd.a:
	cd cephes-double; \
	make -f makefile libmd.a; \
	mv libmd.a ..;\
	cd ..; \

libexpat.la:
	cd expat-1.95.7; \
	./configure; \
	make buildlib; \
	cp libexpat.la ..; \
	cd ..; \

clean:
	rm -f *.o; \
	rm -f *.a; \
	rm -f *.la; \
	cd cephes-double; \
	make -f makefile clean; \
	cd expat-1.95.7; \
	make clean; \
	cd ..; \

srcclean: clean
	rm -f *~; 
