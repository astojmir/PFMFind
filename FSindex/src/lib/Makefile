SHELL=/bin/sh
OBJECTS= misclib.o fastadb.o partition.o FSindex.o hit_list.o blastlib.o \
         randseq.o smatrix.o quadtree.o pmatrix.o\

all: lib

lib: $(OBJECTS)
	ar crv libFSindex.a $(OBJECTS); \
	ranlib libFSindex.a; \

clean:
	rm -f *.o; \
	rm -f *.a; \

srcclean: clean
	rm -f *~; \


