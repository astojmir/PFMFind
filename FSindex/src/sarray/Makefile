SHELL=/bin/sh

INCLUDES=sarray.h

OBJECTS= sarray.o scode.o lcp.o ssarray.o 


all: libsarray.a

libsarray.a: libsarray.a($(OBJECTS)) $(INCLUDES)
	ranlib libsarray.a; 

clean:
	rm -f *.o; \
	rm -f *.a; \
	cd ..; \

srcclean: clean
	rm -f *~; 
