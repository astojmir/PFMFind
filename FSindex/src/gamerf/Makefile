SHELL=/bin/sh
OBJECTS = cdgamma.o  dcbrt.o derf.o derfc.o dgamma.o dierfc.o dlgamma.o


all: lib

lib: $(OBJECTS)
	ar rcs libgamerf.a *.o; \
	mv libgamerf.a ../lib/;

clean:
	rm -f *.o; \

