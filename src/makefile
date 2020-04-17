include ../make.inc

OBJECTS= ../interface/pfaffine.o skpfa.o skr2k.o

include ../config/$(ARCH).mk

all: $(OBJECTS) libpfaffine.a

libpfaffine.a: $(OBJECTS)
	$(AR) rvu libpfaffine.a $(OBJECTS)

clean:
	rm -f *.o
	rm -f kernel/*.o
	rm -f ../interface/*.o
	rm -f libpfaffine.a

install: libpfaffine.a
	mkdir -p ../dest/include
	mkdir -p ../dest/lib
	cp skpfa.hh ../dest/include
	cp skr2k.hh ../dest/include
	cp libpfaffine.a ../dest/lib

.SUFFIXES: .o .cc

.cc.o:
	$(CXX) -I. $(CXXFLAGS) -c $< -o $@

