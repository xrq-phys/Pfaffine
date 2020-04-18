include ../make.inc

OBJECTS= skpfa.o skr2k.o \
		 ../interface/pfaffine.o \
		 ../interface/pfapack.o
HEADERS= blalink.hh kersel.hh skpfa.hh skr2k.hh
CXXFLAGS += -I.

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

.SUFFIXES: .o .cc .tcc

%.o : %.cc %.tcc $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

