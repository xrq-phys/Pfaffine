include ../make.inc

OBJECTS= skpfa.o skr2k.o

all: $(OBJECTS) libpfaffine.a

libpfaffine.a: $(OBJECTS)
	$(AR) rvu libpfaffine.a $(OBJECTS)

clean:
	rm -f *.o
	rm -f libpfaffine.a

install: libpfaffine.a
	mkdir -p ../dest/include
	mkdir -p ../dest/lib
	cp skpfa.hh ../dest/include
	cp skr2k.hh ../dest/include
	cp libpfaffine.a ../dest/lib

.SUFFIXES: .o .cc

.cc.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

