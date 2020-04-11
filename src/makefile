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

sandy: libpfaffine.a kernel/mgemm_sandy.o
	# Relink
	rm libpfaffine.a
	$(AR) rvu libpfaffine.a $(OBJECTS) kernel/mgemm_sandy.o

.SUFFIXES: .o .cc

.cc.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@


