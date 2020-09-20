include ../make.inc
CONFIG= ../make.inc

OBJECTS= skpfa.o \
		 ../interface/pfaffine.o \
		 ../interface/pfapack.o
HEADERS= blalink.hh skpfa.hh
CXXFLAGS += -I. -I$(BLIS_ROOT)/include -I$(BLIS_ROOT)/include/blis
LIBFLAGS += -L$(BLIS_ROOT)/lib
BLIS_STATIC = $(BLIS_ROOT)/lib/libblis.a

all: $(OBJECTS) libpfaffine.a

dll: $(OBJECTS) libpfaffine.so

libpfaffine.so: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@ -shared $(LIBFLAGS) -lblis

libpfaffine.a: $(OBJECTS)
	$(AR) rvu $@ $(OBJECTS)

clean:
	rm -f *.o
	rm -f ../interface/*.o
	rm -f libpfaffine.a

install: libpfaffine.a
	mkdir -p ../dest/include
	mkdir -p ../dest/lib
	cp skpfa.hh ../dest/include
	cp libpfaffine.a ../dest/lib

.SUFFIXES: .o .cc .tcc

%.o : %.cc %.tcc $(HEADERS) $(CONFIG)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.SUFFIXES: .o .c

%.o : %.c $(CONFIG)
	$(CC) $(CFLAGS) -c $< -o $@

.SUFFIXES: .o .s

%.o : %.s
	$(CC) $(CFLAGS) -c $< -o $@

