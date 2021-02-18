include ../make.inc
CONFIG= ../make.inc

ifeq ($(BLIS_ROOT),)
	BLIS_ROOT = ./dep
endif

OBJECTS= skpfa.o \
		 sktdf.o \
		 sktdi.o \
		 ../interface/pfaffine.o \
		 ../interface/pfapack.o
HEADERS= blalink.hh skpfa.hh
CXXFLAGS += -I. -I$(BLIS_ROOT)/include -I$(BLIS_ROOT)/include/blis
LDFLAGS += -L$(BLIS_ROOT)/lib -lblis

all: $(OBJECTS) libpfaffine.a libpfaffine.so

libpfaffine.so: $(OBJECTS)
	$(CXX) $(LDFLAGS) -shared $(OBJECTS) -o libpfaffine.so

libpfaffine.a: $(OBJECTS)
	$(AR) rvu libpfaffine.a $(OBJECTS)

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

