include ../make.inc

CXXFLAGS += -I. -I$(BLIS_ROOT)/include -I$(BLIS_ROOT)/include/blis
BLIS_LIB= -L$(BLIS_ROOT)/lib -lblis
PFAFFINE_LIB= ../src/libpfaffine.a

TEXT_EXE= \
	./small_dskpfa_call.x \
	./large_dskpfa_call.x \
	./large_dskpfa_inv_call.x \
	./small_zskpfa_call.x

all: $(TEXT_EXE)

clean:
	rm -f *.x
	rm -f *.o
	rm -rf *.dSYM

.SUFFIXES: .x .cc

%.x : %.cc $(PFAFFINE_LIB)
	$(CXX) $(CXXFLAGS) $< $(PFAFFINE_LIB) $(BLIS_LIB) -o $@
