include ../make.inc

CXXFLAGS += -I. -I$(BLIS_ROOT)/include -I$(BLIS_ROOT)/include/blis
BLIS_LIB= -Wl,-rpath,$(BLIS_ROOT)/lib -L$(BLIS_ROOT)/lib -lblis
PFAFFINE_LIB= ../src/libpfaffine.a

TEXT_EXE= \
	./small_dskpfa_call.x \
	./large_dskpfa_call.x \
	./large_dskpfa_inv_call.x \
	./multi_dskpfa_call.x \
	./multi_dskpfa_inv_call.x \
	./small_zskpfa_call.x

all: $(TEXT_EXE)

clean:
	rm *.x
	rm *.o
	rm -r *.dSYM

.SUFFIXES: .x .cc

%.x : %.cc $(PFAFFINE_LIB)
	$(CXX) $(CXXFLAGS) $< $(PFAFFINE_LIB) $(BLAS_LIB) $(BLIS_LIB) -o $@
