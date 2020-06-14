include ../make.inc

PFAFFINE_LIB= ../src/libpfaffine.a
TEXT_EXE= \
	./small_dskpfa_call.x \
	./large_dskpfa_call.x \
	./large_dskpfa_inv_call.x \
	./small_zskpfa_call.x \
	./small_dskr2k_call.x \
	./small_zskr2k_call.x
#	./small_dogemm_call.x

all: $(TEXT_EXE)

clean:
	rm *.x
	rm *.o
	rm -r *.dSYM

.SUFFIXES: .x .cc

%.x : %.cc
	$(CXX) $(CXXFLAGS) $< $(PFAFFINE_LIB) $(BLAS_LIB) -o $@
