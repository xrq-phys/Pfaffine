include ../make.inc

OBJ = ltl2inv.o blalink_gemmt.o ilaenv_lauum.o ilaenv_wrap.o
SRC = ltl2inv.cc blalink_gemmt.c ilaenv_lauum.cc ilaenv_wrap.f90
HDR = invert.tcc pfaffian.tcc \
	  blalink.hh blalink_fort.h colmaj.hh
PFA2 = ../ext/pfapack2/fortran/libpfapack.a

all: lib dll

lib: libpfaffine.a

dll: libpfaffine.so

libpfaffine.a: $(OBJ)
	$(AR) rvu $@ $(OBJ)

libpfaffine.so: $(OBJ) $(PFA2)
	$(CXX) $(CXXFLAGS) -shared -o $@ $(OBJ) $(PFA2) $(LDFLAGS)

clean:
	rm -f $(OBJ)

SUFFIXES: .o .c

# -Wno-incompatible-pointer-types-discards-qualifiers is for BLAS declaration incompatibilities
.c.o: $(HDR)
	$(CC) $(OPTION) $(CFLAGS) -c $< -o $@ -Wno-incompatible-pointer-types-discards-qualifiers

SUFFIXES: .o .cc

.cc.o: $(HDR)
	$(CXX) $(CXXFLAGS) -D_CC_IMPL -c $< -o $@ -Wno-incompatible-pointer-types-discards-qualifiers

SUFFIXES: .o .f90

ilaenv_wrap.o: ilaenv_wrap.f90
	$(F90) $(FFLAGS) -c $< -o $@
