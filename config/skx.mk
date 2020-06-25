# NB: This is no example configuration file but architecture epecification.
# Haswell's AVX2 kernels.

# CC +=
# C-flags needs additional define directive.
CFLAGS += -D_SkylakeX -march=skylake-avx512

# CXX +=
CXXFLAGS += -D_Haswell -D_SkylakeX

OBJECTS += kernel/mgemmn_haswell.o \
		   kernel/mgemmt_haswell.o \
		   kernel/dmgemm_skx.o

