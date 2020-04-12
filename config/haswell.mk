# NB: This is no example configuration file but architecture epecification.
# Haswell's AVX2 kernels.

# CC +=
CFLAGS += -march=haswell

# CXX +=
CXXFLAGS += -D_Haswell

OBJECTS += kernel/mgemmn_haswell.o kernel/mgemmt_haswell.o

