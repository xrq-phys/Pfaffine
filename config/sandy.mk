# NB: This is no example configuration file but architecture epecification.
# Sandy Bridge's SSE3 kernels.

# CC +=
CFLAGS += -march=sandybridge

# CXX +=
CXXFLAGS += -D_Sandy

OBJECTS += kernel/mgemm_sandy.o

