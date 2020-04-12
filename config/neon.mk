# NB: This is no example configuration file but architecture epecification.
# Aarch64 Cortex-A NEON kernels.

# CC +=
CFLAGS += -march=armv8-a

# CXX +=
CXXFLAGS += -D_Neon

OBJECTS += kernel/mgemm_neon.o

