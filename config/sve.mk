# NB: This is no example configuration file but architecture epecification.
# Aarch64 Cortex-A SVE kernels.
# Please use this kernel for the Post-K (Fugaku) supercomputer.

# CC +=
CFLAGS += -march=armv8-a+sve

# CXX +=
CXXFLAGS += -D_SVE

OBJECTS += kernel/mgemm_sve_ext.o \
           kernel/mgemm_sve_asm.o \
           kernel/vecln_sve.o \
           kernel/dvecln_sve.o \
           kernel/svecln_sve.o

