# AR executable.
AR=ar

# Compilers. 
# GCC with G++ is the general choice on linux.
# Intel compilers recommended on Intel processors.
CC=icc
CFLAGS=-O2
CXX=icpc
CXXFLAGS=-std=c++11 -O3 -qopenmp #-D_Use_OGemm #-D_Manual_SqBlk
#CXXFLAGS=-std=c++11 -O3 -D_Use_OGemm #-D_Pfaff_Debug

# Architecture specification.
# A list available under config/ directory.
# Use generic if one finds nothing suitable.
# One can create her or his own under config/ as well.
# A mac processor usually has at least SSE.
ARCH=haswell

# BLAS library linker flag.
# This is used only for building tests.
# An apparently good choice for mac is the Accelerate Framework.
#BLAS_LIB=-lmkl_rt
BLAS_LIB=-mkl=parallel
#BLAS_LIB=-mkl=sequential
