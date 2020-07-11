☕️Pfaffine
==========

*Pfaffine aims to be the caffeine that wakes your program up when you work with Geminal wavefunctions.*

### Requirements

#### C++ Compiler
- Compiler supports C++ 11;
- Compiler recognizes `#pragma once`:
- [x] Almost all compilers has this support;
- GNU style inline assembly suport if custom kernel is to be used.
- [x] GCC 4.8 or later;
- [x] Clang and other LLVM-based compilers (e.g. ARM Allinea);
- [x] Intel Compilers with GCC/Clang ABI;
- [x] Fujitsu Compilers on Fugaku;
- [ ] Microsoft Visual C++.

#### Library and Others
- BLAS library;
- (Optional 1) Parallel BLAS library which is thread-safe under OpenMP;
- (Optional 2) OpenMP 3.0+ support of your C++ compiler.

### Build and Linking

Simple `makefile` is provided. One can copy the example `make.darwin` to `make.inc` and modify accordingly. `make install` builds and installs an instantiated library together with several templated headers into `{Repository root}/dest`.

Including directly the `.tcc` files is another way and it's more extensive in some sense.

### Build with OpenMP Threading Support

If one's computing environment satisfies the two optional requirements above, they might be able to build a thread-parallel version of Pfaffine library. To do this, one only needs to add OpenMP compiler options to `CXXFLAGS` in their `make.inc`. e.g. in cases of GCC or Clang, it's `-fopenmp` and in cases of Intel C++ Compiler it's `-qopenmp`.

## Usage

One can select between templated C++ and conventinoal C99 interfaces, which are stored `pfaffine.hh` and `pfaffine.h` correspondingly.

Interfaces are further categorized into 2 types: **full** interface that requires large memory allocation done outside function calls and **simplified** ones that allocates memory automatically for the user.

Full interface looks like this in C++/C99, correspondingly. Note that it can also be called from Fortran:

```cpp
template <typename T>
T skpfa(char uplo, unsigned n,
        T *A, unsigned ldA, unsigned inv,
        T *Sp1, T *Sp2, T *Sp3, T *Sp4, T *Sp5, unsigned npanel);
```

```c
void cal_?skpfa_(fpType *Pfa, char *uplo, unsigned *n, fpType *A, unsigned *ldA, unsigned *inv,
        fpType *Sp1, fpType *Sp2, fpType *Sp3, fpType *Sp4, fpType *Sp5, unsigned *npanel);
```

Here in C99 interface the `fpType` denotes floating point type corresponding to character `?` in function names, i.e. `double` for `dskpfa`, `double complex` for `zskpfa`, etc.

Following full interface is the **simplified** one, defined as:

```cpp
template <typename T>
T skpfa(char uplo, unsigned n, T *A, unsigned ldA, unsigned inv);
```

```c
void ?skpfa(fpType *Pfa, char uplo, unsigned n, fpType *A, unsigned ldA, unsigned inv);
```

In all cases above, `A` will be replaced by its full inverse if `inv=1`, otherwise the upper-half of tri-diagonal form will be kept. Detailed explanation for parameters referred here is available in docstring in `src/skpfa.tcc`.

## Roadmaps (WIP List)

- [x] **Important**: Change Fugaku's core to inline assembly & add clobber declaration;
- [ ] Mix AVX512 and AVX2 kernels on Intel SkylakeX processors;
- [x] Return signed Pfaffian instead of real-part-positive component;
- [x] Due to the new inversion method, some scratchpads in skpfa<T> is no longer used. Interface needs a clean-up;
- [ ] Implement a Pfaffian-inverse object that supports in-place *n*-term fast update;
- [x] Migrate `gemm` kernels from Pfapack and BLIS;
- [x] `needs improvement` SVE kernel for Fugaku;
- [x] `needs improvement` Default memory allocation;
- [x] Adjustable `k`-blocking;
- [x] Adjustable & automatic superblocking;
- [ ] Automatically determine panel size `npanel`;
- [x] Add a low-level **C99**/Fortran interface;
- [x] Add a compatibility interface consistent to [Pfapack](https://michaelwimmer.org/downloads.html);
- [x] Provide document here for all Pfaffine-defined interfaces;

### WIPs on Fugaku's SVE-512 and Other SVE Kernels

- [x] Complex kernels;
- [x] For real double, 14x16 (T-shaped) kernels instead of isotropic size (14x14).
- [x] Frame driver should provide support for boundary kernels.
- Optimize T-shapes other than SVE-512 for double-complex kernels.

## Licensing

Except for some [BSD-licensed](https://opensource.org/licenses/BSD-3-Clause) kernels from the [BLIS](https://github.com/flame/blis) project, these source code forms are provided under the [Mozilla Public License](https://www.mozilla.org/en-US/MPL). As long as one keep modifications of these specific files open, they are is allowed to use it for whatever they like *at her or his own risk*.
