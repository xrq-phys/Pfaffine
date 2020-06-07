Pfaffine
========

*Pfaffine aims to be the caffeine that wakes your program up when you work with Geminal wavefunctions.*

### Requirements

- A working C++ compiler;
- BLAS library.

### Build and Linking

Simple `makefile` is provided. One can copy the example `make.darwin` to `make.inc` and modify accordingly. `make install` builds and installs an instantiated library together with several templated headers into `{Repository root}/dest`.

Including directly the `.tcc` files is another way and it's more extensive in some sense.

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

- **Important**: Change Fugaku's core to inline assembly & add clobber declaration;
- Implement a Pfaffian-inverse object that supports in-place *n*-term fast update;
- ~~Migrate `gemm` kernels from Pfapack and BLIS;~~
- `needs improvement` ~~SVE kernel for Fugaku~~;
- `needs improvement` ~~Default memory allocation~~;
- Adjustable `k`-blocking;
- Automatically determine panel size `npanel`;
- ~~Add a low-level **C99**/Fortran interface~~;
- ~~Add a compatibility interface consistent to [Pfapack](https://michaelwimmer.org/downloads.html)~~;
- ~~Provide document here for all Pfaffine-defined interfaces;~~
- Support CMake.

### WIPs on Fugaku's SVE Kernel

- ~~Complex kernels~~;
- For real double, 14x16 (T-shaped) kernels in order that 16x14 can be fully utilized (currently 14x14 used).
- ~~Frame driver should provide support for boundary kernels.~~

## Licensing

Except for some [BSD-licensed](https://opensource.org/licenses/BSD-3-Clause) kernels from the [BLIS](https://github.com/flame/blis) project, these source code forms are provided under the [Mozilla Public License](https://www.mozilla.org/en-US/MPL). As long as one keeps modifications of these specific files open, she or he is allowed to use it for whatever she or he likes *at her or his own risk*.
