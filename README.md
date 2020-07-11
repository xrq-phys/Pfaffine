☕️Pfaffine
==========

*Pfaffine aims to be the caffeine that wakes your program up when you work with Geminal wavefunctions.*

Pfaffine is an numerical package for fast calculation of Pfaffian.
This branch has its background relinked against an extended version of [BLIS](https://github.com/flame/blis).

### Requirements

#### C++ Compiler
- Compiler supports C++ 11;
- Compiler recognizes `#pragma once`;

#### Library and Others
- My [BLIS library fork](https://github.com/xrq-phys/blis) with `skr2k`. This API-dependence is explicit;

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

## Licensing

Source code forms under this branch are provided under the [Mozilla Public License](https://www.mozilla.org/en-US/MPL). As long as one keep modifications of these specific files open, they are is allowed to use it for whatever they like *at her or his own risk*.
