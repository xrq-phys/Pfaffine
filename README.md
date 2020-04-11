Pfaffine
========

*Pfaffine aims to be the caffeine that wakes your program up when you work with Geminal wavefunctions.*

### Requirements

- A working C++ compiler;
- BLAS library.

### Build and Linking

Simple `makefile` is provided. One can copy the example `make.darwin` to `make.inc` and modify accordingly. `make install` builds and installs an instantiated library together with several templated headers into `{Repository root}/dest`.

Including directly the `.tcc` files is another way and it's more extensive in some sense.

### Using

`skpfa<T>` (where `T` can be `float`, `double`, `complex<float>` or `complex<double>`) is available for doing a blocked-update based calculation of Pfaffian. One can find detailed explanation in docstring in `src/skpfa.tcc`.

## WIP List

- Migrate `gemm` kernels from Pfapack and BLIS;
- Add a low-level C/Fortran interface;
- Support CMake.

## Licensing

These source code forms are provided under the [Mozilla Public License](https://www.mozilla.org/en-US/MPL). As long as one keeps modifications of these specific files open, she or he is allowed to use it for whatever she or he likes *at her or his own risk*.