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

## Roadmaps (WIP List)

- ~~Migrate `gemm` kernels from Pfapack and BLIS;~~
- SVE kernel for Fugaku (Please use `neon` kernels at the moment);
- Default memory allocation;
- Adjustable `k`-blocking;
- Automatically determine panel size `npanel`;
- Add a low-level **C99**/Fortran interface;
- Support CMake.

## Licensing

Except for some [BSD-licensed](https://opensource.org/licenses/BSD-3-Clause) kernels from the [BLIS](https://github.com/flame/blis) project, these source code forms are provided under the [Mozilla Public License](https://www.mozilla.org/en-US/MPL). As long as one keeps modifications of these specific files open, she or he is allowed to use it for whatever she or he likes *at her or his own risk*.
