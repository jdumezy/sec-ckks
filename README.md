secCKKS - CKKS Functional bootstrapping - forked from OpenFHE
===================================

Somewhat Efficient Computations with CKKS - or secCKKS in short - is a fork of the [OpenFHE library](https://github.com/openfheorg/openfhe-development/) that implements functional bootstrapping for CKKS mostly following the blueprint given in Alexandru et al., [paper](https://eprint.iacr.org/2024/1623}) titled General Functional Bootstrapping using CKKS.
The original OpenFHE library provides FHE capabilities for BFV/BGV, CKKS and TFHE scheme, and I extended it by adding:

- Complex number support for CKKS

- Functional bootstrapping for CKKS for cleartexts up to 8 bits

The main benefit of functional bootstrapping is that it allows to evaluate arbitrary Look-Up Tables (LUTs) while bootstrapping and resetting the noises in ciphertexts to a nominal level. As such, it allows to perform a wider range of computations, including evaluating nice functions such as activation functions for neural networks.

This project is still in active development.

## Installation

OpenFHE should be compiled with 64 bits integer, as the 128 bits version is still under development.

You can refer to OpenFHE's guide on installation for your specific operating system:

- [Linux](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/linux.html)

- [MacOS](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/macos.html)

- [Windows](https://openfhe-development.readthedocs.io/en/latest/sphinx_rsts/intro/installation/windows.html)

Some unit test will fail with the current version of secCKKS (serialization).
The project was tested on Ubuntu 24.04.1 LTS and compiled with clang++.

## Usage

I highly recommend that you start by looking at the following example : [CKKS functional bootstrapping](src/pke/examples/ckks-functional-bootstrapping.cpp)

The example will be built automatically (if not set `BUILD_EXAMPLES` to `ON` with `set(BUILD_EXAMPLES ON)` in the main [CMakeLists.txt](CMakeLists.txt)) and the executable will be located in `build/bin/examples/pke/`, so you can run it from the project folder with:
```bash
./build/bin/examples/pke/ckks-functional-bootstrapping
```

The basic syntax of CKKS functional bootstrapping is the following: to evaluate an 8 bits LUT of the identity function with first order Hermite interpolation:
```c++
int p = pow(2, 4);
auto f = [](double x) -> double { return x; };
int hermite_order = 1;
auto ctxtResult = cc->EvalFuncBootstrap(ctxt, f, p, hermite_order);
```

Please pay attention to the fact that `EvalFuncBootstrapSetup` should be called instead of `EvalBootstrapSetup` to use functional bootstrapping.

For performances (especially multi-value bootstrapping), you should use the option:
```bash
export OMP_MAX_ACTIVE_LEVELS=4
```

## Detailed Changelog

- Fixed complex support for encoding

- Fixed complex constant/ciphertext addition/multiplication

- Added complex support for `EvalPoly` (linear)

- Added support for real->complex function interpolation via `EvalChebyshev` (linear and PS)

- Added Hermite coefficients calculation and `EvalHermiteFunction` for Hermite interpolation for up to order 3

- Added StC-first Bootstrapping (experimental)

- Added Functional Bootstrapping with use of both real and imaginary components

- Added Multi-Value Functional bootstrapping

- Added Tree base method for Functional Bootstrapping (WIP)

- Removed noise estimation and noise estimation unit tests

- Added the `FBTS` option to round at decryption

- Parallelization of (MV) Functional Bootstrapping with OpenMP

## License

This project is licensed under the BSD-2 License - see the [LICENSE](LICENSE) file for details.

This project was developed at CEA-List.

