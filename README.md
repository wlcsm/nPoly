# nPoly - Algebraic polynomials

A general polynomial library for the Rust programming language. Focusing on polynomial multiplication algorithms for arbitrary algebras.
This library served as the basis for a talk given at Maple Conf 2020 ["Rust for developing Fast, Parallelized Computer Algebra Systems"](https://www.youtube.com/watch?v=5ji7091sDi8&list=PLlcD7K2JXjTCSM5syc-A4qxd4_O4uV9D7&index=3)

Currently includes the following features:

- [x] Polynomial addition
- [x] Naive polynomial multiplication
- [x] Karatsuba Multiplication
- [x] FFT Multiplication
- [x] F4 algorithm for reducing a polynomial ideal to a Gröbner basis
- [x] Encoding/decoding polynomials as a human readable string

The following algebraic structures are supported:

- [x] Polynomial rings, including quotient rings and over multiple variables.
- [x] Polynomial ideals

Which are both supported over the real numbers, complex numbers, integers, and finite fields.

These operations are mostly supported on a dense polynomial representation but some operation are also supported on an experimental sparse polynomial representation.
