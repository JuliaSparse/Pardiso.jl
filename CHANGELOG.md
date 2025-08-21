# Changelog

## [1.1.0] - 2025-08-21

- added `panua_is_loaded()` and `panua_is_licensed()` methods (#116)
- add finalizers to Pardiso structs (#117)

## [1.0.0] - 2025-01-24

### Breaking

- MKL v2025 dropped 32bit support. If 32bit support is needed, pin MKL to  v2024.

### Features

- Allow for MKL v2025
- Bump version to 1.0
- Remove superfluous loading of libblas and libgomp under linux
- add Changelog

## [0.5.6] - 2024-03-04

### Features

- Adaptations for Panua Pardiso (#75)

## [0.5.5] - 2024-02-24

### Features

- Bump MKL compat
- Try to pin MKL_jll to 2023 for macOS.
- introduces Pardiso.mkl_is_available()
- adaptations to use MKL Pardiso from locally installes oneAPI

## [0.5.4] - 2022-03-01

### Features

- Update to allow MKL 2022 

## [0.5.2] - 2021-07-09

### Features

- Allow StridedVecOrMat for RHS 
- Use MKL_jll if MKLROOT is not set
- Drop support for pardiso version 5
