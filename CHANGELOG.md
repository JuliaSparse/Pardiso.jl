# Changelog

## [1.2.0] - 2026-07-23

### Fixed

- `schur_complement(ps, M, x)` no longer corrupts the sparse input `x`
- `schur_complement` now restores the solver's iparms and phase after the call
  (also when the factorization errors) and validates its inputs (square matrix,
  `0 ≤ n < size(M, 1)`, valid transpose symbol, throwing `ArgumentError`) before
  touching any solver state; it is now restricted to `PardisoSolver` since
  MKL does not support it
- `pardisogetschur` allocates a large enough row-pointer buffer
- `printstats` uses the structure of the passed matrix instead of stale
  (possibly empty) internal buffers
- The release finalizer is registered once per solver object instead of on
  every `pardisoinit` call (previously every `solve!` accumulated finalizers)
- `isstructurallysymmetric` no longer throws a `BoundsError` for matrices
  with stored zeros
- The solution buffer `X` is now validated to be memory-contiguous, and
  phases that compute a solution reject empty dummy `X` buffers
- A failure while releasing memory in the `solve!` positive-definite
  fallback no longer masks the original error
- `set_nprocs!` for `MKLPardisoSolver` checks the return status of MKL
- `set_nprocs!` for `PardisoSolver` throws a descriptive error instead of a
  `MethodError`
- Removed broken internal `phases`/`valid_phases` helpers that referenced
  undefined constants
- `OMP_NUM_THREADS` values that do not parse as a single integer no longer
  throw when creating a `PardisoSolver`

## [1.1.2] - 2025-12-10

- fix Julia v1.13 support (#124)
- fix some typos / dead code (#123)

## [1.1.1] - 2025-11-18

- add missing function `_isnotzero` (#122)

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
