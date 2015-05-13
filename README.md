# Pardiso.jl

The Pardiso.jl package provides an interface for using [PARDISO 5.0](http://www.pardiso-project.org/) from the [Julia language](http://julialang.org). You cannot use `Pardiso.jl` without having a valid license for PARDISO. This package is available free of charge and in no way replaces or alters any functionality of PARDISO.

**Notes**: This package does currently not support MKL's PARDISO solver.


## Installation

`Pardiso.jl` expects the following libraries to be loadable from within Julia with `dlopen`.

* `libpardiso.so` - The PARDISO library.
* `libgfortran.so` - The gfortran library. Should correspond to the same version as PARDISO is compiled against.
* `libgomp.so` - Library for OpenMP

`Pardiso.jl`  has currently only been tested to work on Linux.


## Basic Usage

This section will explain how to use solve equations using `Pardiso.jl` with the default settings of the library. For more advanced usage there is a section further down.

### Setting the matrix type

The matrix type (default 11) should be set before calling the solve functions. This is done with `set_mtype(key)` where the key has the following meaning:

| key   | Matrix type                               |
|----   |-----------------------------------------  |
| 1     | real and structurally symmetric           |
| 2     | real and symmetric positive definite      |
| -2    | real and symmetric indefinite             |
| 3     | complex and structurally symmetric        |
| 4     | complex and Hermitian positive definite   |
| -4    | complex and Hermitian indefinite          |
| 6     | complex and symmetric                     |
| 11    | real and nonsymmetric                     |
| 13    | complex and nonsymmetric                  |


### Setting the number of processors

The number of processors to use is set by defining the environment variable `OMP_NUM_THREADS` before loading the package. If this variable does not exist, the number of cores on the machine will be used.


### Solving

Four different versions are provided. The syntax now is not very beautiful but they are currently named this way to be similar to the Julia Base versions.

* `pA_ldiv_B!(X, A, B)` solves `AX = B` and stores the result in `X`.
* `pA_ldiv_B(A, B)` solves `AX = B` and returns a newly allocated `X`.
* `pAt_ldiv_B!(X, A, B)` solves `A^T X = B` and stores the result in `X`.
* `pAt_ldiv_B(A, B)` solves `A^T X = B` and returns a newly allocated `X`.


## More advanced usage.

For terminology in this section please refer to the [PARDISO manual](http://www.pardiso-project.org/manual/manual.pdf).

`Pardiso.jl` operates like a state machine where the properties of the solver is set before the call to the solve functions. After the solve function has completed, different types of data can be extracted.

### Setting the solver
PARDISO supports direct and iterative solvers. The solver is set with `set_solver(key)` where the key has the following meaning:

| key | Solver                           |
|-----|----------------------------------|
| 0   | sparse direct solver             |
| 1   | multi-recursive iterative solver |


### Setting the phase

Depending on the phase calls to `pardiso` does different things. The phase is set with `set_phase(key::Int)` where key has the meaning:

| key   | Solver Execution Steps                                         |
|-------|----------------------------------------------------------------|
| 11    | Analysis                                                       |
| 12    | Analysis, numerical factorization                              |
| 13    | Analysis, numerical factorization, solve, iterative refinement |
| 22    | Numerical factorization                                        |
| -22   | Selected Inversion                                             |
| 23    | Numerical factorization, solve, iterative refinement           |
| 33    | Solve, iterative refinement                                    |
| 0     | Release internal memory for L and U matrix number MNUM         |
| -1    | Release all internal memory for all matrices                   |

### Setting `IPARM` and `DPARM` explicitly
Advanced users might want to explicitly set and retrieve the `DPARM` and `IPARM` settings.
This can be done with the getters `get_iparm()`, `get_dparm()` and the setters `set_iparm(v::Int, i::Int)`, `set_dparm(v::FloatingPoint, i::Int)`, where the first argument is the value to set and the second is the index at which to set it.

To set the default values of `IPARM` and `DPARM` for a given matrix type and solver call `init_pardiso()`.

When setting `IPARM` and `DPARM` explicitly, calls should now be made directly to
```
pardiso(X, A, B)
```
which will not modify the `IPARM` and `DPARM` values.

# Contributions

If you have suggestions or idea of improving this package, please file an issue or even better, create a PR!
