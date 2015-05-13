# Pardiso.jl

The Pardiso.jl package provides an interface for using [PARDISO 5.0](http://www.pardiso-project.org/) from the [Julia language](http://julialang.org). You cannot use `Pardiso.jl` without having a valid licence for PARDISO. This package is available free of charge and in no way replaces or alters any functionality of PARDISO.

**Notes**: This package does currently not support MKL's PARDISO solver.


## Installation

`Pardiso.jl` expectes the following libraries to be loadable from within Julia with `dlopen`.
    * `libpardiso.so` - The PARDISO library.
    * `libgfortran.so` - The gfortran library. Should correspond to the same version as PARDISO is compiled against.
    * `libgomp.so` - Library for OpenMP

`Pardiso.jl`  has currently only been tested to work on Linux.


## Basic Usage

`Pardiso.jl` operates a bit like a state machine where you must first set the matrix type and subsequent calls to the solver routines will use the set type.


### Setting the matrix type

The matrix type (default 11) is set with `set_mtype(key)` where the key has the following meaning:
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

Note that currently, only real matrices is supported.

Four different versions are provided.

## More advanced usage.

To call the pardiso function directly

### Setting the solver
PARDISO also supports solving iteratively. The solver is set with `set_solver(key)` where the key has the following meaning:
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
Advanced users might want to explicitly set and retireve the `DPARM` and `IPARM` settings.
This can be done with the getters `get_iparm()`, `get_dparm()` and the setters `set_iparm(v::Int, i::Int)`, `set_dparm(v::FloatingPoint, i::Int)`, where the first argument is the value to set and the second is the index at which to set it.

To set the default values of `IPARM` and `DPARM` for a given matrix type and solver call `init_pardiso()`

When setting `IPARM` and `DPARM` explicitly, calls should now be made directly to
```
pardiso(X::VecOrMat{Float64}, A::SparseMatrixCSC, B::VecOrMat{Float64})
```
which will not modify the `IPARM` and `DPARM` values.

## Current limitations
    * Only Float64 matrices supported
    * No way to set a user defined fill-in reducing ordering.
 - There is currently no way of changing the following:


Documentation can be found in the [PARDISO manual](http://www.pardiso-project.org/manual/manual.pdf)
