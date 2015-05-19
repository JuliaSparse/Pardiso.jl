# Pardiso.jl

The Pardiso.jl package provides an interface for using [PARDISO 5.0](http://www.pardiso-project.org/) and [Intel MKL PARDISO](https://software.intel.com/en-us/node/470282)from the [Julia language](http://julialang.org). You cannot use `Pardiso.jl` without either having a valid license for PARDISO or having the MKL library installed. This package is available free of charge and in no way replaces or alters any functionality of the linked libraries.

## Installation

### For MKL PARDISO

To use the MKL PARDISO the `MKLROOT` environment variable should be set. How to do this is shown [here](https://software.intel.com/en-us/articles/intel-mkl-103-getting-started).

### PARDISO 5.0

For PARDISO 5.0 the following libraries should be loadable from within Julia with `dlopen`.

* `libpardiso.so` - The PARDISO library.
* `libblas.so` - A (fast) BLAS library.
* `libgfortran.so` - The gfortran library. Should correspond to the same version as PARDISO is compiled against.
* `libgomp.so` - Library for OpenMP

**Note** The BLAS library should run in a single thread.

`Pardiso.jl`  has currently only been tested to work on Linux.


## Basic Usage

This section will explain how solve equations using `Pardiso.jl` with the default settings of the library. For more advanced usage there is a section further down.

## Creating the ParadisoSolver

A `ParadisoSolver` is created with `ParadisoSolver()` for solving with PARDISO 5.0 or `MKLPardisoSolver()` for solving with MKL PARDISO. This object will hold the settings of the solver and will be passed into the solve functions. In the following sections an instance of a `ParadisoSolver` or a `MKLPardisoSolver()`.will be referred to as `ps`.


### Setting the number of threads

The number of threads to use when solving is set in different ways for MKL PARDISO and PARDISO 5.0.

#### MKL PARDISO

```julia
set_nprocs(ps, i) # Sets the number of threads
get_nprocs(ps) # Gets the number of threads
```

#### PARDISO 5.0

The number of threads to use when solving is set at the creation of the `PardisoSolver` by looking for the environment variable `OMP_NUM_THREADS`. This can be done in Julia with `ENV["OMP_NUM_THREADS"] = 2`. If this variable does not exist, an exception is thrown.

The number of threads used by a `PardisoSolver` can be retrieved with `get_nprocs(ps)`

### Solving

Solving equations is done with the `solve` and `solve!` functions. They have the following signatures:

* `solve(ps, A, B)` solves `AX=B` and returns `X`
* `solve!(ps, X, A, B)` solves `AX=B` and stores it in `X`

If instead one wants to solve`A^T X = B`, the symbol `:T` should be passed as an extra last argument to the functions.

**Note**: The transposed versions are **not** the conjugate transpose in cases where `A` is complex.

Here is a contrived example of solving a system of real equations with two right hand sides:

```
ps = PardisoSolver()
set_mtype(ps, 11)

A = sparse(rand(10, 10))
B = rand(10, 2)
X = zeros(10, 2)
solve!(ps, X, A, B)
```

which happened to give the result

```julia
julia> X
10x2 Array{Float64,2}:
 -0.487361  -0.715372
 -0.644219  -3.38342
  0.465575   4.4838
  1.14448   -0.103854
  2.00892   -7.04965
  0.870507   1.7014
  0.590723  -5.74338
 -0.843841  -0.903796
 -0.279381   7.24754
 -1.17295    8.47922
```

## More advanced usage.

This section discusses some more advanced usage of `Pardiso.jl`

For terminology in this section please refer to the [PARDISO 5.0 manual](http://www.pardiso-project.org/manual/manual.pdf) and the [MKL PARDISO section](https://software.intel.com/en-us/node/470282).

### Setting the matrix type

The matrix type can be explicitly set with `set_mtype(ps, key)` where the key has the following meaning:

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

The matrix type for a solver can be retrieved with `get_mtype(ps)`.

### Setting the solver (5.0 only)
PARDISO 5.0 supports direct and iterative solvers. The solver is set with `set_solver(ps, key)` where the key has the following meaning:

| key | Solver                           |
|-----|----------------------------------|
| 0   | sparse direct solver             |
| 1   | multi-recursive iterative solver |


### Setting the phase

Depending on the phase calls to `solve` (and `pardiso` which is mentioned later) does different things. The phase is set with `set_phase(ps, key)` where key has the meaning:

| key   | Solver Execution Steps                                         |
|-------|----------------------------------------------------------------|
| 11    | Analysis                                                       |
| 12    | Analysis, numerical factorization                              |
| 13    | Analysis, numerical factorization, solve, iterative refinement |
| 22    | Numerical factorization                                        |
| -22   | Selected Inversion                                             |
| 23    | Numerical factorization, solve, iterative refinement           |
| 33    | Solve, iterative refinement                                    |
|331    | MKL only, like phase=33, but only forward substitution         |
|332    | MKL only, like phase=33, but only diagonal substitution (if available) |
|333    | MKL only,like phase=33, but only backward substitution
| 0     | Release internal memory for L and U matrix number MNUM         |
| -1    | Release all internal memory for all matrices                   |

### Setting `IPARM` and `DPARM` explicitly
Advanced users likely want to explicitly set and retrieve the `DPARM` (5.0 only) and `IPARM` settings.
This can be done with the getters and setters:

```julia
get_iparm(ps, i) # Gets IPARM[i]
get_iparms(ps) # Gets IPARM
set_iparm(ps, i, v) # Sets IPARM[i] = v

# 5.0 only
get_dparm(ps, i) # Gets DPARM[i]
get_dparms(ps) # Gets DPARM
set_dparm(ps, i, v) # Sets DPARM[i] = v
```

To set the default values of the `IPARM` and `DPARM` call `pardisoinit(ps)`. The default values depend on what solver and matrix type is set.

After setting `IPARM` and `DPARM` explicitly, calls should be made directly to the function
```
pardiso(ps, X, A, B)
```
which will not modify the `IPARM` and `DPARM` values.

### MNUM, MAXFCT, PERM

These are set and retrieved with the functions

```julia
set_mnum(ps, i)
get_mnum(ps)

set_maxfct(ps, i)
get_maxfct(ps)

get_perm(ps)
set_perm(ps, perm) # Perm is a Vector{Integer}
```

### PARDISO checkers (5.0 only)

PARDISO comes with a few matrix and vector checkers to check the consistency and integrity of the input data. These can be called with the functions:

```julia
printstats(ps, A, B)
checkmatrix(ps, A, B)
checkvec(B)
```

### Potential "gotchas"

* Julia uses CSC sparse matrices while PARDISO expects a CSR matrix. These can be seen as transposes of each other so to solve `AX = B` the transpose flag (`IPARAM[12]`) should be set to 1.
* For **symmetric** matrices, PARDISO needs to have the diagonal stored in the sparse structure even if the diagonal element happens to be 0. The manual recommends to add an `eps` to the diagonal when you suspect you might have 0 values diagonal elements that are not stored in the sparse structure.
* Unless `IPARM[1] = 1`, all values in `IPARM` will be ignored and default values are used.

# Contributions

If you have suggestions or idea of improving this package, please file an issue or even better, create a PR!
