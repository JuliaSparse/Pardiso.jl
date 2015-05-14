# Pardiso.jl

The Pardiso.jl package provides an interface for using [PARDISO 5.0](http://www.pardiso-project.org/) from the [Julia language](http://julialang.org). You cannot use `Pardiso.jl` without having a valid license for PARDISO. This package is available free of charge and in no way replaces or alters any functionality of PARDISO.

**Note**: This package does currently not support MKL's PARDISO solver.


## Installation

`Pardiso.jl` expects the following libraries to be loadable from within Julia with `dlopen`.

* `libpardiso.so` - The PARDISO library.
* `libgfortran.so` - The gfortran library. Should correspond to the same version as PARDISO is compiled against.
* `libgomp.so` - Library for OpenMP

`Pardiso.jl`  has currently only been tested to work on Linux.


## Basic Usage

This section will explain how solve equations using `Pardiso.jl` with the default settings of the library. For more advanced usage there is a section further down.

## Creating the ParadisoSolver

A `ParadisoSolver` is created with `ParadisoSolver()`. This object will hold the settings of the solver and will be passed into the solve functions. In the following sections an instance of a `ParadisoSolver` will be referred to as `ps` as if it was created like this:

```julia
julia> ps = PardisoSolver()
```

### Setting the matrix type

The matrix type should be set before calling the solve functions. This is done with `set_mtype(ps, key)` where the key has the following meaning:

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

### Setting the number of processors

The number of processors is set at the creation of the `PardisoSolver` by looking for the environment variable `OMP_NUM_THREADS`. This can be done in Julia with `ENV["OMP_NUM_THREADS"] = 2`. If this variable does not exist, the number of cores on the machine will be used.

The number of processors used by a solver can be retrieved with `get_nprocs(ps)`

### Solving

Solving equations is done with the `solve` and `solve!` functions. They have the following signatures:

* `solve(ps, A, B)` solves `AX=B` and returns X
* `solve!(ps, X, A, B)` solves `AX=B` and stores it in X

If instead one wants to solve`A^T X = B`, the symbol `:T` should be passed as an extra last argument to the functions.

**Note**: The transposed versions are **not** the conjugate transpose in cases where `A` is complex.

Here is a contrived example of solving a system of real equations with two right hand sides:

```
ps = PardisoSolver()
set_mtype(ps, 11)

A = sparse(rand(10, 10))
b = rand(10, 2)
x = zeros(10, 2)
solve!(ps, x, A, b)
````

which happened to give the result

```julia
julia> x
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

For terminology in this section please refer to the [PARDISO manual](http://www.pardiso-project.org/manual/manual.pdf).

### Setting the solver
PARDISO supports direct and iterative solvers. The solver is set with `set_solver(ps, key)` where the key has the following meaning:

| key | Solver                           |
|-----|----------------------------------|
| 0   | sparse direct solver             |
| 1   | multi-recursive iterative solver |


### Setting the phase

Depending on the phase calls to `pardiso` does different things. The phase is set with `set_phase(ps, key)` where key has the meaning:

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
This can be done with the getters `get_iparm(ps, i)`, `get_dparm(ps, i)` and the setters `set_iparm(ps, i, v)`, `set_dparm(ps, i, v)`, where the second argument is the value to set and the first is the index at which to set it.

To set the default values of the `IPARM` and `DPARM` states for a set state of matrix type and solver call `init_pardiso(ps)`.

When setting `IPARM` and `DPARM` explicitly, calls should now be made directly to
```
pardiso(ps, X, A, B)
```
which will not modify the `IPARM` and `DPARM` values.

Some potential "gotchas":

* Julia uses CSC sparse matrices while PARDISO expects a CSR matrix. These can be seen as transposes of each other so to solve `AX = B` the transpose flag (`IPARAM[12]`) should be set to 1.
* For symmetric matrices, PARDISO needs to have the diagonal stored in the sparse structure even if the diagonal element happens to be 0. The manual recommends to add an `eps` to the diagonal when you suspect you might have 0 values diagonal elements that are not stored in the sparse structure.

# Contributions

If you have suggestions or idea of improving this package, please file an issue or even better, create a PR!
