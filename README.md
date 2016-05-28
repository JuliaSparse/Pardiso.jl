# Pardiso.jl

The Pardiso.jl package provides an interface for using [PARDISO 5.0](http://www.pardiso-project.org/) and [Intel MKL PARDISO](https://software.intel.com/en-us/node/470282) from the [Julia language](http://julialang.org). You cannot use `Pardiso.jl` without either having a valid license for PARDISO or having the MKL library installed. This package is available free of charge and in no way replaces or alters any functionality of the linked libraries.

## Installation

The package itself is installed with `Pkg.add("Pardiso")` but you also need to follow the installation instructions below to install a working PARDISO library.

### MKL PARDISO

* Set the `MKLROOT` environment variable. See the [MKL getting started manual](https://software.intel.com/en-us/articles/intel-mkl-103-getting-started) for a thorough guide how to set this variable correctly.
* Run `Pkg.build("Pardiso")`

### PARDISO 5.0

#### Windows

* Put the PARDISO library `libpardiso500-WIN-X86-64.dll` in the `deps` folder.
* Run `Pkg.build("Pardiso")`

#### UNIX systems

* Put the PARDISO library `libpardiso500-GNUXXX-X86-64.so` in the `deps` folder.
* Install a (fast) installation of a BLAS and LAPACK (this should preferably be single threaded since PARDISO handles threading itself).
* Make sure OpenMP is installed.
* Make sure that the version of `gfortran` corresponding to the pardiso library is installed.
* Run `Pkg.build("Pardiso")``

## Basic Usage

This section will explain how to solve equations using `Pardiso.jl` with the default settings of the library. For more advanced usage there is a section further down.

## Creating the PardisoSolver

A `PardisoSolver` is created with `PardisoSolver()` for solving with PARDISO 5.0 or `MKLPardisoSolver()` for solving with MKL PARDISO. This object will hold the settings of the solver and will be passed into the solve functions. In the following sections an instance of a `PardisoSolver` or a `MKLPardisoSolver()` will be referred to as `ps`.

### Solving

Solving equations is done with the `solve` and `solve!` functions. They have the following signatures:

* `solve(ps, A, B)` solves `AX=B` and returns `X`
* `solve!(ps, X, A, B)` solves `AX=B` and stores it in `X`

The symbols `:T` or `:C` can be added as an extra argument to solve the transposed or the conjugate transposed system of equations, respectively.

Here is an example of solving a system of real equations with two right hand sides:

```jl
ps = PardisoSolver()

A = sparse(rand(10, 10))
B = rand(10, 2)
X = zeros(10, 2)
solve!(ps, X, A, B)
```

which happened to give the result

```jl
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

### Setting the number of threads

The number of threads to use are set in different ways for MKL PARDISO and PARDISO 5.0.

#### MKL PARDISO

```jl
set_nprocs!(ps, i) # Sets the number of threads to use
get_nprocs(ps) # Gets the number of threads being used
```

#### PARDISO 5.0

The number of threads are set at the creation of the `PardisoSolver` by looking for the environment variable `OMP_NUM_THREADS`. This can be done in Julia with for example `ENV["OMP_NUM_THREADS"] = 2`. **Note:** `OMP_NUM_THREADS` must be set *before* `Pardiso` is loaded and can not be changed during runtime.

The number of threads used by a `PardisoSolver` can be retrieved with `get_nprocs(ps)`

## More advanced usage.

This section discusses some more advanced usage of `Pardiso.jl`.

For terminology in this section please refer to the [PARDISO 5.0 manual](http://www.pardiso-project.org/manual/manual.pdf) and the [MKL PARDISO section](https://software.intel.com/en-us/node/470282).

After using functionality in this section, calls should no longer be made to the `solve` functions but instead directly to the function

```jl
pardiso(ps, X, A, B)
```

This will ensure that the properties you set will not be overwritten. 

If you want, you can use `get_matrix(ps, A, T)` to return a matrix that is suitable to use with `pardiso` depending on the matrix type that `ps` has set. The parameter `T` is a symbol representing if you will solve the normal, transposed or conjugated system. These are represented by `:N, :T, :C)` respectively.

For ease of use, `Pardiso.jl` provides enums for most options. These are not exported so has to either be explicitly imported or qualified with the module name first. It is possible to both use the enum as an input key to the options or the corresponding integer as given in the manuals.

### Setting the matrix type

The matrix type can be explicitly set with `set_matrixtype!(ps, key)` where the key has the following meaning:

| enum                 | integer | Matrix type                               |
|--------------------- |---------| ----------------------------------------  |
| REAL_SYM             | 1       | real and structurally symmetric           |
| REAL_SYM_POSDEF      | 2       | real and symmetric positive definite      |
| REAL_SYM_INDEF       | -2      | real and symmetric indefinite             |
| COMPLEX_STRUCT_SYM   | 3       | complex and structurally symmetric        |
| COMPLEX_HERM_POSDEF  | 4       | complex and Hermitian positive definite   |
| COMPLEX_HERM_INDEF   | -4      | complex and Hermitian indefinite          |
| COMPLEX_SYM          | 6       | complex and symmetric                     |
| REAL_NONSYM          | 11      | real and nonsymmetric                     |
| COMPLEX_NONSYM       | 13      | complex and nonsymmetric                  |

The matrix type for a solver can be retrieved with `get_matrixtype(ps)`.

### Setting the solver (5.0 only)
PARDISO 5.0 supports direct and iterative solvers. The solver is set with `set_solver!(ps, key)` where the key has the following meaning:

| enum               | integer | Solver                           |
|--------------------|---------|----------------------------------|
| DIRECT_SOLVER      | 0       | sparse direct solver             |
| ITERATIVE_SOLVER   | 1       | multi-recursive iterative solver |


### Setting the phase

Depending on the phase calls to `solve` (and `pardiso` which is mentioned later) does different things. The phase is set with `set_phase!(ps, key)` where key has the meaning:

| enum                                  | integer |  Solver Execution Steps                                         |
| --------------------------------------|---------|----------------------------------------------------------------|
| ANALYSIS                              | 11      | Analysis                                                       |
| ANALYSIS_NUM_FACT                     | 12      | Analysis, numerical factorization                              |
| ANALYSIS_NUM_FACT_SOLVE_REFINE        | 13      | Analysis, numerical factorization, solve, iterative refinement |
| NUM_FACT                              | 22      | Numerical factorization                                        |
| SELECTED_INVERSION                    | -22     | Selected Inversion                                             |
| NUM_FACT_SOLVE_REFINE                 | 23      | Numerical factorization, solve, iterative refinement           |
| SOLVE_ITERATIVE_REFINE                | 33      | Solve, iterative refinement                                    |
| SOLVE_ITERATIVE_REFINE_ONLY_FORWARD   | 331     | MKL only, like phase=33, but only forward substitution         |
| SOLVE_ITERATIVE_REFINE_ONLY_DIAG      | 332     | MKL only, like phase=33, but only diagonal substitution (if available) |
| SOLVE_ITERATIVE_REFINE_ONLY_BACKWARD  | 333     | MKL only, like phase=33, but only backward substitution
| RELEASE_LU_MNUM                       | 0       | Release internal memory for L and U matrix number MNUM         |
| RELEASE_ALL                           | -1      | Release all internal memory for all matrices                   |

### Setting `IPARM` and `DPARM` explicitly
Advanced users likely want to explicitly set and retrieve the `IPARM` and `DPARM` (5.0 only) parameters.
This can be done with the getters and setters:

```jl
get_iparm(ps, i) # Gets IPARM[i]
get_iparms(ps) # Gets IPARM
set_iparm!(ps, i, v) # Sets IPARM[i] = v

# 5.0 only
get_dparm(ps, i) # Gets DPARM[i]
get_dparms(ps) # Gets DPARM
set_dparm!(ps, i, v) # Sets DPARM[i] = v
```

To set the default values of the `IPARM` and `DPARM` call `pardisoinit(ps)`. The default values depend on what solver and matrix type is set.

### Setting message level

It is possible for Pardiso to print out timings and statistics when solving. This is done by `set_msglvl!(ps, key)` where key has the meaning:

| enum               | integer | Solver                           |
|--------------------|---------|----------------------------------|
| MESSAGE_LEVEL_OFF  | 0       | no statistics printed            |
| MESSAGE_LEVEL_ON   | 1       | statistics printed               |

In MKL PARDISO this is instead done by setting `IPARM[27]` to 1 before calling `pardiso`.

### Matrix and vector checkers

PARDISO 5.0 comes with a few matrix and vector checkers to check the consistency and integrity of the input data. These can be called with the functions:

```jl
printstats(ps, A, B)
checkmatrix(ps, A)
checkvec(ps, B)
```

### MNUM, MAXFCT, PERM

These are set and retrieved with the functions

```jl
set_mnum!(ps, i)
get_mnum(ps)

set_maxfct!(ps, i)
get_maxfct(ps)

get_perm(ps)
set_perm!(ps, perm) # Perm is a Vector{Int}
```

### Potential "gotchas"

* Julia uses CSC sparse matrices while PARDISO expects a CSR matrix. These can be seen as transposes of each other so to solve `AX = B` the transpose flag (`IPARAM[12]`) should be set to 1.
* For **symmetric** matrices, PARDISO needs to have the diagonal stored in the sparse structure even if the diagonal element happens to be 0. The manual recommends to add an `eps` to the diagonal when you suspect you might have 0 values diagonal elements that are not stored in the sparse structure.
* Unless `IPARM[1] = 1`, all values in `IPARM` will be ignored and default values are used.
* When solving a symmetric matrix, Pardiso expects only the upper triangular part. Since Julia has CSC matrices this means you should pass in `tril(A)` to the `pardiso` function. Use `checkmatrix` to see that you managed to get the matrix in a valid format.

# Contributions

If you have suggestions or idea of improving this package, please file an issue or even better, create a PR!
