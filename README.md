# Pardiso.jl

[![CI Testing](https://github.com/JuliaSparse/Pardiso.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaSparse/Pardiso.jl/actions/workflows/CI.yml)

The Pardiso.jl package provides an interface for using [Panua Pardiso](https://panua.ch/pardiso), it's predecessors from
[pardiso-project.org](http://www.pardiso-project.org/), and [Intel MKL
PARDISO](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2024-0/onemkl-pardiso-parallel-direct-sparse-solver-iface.html) from the [Julia
language](http://julialang.org).

You cannot use `Pardiso.jl` without either having a valid license for Panua Pardiso or
having the MKL library installed. This
package is available free of charge and in no way replaces or alters any
functionality of the linked libraries.

## Installation

The package itself is installed with `Pkg.add("Pardiso")` but you also need to
follow the installation instructions below to install a working PARDISO
library.

### MKL PARDISO

By default, when adding "Pardiso.jl" to the active environmnent, Julia will automatically install a suitable MKL for your platform by loading `MKL_jll.jl`.
Note that if you use a mac you will need to pin `MKL_jll` to version 2023.

If you instead use a self installed MKL, follow these instructions:

* Set the `MKLROOT` environment variable. See the [MKL set environment variables
  manual](https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-linux/2024-0/scripts-to-set-environment-variables.html)
  for a thorough guide how to set this variable correctly, typically done by
  executing something like `source /opt/intel/oneapi/setvars.sh intel64` or
  running `"C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\bin\mklvars.bat" intel64`
* Run `Pkg.build("Pardiso", verbose=true)`
* Eventually, run `Pardiso.show_build_log()` to see the build log for additional information.
* Note that the `MKLROOT` environment variable must be set, and `LD_LIBRARY_PATH` must contain `$MKLROOT/lib` whenever using the library this way.

### PARDISO from [panua.ch](https://panua.ch) ("PanuaPardiso", formerly "ProjectPardiso")

* Unzip the download file `panua-pardiso-yyyymmdd-os.zip` to some folder and set the environment variable `JULIA_PARDISO` to the `lib` subdirectory of this folder.  For example, create an entry `ENV["JULIA_PARDISO"] = "/Users/Someone/panua-pardiso-yyyymmdd-os/lib"` in `.julia/config/startup.jl`. If you have a valid license for the predecessor from pardiso-project.org, put the PARDISO library to a subdirectory denoted by `ENV["JULIA_PARDISO"]` and
  evenutally rename it to `libpardiso.so`.
* Perform the platform specific steps described below
* Run `Pkg.build("Pardiso", verbose=true)`
* Eventually, run `Pardiso.show_build_log()` to see the build log for additional information.

Note: In the past, weird errors and problems with MKL Pardiso had been observed when PanuaPardiso is enabled
(likely because some library that is needed by  PanauaPardiso was problematic with MKL).
In that case, if you want to use MKL Pardiso it is better to just disable  PanuaPardiso by not setting
the environment variable `JULIA_PARDISO` (and rerunning `Pkg.build("Pardiso")`).


## Basic Usage

This section will explain how to solve equations using `Pardiso.jl` with the default settings of the library. For more advanced users there is a section further down.

## Creating the PardisoSolver

A `PardisoSolver` is created with `PardisoSolver()` for solving with PanuaPardiso or `MKLPardisoSolver()` for solving with MKL PARDISO. This object will hold the settings of the solver and will be passed into the solve functions. In the following sections an instance of a `PardisoSolver` or an `MKLPardisoSolver()` will be referred to as `ps`.

### Solving

Solving equations is done with the `solve` and `solve!` functions. They have the following signatures:

* `solve(ps, A, B)` solves `AX=B` and returns `X`
* `solve!(ps, X, A, B)` solves `AX=B` and stores it in `X`

The symbols `:T` or `:C` can be added as an extra argument to solve the transposed or the conjugate transposed system of equations, respectively.

Here is an example of solving a system of real equations with two right-hand sides:

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

### Schur Complement (PanuaPardiso only)

Given a partitioned matrix `M = [A B; C D]`, the Schur complement of `A` in `M` is `S = D-CA⁻¹B`.
This can be found with the function `schur_complement` with the following signatures:

* `schur_complement(ps, M, n)` returns Schur complement of submatrix `A` in `M`, where `n` is the size of submatrix `D` (and therefore also of Schur complement)
* `schur_complement(ps, M, x)` returns Schur complement of submatrix `A` in `M`, where submatrix `D` is defined by nonzero rows of `SparseVector` or `SparseMatrix` `x`.

The symbols `:T` or `:C` can be added as an extra argument to solve the transposed or the conjugate transposed system of equations, respectively.

Here is an example of finding the Schur complement:

```jl
ps = PardisoSolver()
m = 100; n = 5; p = .5; T = Float64
rng = MersenneTwister(1234);
A = I + sprand(rng,T,m,m,p)
A⁻¹ = inv(Matrix(A))
B = sprand(rng,T,m,n,p)
C = sprand(rng,T,n,m,p)
D = sprand(rng,T,n,n,p)
M = [A B; C D]
S = schur_complement(ps,M,n)
```

which gives

```jl
julia> S
5×5 Array{Float64,2}:
  -0.121404    1.49473  -1.25965    7.40326    0.571538
 -19.4928     -7.71151  12.9496    -7.13646  -20.4194
   9.88029     3.35502  -7.2346     1.70651   13.9759
  -9.06094    -5.86454   7.44917   -2.54985   -9.17327
 -33.7006    -17.8323   20.2588   -19.5863   -37.6132
```

We can check the validity by comparing to explicity form:
```jl
julia> norm(D - C*A⁻¹*B - S)
5.033075778861378e-13
```

At present there seems to be an instability in the Schur complement computation for complex matrices.

### Setting the number of threads

The number of threads to use is set in different ways for MKL PARDISO and PanuaPardiso.

#### MKL PARDISO

```jl
set_nprocs!(ps, i) # Sets the number of threads to use
get_nprocs(ps) # Gets the number of threads being used
```

#### PanuaPardiso

The number of threads are set at the creation of the `PardisoSolver` by looking for the environment variable `OMP_NUM_THREADS`. This can be done in Julia with for example `ENV["OMP_NUM_THREADS"] = 2`. **Note:** `OMP_NUM_THREADS` must be set *before* `Pardiso` is loaded and can not be changed during runtime.

The number of threads used by a `PardisoSolver` can be retrieved with `get_nprocs(ps)`

## More advanced usage.

This section discusses some more advanced usage of `Pardiso.jl`.

For terminology in this section please refer to the [PanuaPardiso manual](http://panua.ch/manual/manual.pdf) and the [oneMKL PARDISO  manual](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2024-0/onemkl-pardiso-parallel-direct-sparse-solver-iface.html).

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
| REAL_STRUCT_SYM      | 1       | real and structurally symmetric           |
| REAL_SYM_POSDEF      | 2       | real and symmetric positive definite      |
| REAL_SYM_INDEF       | -2      | real and symmetric indefinite             |
| COMPLEX_STRUCT_SYM   | 3       | complex and structurally symmetric        |
| COMPLEX_HERM_POSDEF  | 4       | complex and Hermitian positive definite   |
| COMPLEX_HERM_INDEF   | -4      | complex and Hermitian indefinite          |
| COMPLEX_SYM          | 6       | complex and symmetric                     |
| REAL_NONSYM          | 11      | real and nonsymmetric                     |
| COMPLEX_NONSYM       | 13      | complex and nonsymmetric                  |

The matrix type for a solver can be retrieved with `get_matrixtype(ps)`.

### Setting the solver (PanuaPardiso only)
PanuatPardiso supports direct and iterative solvers. The solver is set with `set_solver!(ps, key)` where the key has the following meaning:

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
Advanced users likely want to explicitly set and retrieve the `IPARM` and `DPARM` (PanuaPardiso only) parameters.
This can be done with the getters and setters:

```jl
get_iparm(ps, i) # Gets IPARM[i]
get_iparms(ps) # Gets IPARM
set_iparm!(ps, i, v) # Sets IPARM[i] = v

# PanuaPardiso only
get_dparm(ps, i) # Gets DPARM[i]
get_dparms(ps) # Gets DPARM
set_dparm!(ps, i, v) # Sets DPARM[i] = v
```

To set the default values of the `IPARM` and `DPARM` call `pardisoinit(ps)`. The default values depend on what solver and matrix type is set.

### Setting message level

It is possible for Pardiso to print out timings and statistics when solving. This is done by `set_msglvl!(ps, key)` where `key` has the meaning:

| enum               | integer | Solver                           |
|--------------------|---------|----------------------------------|
| MESSAGE_LEVEL_OFF  | 0       | no statistics printed            |
| MESSAGE_LEVEL_ON   | 1       | statistics printed               |

### Matrix and vector checkers

PanuaPardiso comes with a few matrix and vector checkers to check the consistency and integrity of the input data. These can be called with the functions:

```jl
printstats(ps, A, B)
checkmatrix(ps, A)
checkvec(ps, B)
```

In MKL PARDISO this is instead done by setting `IPARM[27]` to 1 before calling `pardiso`.

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

### Schur Complement (PanuaPardiso only)

The `pardiso(ps,...)` syntax can be used to compute the Schur compelement (as described below). The answer can be retrieved with `pardisogetschur(ps)`.

To use the low-level API to compute the Schur complement:
  * use custom IPARMS (`set_iparm!(ps,1,1)`), set the Schur complement block size to `n` (`set_iparm!(ps,38,n)`), and set the phase to analyze & factorize (`set_phase!(ps,12)`).
  * compute the Schur complement by calling `pardiso(ps,X,M,X)`, where `B` is a dummy vector with `length(X)=size(M,1)` that shares element type with `M`.
  * retrieve with `pardisogetschur(ps)`


### Potential "gotchas"

* Julia uses CSC sparse matrices while PARDISO expects a CSR matrix. These can be seen as transposes of each other so to solve `AX = B` the transpose flag (`IPARAM[12]`) should be set to 1.
* For **symmetric** matrices, PARDISO needs to have the diagonal stored in the sparse structure even if the diagonal element happens to be 0. The manual recommends adding an `eps` to the diagonal when you suspect you might have 0 values diagonal elements that are not stored in the sparse structure.
* Unless `IPARM[1] = 1`, all values in `IPARM` will be ignored and default values are used.
* When solving a symmetric matrix, Pardiso expects only the upper triangular part. Since Julia has CSC matrices this means you should pass in `tril(A)` to the `pardiso` function. Use `checkmatrix` to see that you managed to get the matrix in a valid format.

# Contributions

If you have suggestions or idea of improving this package, please file an issue or even better, create a PR!
