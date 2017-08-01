# This is an example script demonstrating how PARDISO works on a small,
# sparse, real symmetric matrix. It computes the m solutions X to the
# collection of linear systems
#
#    A * X = B
#
# using the PARDISO solver, where A is a symmetric n x n matrix, B is an
# n x m matrix, and X is another n x m matrix.
using Pardiso
verbose = false

n = 4  # The number of equations.
m = 3  # The number of right-hand sides.

A = sparse([ 1. 0 -2  3
             0  5  1  2
            -2  1  4 -7
             3  2 -7  5 ])

# Generate a random collection of right-hand sides.
B = rand(n,m)

# Initialize the PARDISO internal data structures.
ps = PardisoSolver()

if verbose
    set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
end

# If we want, we could just solve the system right now.
# Pardiso.jl will automatically detect the correct matrix type,
# solve the system and free the data
X1 = solve(ps, A, B)

# We also show how to do this in incremental steps.

# First set the matrix type to handle general real symmetric matrices
set_matrixtype!(ps, Pardiso.REAL_SYM_INDEF)

# Initialize the default settings with the current matrix type
pardisoinit(ps)

# Get the correct matrix to be sent into the pardiso function.
# :N for normal matrix, :T for transpose, :C for conjugate
A_pardiso = get_matrix(ps, A, :N)

# Analyze the matrix and compute a symbolic factorization.
set_phase!(ps, Pardiso.ANALYSIS)
set_perm!(ps, randperm(n))
pardiso(ps, A_pardiso, B)
@printf("The factors have %d nonzero entries.\n", get_iparm(ps, 18))

# Compute the numeric factorization.
set_phase!(ps, Pardiso.NUM_FACT)
pardiso(ps, A_pardiso, B)
@printf("The matrix has %d positive and %d negative eigenvalues.\n",
        get_iparm(ps, 22), get_iparm(ps, 23))

# Compute the solutions X using the symbolic factorization.
set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)
X = similar(B) # Solution is stored in X
pardiso(ps, X, A_pardiso, B)
@printf("PARDISO performed %d iterative refinement steps.\n", get_iparm(ps, 7))

# Compute the residuals.
R = maximum(abs.(A*X - B))
@printf("The maximum residual for the solution X is %0.3g.\n", R)

# Free the PARDISO data structures.
set_phase!(ps, Pardiso.RELEASE_ALL)
pardiso(ps)
