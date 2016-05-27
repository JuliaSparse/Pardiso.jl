# This is an example script demonstrating how PARDISO works on a
# medium-sized Hermitian positive definite matrix.
using Pardiso

# Script parameters.
# -----------------
verbose = false
n       = 100
lambda  = 3

# Create the Hermitian positive definite matrix A and the vector b in the
# linear system Ax = b.
e = ones(n,1)
e2 = ones(n-1,1)
A = spdiagm((im*e2, lambda*e, -im*e2), (-1,0,1))
b = rand(n,1) + im * zeros(n,1)

# Initialize the PARDISO internal data structures.
ps = PardisoSolver()

if verbose
    set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
end

# If we want, we could just solve the system right now.
# Pardiso.jl will automatically detect the correct matrix type,
# solve the system and free the data
X1 = solve(ps, A, b)

# First set the matrix type to handle general complex
# hermitian positive definite matrices
set_matrixtype!(ps, Pardiso.COMPLEX_HERM_POSDEF)

# Initialize the default settings with the current matrix type
pardisoinit(ps)

# Get the correct matrix to be sent into the pardiso function.
# :N for normal matrix, :T for transpose, :C for conjugate
A_pardiso = get_matrix(ps, A, :N)

# Analyze the matrix and compute a symbolic factorization.
set_phase!(ps, Pardiso.ANALYSIS)
pardiso(ps, A_pardiso, b)
@printf("The factors have %d nonzero entries.\n", get_iparm(ps, 18))

# Compute the numeric factorization.
set_phase!(ps, Pardiso.NUM_FACT)
pardiso(ps, A_pardiso, b)

# Compute the solutions X using the symbolic factorization.
set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)
x = similar(b) # Solution is stored in X
pardiso(ps, x, A_pardiso, b)
@printf("PARDISO performed %d iterative refinement steps.\n", get_iparm(ps, 7))

# Compute the residuals.
r = abs(A*x - b)
@printf("The maximum residual for the solution is %0.3g.\n",maximum(r))

# Free the PARDISO data structures.
set_phase!(ps, Pardiso.RELEASE_ALL)
pardiso(ps)
