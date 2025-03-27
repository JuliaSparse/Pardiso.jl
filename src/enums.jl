# Matrix type
@enum(MatrixType::Int32,
    REAL_STRUCT_SYM     = 1,
    REAL_SYM_POSDEF     = 2,
    REAL_SYM_INDEF      = -2,
    COMPLEX_STRUCT_SYM  = 3,
    COMPLEX_HERM_POSDEF = 4,
    COMPLEX_HERM_INDEF  = -4,
    COMPLEX_SYM         = 6,
    REAL_NONSYM         = 11,
    COMPLEX_NONSYM      = 13,
)

Base.isreal(v::MatrixType) = v in (REAL_STRUCT_SYM, REAL_SYM_POSDEF, REAL_SYM_INDEF, REAL_NONSYM)
LinearAlgebra.issymmetric(v::MatrixType) = v in (REAL_SYM_POSDEF, REAL_SYM_INDEF,
                                        COMPLEX_HERM_POSDEF, COMPLEX_HERM_INDEF, COMPLEX_SYM)
LinearAlgebra.ishermitian(v::MatrixType) = v in (REAL_SYM_POSDEF, COMPLEX_HERM_POSDEF, COMPLEX_HERM_INDEF)
isposornegdef(v::MatrixType) = v in (REAL_SYM_POSDEF, REAL_SYM_INDEF, COMPLEX_HERM_POSDEF, COMPLEX_HERM_INDEF)

const MATRIX_STRING = Dict{MatrixType, String}(
    REAL_STRUCT_SYM => "Real structurally symmetric",
    REAL_SYM_POSDEF => "Real symmetric positive definite",
    REAL_SYM_INDEF => "Real symmetric indefinite",
    COMPLEX_STRUCT_SYM => "Complex structurally symmetric",
    COMPLEX_HERM_POSDEF => "Complex Hermitian postive definite",
    COMPLEX_HERM_INDEF => "Complex Hermitian indefinite",
    COMPLEX_SYM  => "Complex symmetric",
    REAL_NONSYM  => "Real nonsymmetric",
    COMPLEX_NONSYM  => "Complex nonsymmetric"
)

const REAL_MATRIX_TYPES = [REAL_STRUCT_SYM, REAL_SYM_POSDEF, REAL_SYM_INDEF, REAL_NONSYM]
const COMPLEX_MATRIX_TYPES = [COMPLEX_STRUCT_SYM, COMPLEX_HERM_POSDEF, COMPLEX_HERM_INDEF, COMPLEX_NONSYM, COMPLEX_SYM]

# Messages
@enum(MessageLevel::Int32,
    MESSAGE_LEVEL_OFF = 0,
    MESSAGE_LEVEL_ON = 1
)

# Solver
@enum(Solver::Int32,
    DIRECT_SOLVER = 0,
    ITERATIVE_SOLVER = 1
)

const SOLVER_STRING = Dict{Solver, String}(
    DIRECT_SOLVER => "Direct solver",
    ITERATIVE_SOLVER => "Iterative solver"
)

# Phase
@enum(Phase::Int32,
    ANALYSIS                             = 11,
    ANALYSIS_NUM_FACT                    = 12,
    ANALYSIS_NUM_FACT_SOLVE_REFINE       = 13,
    NUM_FACT                             = 22,
    SELECTED_INVERSION                   = -22,
    NUM_FACT_SOLVE_REFINE                = 23,
    SOLVE_ITERATIVE_REFINE               = 33,
    SOLVE_ITERATIVE_REFINE_ONLY_FORWARD  = 331,
    SOLVE_ITERATIVE_REFINE_ONLY_DIAG     = 332,
    SOLVE_ITERATIVE_REFINE_ONLY_BACKWARD = 333,
    RELEASE_LU_MNUM                      = 0,
    RELEASE_ALL                           = -1
)

const PHASE_STRING = Dict{Phase, String}(
    ANALYSIS                             => "Analysis",
    ANALYSIS_NUM_FACT                    => "Analysis, numerical factorization",
    ANALYSIS_NUM_FACT_SOLVE_REFINE       => "Analysis, numerical factorization, solve, iterative refinement",
    NUM_FACT                             => "Numerical factorization",
    SELECTED_INVERSION                   => "Selected Inversion",
    NUM_FACT_SOLVE_REFINE                => "Numerical factorization, solve, iterative refinement",
    SOLVE_ITERATIVE_REFINE               => "Solve, iterative refinement",
    RELEASE_LU_MNUM                      => "Release internal memory for L and U matrix number MNUM",
    RELEASE_ALL                           => "Release all internal memory for all matrices",
    SOLVE_ITERATIVE_REFINE_ONLY_FORWARD  => "like phase=33, but only forward substitution",
    SOLVE_ITERATIVE_REFINE_ONLY_DIAG     => "like phase=33, but only diagonal substitution (if available)",
    SOLVE_ITERATIVE_REFINE_ONLY_BACKWARD => "like phase=33, but only backward substitution",
)
