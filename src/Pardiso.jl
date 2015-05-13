module Pardiso

using Base.LinAlg
using Base.SparseMatrix

export set_iparm, set_dparm, set_mtype, set_solver, set_phase, set_msglvl
export get_iparm, get_dparm, get_mtype, get_solver, get_phase, get_msglvl
export get_nprocs, pardiso
export check_matrix, check_vec, print_stats, init_pardiso
export pA_ldiv_B!, pA_ldiv_B, pAt_ldiv_B!, pAt_ldiv_B

# Libraries
const libblas= Libdl.dlopen("libblas", Libdl.RTLD_GLOBAL)
const libgfortran = Libdl.dlopen("libgfortran", Libdl.RTLD_GLOBAL)
const libgomp = Libdl.dlopen("libgomp", Libdl.RTLD_GLOBAL)
const libpardiso = Libdl.dlopen("libpardiso", Libdl.RTLD_GLOBAL)

# Pardiso functions
const init = Libdl.dlsym(libpardiso, "pardisoinit")
const pardiso_f = Libdl.dlsym(libpardiso, "pardiso")
const pardiso_chkmatrix = Libdl.dlsym(libpardiso, "pardiso_chkmatrix")
const pardiso_chkmatrix_z = Libdl.dlsym(libpardiso, "pardiso_chkmatrix_z")
const pardiso_printstats = Libdl.dlsym(libpardiso, "pardiso_printstats")
const pardiso_printstats_z = Libdl.dlsym(libpardiso, "pardiso_printstats_z")
const pardiso_chkvec = Libdl.dlsym(libpardiso, "pardiso_chkvec")
const pardiso_chkvec_z = Libdl.dlsym(libpardiso, "pardiso_chkvec_z")



# PARDISO states
const IPARM = zeros(Int32, 64)
# Set numper of processors to CPU_CORES unless "OMP_NUM_THREADS" is set
if ("OMP_NUM_THREADS" in keys(ENV))
    IPARM[3] = parse(Int, ENV["OMP_NUM_THREADS"])
else
     IPARM[3]= CPU_CORES
end
get_nprocs() = IPARM[3]


const DPARM = zeros(Float64, 64)
const PT = zeros(Int, 64)
const MTYPE = Int32[11]
const SOLVER = Int32[0]
const PHASE = Int32[33] # Default to solve and iterative refinement
const MSGLVL = Int32[0]


const VALID_MTYPES = [1, 2, -2, 3, 4, -4, 6, 11, 13]
const REAL_MTYPES = [1, 2, -2, 11]
const COMPLEX_MTYPES = [3, 4, -4, 6, 13]
const VALID_SOLVERS = [0, 1]
const VALID_PHASES = [11, 12, 13, 22, -22, 23, 33, 0, -1]
const VALID_MSGLVLS = [0, 1]

# Getters and setters
function set_solver(v::Int)
    v in VALID_SOLVERS || throw(ArgumentError(string("Invalid solver, valid solvers are 0 for",
                        " sparse direct solver, 1 for multi-recursive iterative solver")))
    SOLVER[1] = v
end
get_solver() = SOLVER[1]

function set_mtype(v::Int)
    v in VALID_MTYPES || throw(ArgumentError(string(
                                    "Invalid matrix type, valid matrix ",
                                    "types are $VALID_MTYPES.")))
    MTYPE[1] = v
end
get_mtype() = MTYPE[1]

set_iparm(v::Int, i::Int) = IPARM[i] = v
set_dparm(v::FloatingPoint, i::Int) = DPARM[i] = v

get_iparm() = IPARM
get_dparm() = DPARM

get_phase() = PHASE[1]
function set_phase(v::Int)
    v in VALID_PHASES || throw(ArgumentError(string(
                                    "Invalid phase, valid phases ",
                                    "are $VALID_PHASES.")))
    PHASE[1] = v
end

get_msglvl() = MSGLVL
function set_msglvl(v::Int)
    v in VALID_MSGLVLS || throw(ArgumentError(string(
                                "Invalid message level, valid message levels ",
                                "are $VALID_MSGLVLS.")))
    MSGLVL[1] = v
end


function init_pardiso()
    ERR = Int32[0]
    ccall(init, Void,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Float64}, Ptr{Int32}),
          PT, MTYPE, SOLVER, IPARM, DPARM, ERR)
    error_check(ERR)
    return 0
end
init_pardiso(MTYPE::Int, SOLVER::Int) = init_pardiso(convert(Int32, MTYPE),
                                                     convert(Int32, SOLVER))

# Solvers
function pA_ldiv_B{Ti, Tv}(A::SparseMatrixCSC{Tv, Ti}, B::VecOrMat{Tv})
    dim_check(B, A, B)
    X = similar(B)
    pA_ldiv_B!(X, A, B)
    return X
end

function pAt_ldiv_B{Ti, Tv}(A::SparseMatrixCSC{Tv, Ti}, B::VecOrMat{Tv})
    dim_check(B, A, B)
    X = similar(B)
    pAt_ldiv_B!(X, A, B)
    return X
end

function pA_ldiv_B!{Ti, Tv}(X::VecOrMat{Tv}, A::SparseMatrixCSC{Tv, Ti}, B::VecOrMat{Tv})
    dim_check(X, A, B)
    init_pardiso()
    IPARM[1] = 1
    IPARM[12] = 1
    pardiso(X, A, B)
    return X
end

function pAt_ldiv_B!{Ti, Tv}(X::VecOrMat{Tv}, A::SparseMatrixCSC{Tv, Ti}, B::VecOrMat{Tv})
    dim_check(X, A, B)
    init_pardiso()
    IPARM[1] = 1
    IPARM[12] = 0
    pardiso(X, A, B)
    return X
end



function pardiso{Ti, Tv}(X::VecOrMat{Tv}, A::SparseMatrixCSC{Tv, Ti}, B::VecOrMat{Tv})

    if Tv <: Complex && MTYPE[1] in REAL_MTYPES
        throw(ErrorException("Complex matrix and real matrix type set"))
    end

    if Tv <: Real && MTYPE[1] in COMPLEX_MTYPES
        throw(ErrorException("Real matrix and complex matrix type set"))
    end

    # For now only support one factorization
    MAXFCT = Int32(1)
    MNUM = Int32(1)

    N = Int32(size(A, 2))

    AA = A.nzval
    IA = convert(Vector{Int32}, A.colptr)
    JA = convert(Vector{Int32}, A.rowval)

    # For now don't support user fill-in reducing order
    PERM = Int32[]

    NRHS = Int32(size(B, 2))

    # For now disable messages
    MSGLVL = Int32(0)
    ERR = Int32[0]
    ccall(pardiso_f, Void,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Tv}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Tv},
           Ptr{Int32}, Ptr{Float64}),
          PT, &MAXFCT, &MNUM, MTYPE, PHASE,
          &N, AA, IA, JA, PERM,
          &NRHS, IPARM, &MSGLVL, B, X,
          ERR, DPARM)

    error_check(ERR)
    return
end

# Different checks
function print_stats{Ti, Tv}(X::VecOrMat{Tv}, A::SparseMatrixCSC{Tv, Ti},
                     B::VecOrMat{Tv})
    N = Int32(size(A, 2))
    AA = A.nzval
    IA = convert(Vector{Int32}, A.colptr)
    JA = convert(Vector{Int32}, A.rowval)
    NRHS = Int32(size(B, 2))
    ERR = Int32[0]
    if Tv <: Complex
        ccall(pardiso_printstats_z, Void,
              (Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Int32},
               Ptr{Int32}, Ptr{Int32}, Ptr{Tv},
               Ptr{Int32}),
              MTYPE, &N, AA, IA, JA, &NRHS, B, ERR)
    elseif Tv <: Real
            ccall(pardiso_printstats, Void,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}, Ptr{Tv},
           Ptr{Int32}),
          MTYPE, &N, AA, IA, JA, &NRHS, B, ERR)
    end
    error_check(ERR)
    return
end

function check_matrix{Ti, Tv}(A::SparseMatrixCSC{Tv, Ti}, B::VecOrMat{Tv})
    N = Int32(size(A, 1))
    AA = A.nzval
    IA = convert(Vector{Int32}, A.colptr)
    JA = convert(Vector{Int32}, A.rowval)
    ERR = Int32[0]

    if Tv <: Complex
        ccall(pardiso_chkmatrix_z, Void,
             (Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Int32},
              Ptr{Int32}, Ptr{Int32}),
             MTYPE, &N, AA, IA, JA, ERR)
    elseif Tv <: Real
        ccall(pardiso_chkmatrix, Void,
             (Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Int32},
              Ptr{Int32}, Ptr{Int32}),
             MTYPE, &N, AA, IA, JA, ERR)
    end
    error_check(ERR)
    return
end

function check_vec{Tv}(B::VecOrMat{Tv})
    N = Int32(size(B, 1))
    NRHS = Int32(size(B, 2))
    ERR = Int32[0]

    if Tv <: Complex
        ccall(pardiso_chkvec_z, Void,
             (Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Int32}),
             &N, &NRHS, B, ERR)
    elseif Tv <: Real
        ccall(pardiso_chkvec_z, Void,
             (Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Int32}),
             &N, &NRHS, B, ERR)
    end

    error_check(ERR)
    return
end

# Error checks
function dim_check(X, A, B)
    size(X) == size(B) || throw(DimensionMismatch(string(
                                 "Solution has $(size(X)), ",
                                 "RHS has size as $(size(B)).")))
    size(A,1) == size(B,1) || throw(DimensionMismatch(string(
                                    "Matrix has $(size(A,1)) ",
                                    "rows, RHS has $(size(B,1)) rows.")))
end


function error_check(err::Vector{Int32})
    err = err[1]
    if err == -1  ; error("Input inconsistent."); end
    if err == -2  ; error("Not enough memory."); end
    if err == -3  ; error("Reordering problem."); end
    if err == -4  ; error("Zero pivot, numerical fact. or iterative refinement problem."); end
    if err == -5  ; error("Unclassified (internal) error."); end
    if err == -6  ; error("Preordering failed (matrix types 11, 13 only)."); end
    if err == -7  ; error("Diagonal matrix problem."); end
    if err == -8  ; error("32-bit integer overflow problem."); end
    if err == -10 ; error("No license file pardiso.lic found."); end
    if err == -11 ; error("License is expired."); end
    if err == -12 ; error("Wrong username or hostname."); end
    if err == -100; error("Reached maximum number of Krylov-subspace iteration in iterative solver."); end
    if err == -101; error("No sufficient convergence in Krylov-subspace iteration within 25 iterations."); end
    if err == -102; error("Error in Krylov-subspace iteration."); end
    if err == -103; error("Break-Down in Krylov-subspace iteration."); end
    return
end

end # module
-