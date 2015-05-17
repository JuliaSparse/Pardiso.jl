module Pardiso

using Base.LinAlg
using Base.SparseMatrix

import Base.show

export PardisoSolver
export set_iparm, set_dparm, set_mtype, set_solver, set_phase, set_msglvl
export get_iparm, get_iparms, get_dparm, get_dparms
export get_mtype, get_solver, get_phase, get_msglvl, get_nprocs
export set_maxfct, set_perm, set_mnum
export get_maxfct, get_perm, get_mnum
export checkmatrix, checkvec, printstats, pardisoinit, pardiso
export solve, solve!

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

const VALID_MTYPES = [1, 2, -2, 3, 4, -4, 6, 11, 13]
const REAL_MTYPES = [1, 2, -2, 11]
const COMPLEX_MTYPES = [3, 4, -4, 6, 13]
const VALID_SOLVERS = [0, 1]
const VALID_PHASES= [11, 12, 13, 22, -22, 23, 33, 0, -1]
const VALID_MSGLVLS = [0, 1]

const SOLVERS = Dict{Int, ASCIIString}()
SOLVERS[0] = "Direct"
SOLVERS[1] = "Iterative"

const MTYPES = Dict{Int, ASCIIString}(
  1 => "Real structurally symmetric",
  2 => "Real symmetric positive definite",
 -2 => "Real symmetric indefinite",
  3 => "Complex structurally symmetric",
  4 => "Complex Hermitian postive definite",
 -4 => "Complex Hermitian indefinite",
  6 => "Complex symmetric",
 11 => "Real nonsymmetric",
 13 => "Complex nonsymmetric")

const PHASES = Dict{Int, ASCIIString}(
 11  => "Analysis",
 12  => "Analysis, numerical factorization",
 13  => "Analysis, numerical factorization, solve, iterative refinement",
 22  => "Numerical factorization",
-22  => "Selected Inversion",
 23  => "Numerical factorization, solve, iterative refinement",
 33  => "Solve, iterative refinement",
  0  => "Release internal memory for L and U matrix number MNUM",
 -1  => "Release all internal memory for all matrices")

typealias PardisoTypes Union(Float64, Complex128)

type PardisoSolver
    pt::Vector{Int}
    iparm:: Vector{Int32}
    dparm::Vector{Float64}
    mtype::Int32
    solver::Int32
    phase::Int32
    msglvl::Int32
    maxfct::Int32
    mnum::Int32
    perm::Vector{Int32}
end

function PardisoSolver()
    pt = zeros(Int, 64)
    iparm = zeros(Int32, 64)
    dparm = zeros(Float64, 64)
    mtype = 11 # Default to real unsymmetric matrices
    solver = 0 # Default to direct solver
    phase = 13 # Default to analysis + fact + solve + refine
    msglvl = 0

    # Set numper of processors to CPU_CORES unless "OMP_NUM_THREADS" is set
    if ("OMP_NUM_THREADS" in keys(ENV))
        iparm[3] = parse(Int, ENV["OMP_NUM_THREADS"])
    else
        iparm[3]= CPU_CORES
    end
    mnum = 1
    maxfct = 1
    perm = Int32[]
    PardisoSolver(pt, iparm, dparm, mtype, solver,
                  phase, msglvl, maxfct, mnum, perm)
end
show(io::IO, ps::PardisoSolver) = print(io, string("PardisoSolver:\n",
                                  "\tSolver: $(SOLVERS[get_solver(ps)])\n",
                                  "\tMatrix type: $(MTYPES[get_mtype(ps)])\n",
                                  "\tPhase: $(PHASES[get_phase(ps)])\n",
                                  "\tNum processors: $(get_nprocs(ps))"))
get_nprocs(ps::PardisoSolver) = ps.iparm[3]


# Getters and setters
function set_solver(ps::PardisoSolver, v::Integer)
    v in VALID_SOLVERS || throw(ArgumentError(string("invalid solver, valid solvers are 0 for",
                        " sparse direct solver, 1 for multi-recursive iterative solver")))
    ps.solver = v
end
get_solver(ps::PardisoSolver) = ps.solver

function set_mtype(ps::PardisoSolver, v::Integer)
    v in VALID_MTYPES || throw(ArgumentError(string(
                                    "invalid matrix type, valid matrix ",
                                    "types are $VALID_MTYPES")))
    ps.mtype = v
end
get_mtype(ps::PardisoSolver) = ps.mtype

get_iparm(ps::PardisoSolver, i::Integer) = ps.iparm[i]
get_iparms(ps::PardisoSolver) = ps.iparm
set_iparm(ps::PardisoSolver, i::Integer, v::Integer) = ps.iparm[i] = v

get_dparm(ps::PardisoSolver, i::Integer) = ps.dparm[i]
get_dparms(ps::PardisoSolver) = ps.dparm
set_dparm(ps::PardisoSolver, i::Integer, v::FloatingPoint) = ps.dparm[i] = v

get_mnum(ps::PardisoSolver) = ps.mnum
set_mnum(ps::PardisoSolver, mnum::Integer) = ps.mnum = mnum

get_maxfct(ps::PardisoSolver) = ps.maxfct
set_maxfct(ps::PardisoSolver, maxfct::Integer) = ps.maxfct = maxfct

get_perm(ps::PardisoSolver) = ps.perm
set_perm{T <: Integer}(ps::PardisoTypes, perm::Vector{T}) = ps.perm = convert(Vector{Int32}, perm)

get_phase(ps::PardisoSolver) = ps.phase

function set_phase(ps::PardisoSolver, v::Integer)
    v in VALID_PHASES|| throw(ArgumentError(string(
                                    "invalid phase, valid phases ",
                                    "are $VALID_PHASES")))
    ps.phase = v
end

get_msglvl(ps::PardisoSolver) = ps.msglvl
function set_msglvl(ps::PardisoSolver, v::Integer)
    v in VALID_MSGLVLS || throw(ArgumentError(string(
                                "invalid message level, valid message levels ",
                                "are $VALID_MSGLVLS")))
    ps.msglvl = v
end


function pardisoinit(ps::PardisoSolver)
    ERR = Int32[0]
    ccall(init, Void,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Float64}, Ptr{Int32}),
          ps.pt, &ps.mtype, &ps.solver, ps.iparm, ps.dparm, ERR)
    error_check(ERR)
    return
end

function solve{Ti, Tv <: PardisoTypes}(ps::PardisoSolver, A::SparseMatrixCSC{Tv, Ti},
                                       B::VecOrMat{Tv}, T::Symbol=:N)
  X = copy(B)
  solve!(ps, X, A, B, T)
  return X
end

function solve!{Ti, Tv <: PardisoTypes}(ps::PardisoSolver, X::VecOrMat{Tv},
                                        A::SparseMatrixCSC{Tv, Ti}, B::VecOrMat{Tv},
                                        T::Symbol=:N)
    pardisoinit(ps)

    if (T != :N && T != :T)
        throw(ArgumentError("only :T and :N are valid transpose symbols"))
    end

    # We need to set the transpose flag in PARDISO when we DON*T want
    # a transpose in Julia because we are passing a CSC formatted
    # matrix to PARDISO which expects a CSR matrix.
    if T == :N
      set_iparm(ps, 12, 1)
    end

    original_type = get_phase(ps)
    pardiso(ps, X, A, B)

    # Release memory, TODO: We are running the convert on IA and JA here
    # again which is unnecessary.
    set_phase(ps, -1)
    pardiso(ps, X, A, B)
    set_phase(ps, original_type)
    return X
end

function pardiso{Ti, Tv <: PardisoTypes}(ps::PardisoSolver, X::VecOrMat{Tv},
                                         A::SparseMatrixCSC{Tv, Ti}, B::VecOrMat{Tv})

    dim_check(X, A, B)

    if Tv <: Complex && get_mtype(ps) in REAL_MTYPES
        throw(ErrorException(string("input matrix is complex while PardisoSolver ",
                                    "has a real matrix type set")))
    end

    if Tv <: Real && get_mtype(ps) in COMPLEX_MTYPES
        throw(ErrorException(string("input matrix is real while PardisoSolver ",
                                    "has a complex matrix type set")))
    end

    N = Int32(size(A, 2))

    AA = A.nzval
    IA = convert(Vector{Int32}, A.colptr)
    JA = convert(Vector{Int32}, A.rowval)

    NRHS = Int32(size(B, 2))

    ERR = Int32[0]
    ccall(pardiso_f, Void,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Tv}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Tv},
           Ptr{Int32}, Ptr{Float64}),
          ps.pt, &ps.maxfct, &ps.mnum, &ps.mtype, &ps.phase,
          &N, AA, IA, JA, ps.perm,
          &NRHS, ps.iparm, &ps.msglvl, B, X,
          ERR, ps.dparm)

    error_check(ERR)
    # Return X here or not? For now, return.
    return X
end


# Different checks
function printstats{Ti, Tv <: PardisoTypes}(ps::PardisoSolver, A::SparseMatrixCSC{Tv, Ti},
                                            B::VecOrMat{Tv})
    N = Int32(size(A, 2))
    AA = A.nzval
    IA = convert(Vector{Int32}, A.colptr)
    JA = convert(Vector{Int32}, A.rowval)
    NRHS = Int32(size(B, 2))
    ERR = Int32[0]
    if Tv <: Complex
        f = pardiso_printstats_z
      else
        f = pardiso_printstats
    end
    ccall(f, Void,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}, Ptr{Tv},
           Ptr{Int32}),
          &ps.mtype, &N, AA, IA, JA, &NRHS, B, ERR)

    error_check(ERR)
    return
end

function checkmatrix{Ti, Tv <: PardisoTypes}(ps::PardisoSolver, A::SparseMatrixCSC{Tv, Ti})
    N = Int32(size(A, 1))
    AA = A.nzval
    IA = convert(Vector{Int32}, A.colptr)
    JA = convert(Vector{Int32}, A.rowval)
    ERR = Int32[0]

    if Tv <: Complex
        f = pardiso_chkmatrix_z
    else
        f = pardiso_chkmatrix
    end

    ccall(f, Void,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}),
          &ps.mtype, &N, AA, IA,
          JA, ERR)

    error_check(ERR)
    return
end

function checkvec{Tv <: PardisoTypes}(B::VecOrMat{Tv})
    N = Int32(size(B, 1))
    NRHS = Int32(size(B, 2))
    ERR = Int32[0]

    if Tv <: Complex
        f = pardiso_chkvec_z
    else
        f = pardiso_chkvec
    end
    ccall(f, Void,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Int32}),
          &N, &NRHS, B, ERR)

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
    err != -1  || throw(ErrorException("Input inconsistent."))
    err != -2  || throw(ErrorException("Not enough memory."))
    err != -3  || throw(ErrorException("Reordering problem."))
    err != -4  || throw(ErrorException("Zero pivot, numerical fact. or iterative refinement problem."))
    err != -5  || throw(ErrorException("Unclassified (internal) error."))
    err != -6  || throw(ErrorException("Preordering failed (matrix types 11, 13 only)."))
    err != -7  || throw(ErrorException("Diagonal matrix problem."))
    err != -8  || throw(ErrorException("32-bit integer overflow problem."))
    err != -10 || throw(ErrorException("No license file pardiso.lic found."))
    err != -11 || throw(ErrorException("License is expired."))
    err != -12 || throw(ErrorException("Wrong username or hostname."))
    err != -100|| throw(ErrorException("Reached maximum number of Krylov-subspace iteration in iterative solver."))
    err != -101|| throw(ErrorException("No sufficient convergence in Krylov-subspace iteration within 25 iterations."))
    err != -102|| throw(ErrorException("Error in Krylov-subspace iteration."))
    err != -103|| throw(ErrorException("Break-Down in Krylov-subspace iteration."))
    return
end

end # module

