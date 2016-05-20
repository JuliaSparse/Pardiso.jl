module Pardiso

using Compat

import Compat.String

if VERSION < v"0.5.0-dev+2915"
    import Compat.issymmetric
else
    import Base.issymmetric
end

using Base.LinAlg

import Base.show

export PardisoSolver, MKLPardisoSolver
export set_iparm!, set_dparm!, set_matrixtype!, set_solver!, set_phase!, set_msglvl!, set_nprocs!
export get_iparm, get_iparms, get_dparm, get_dparms
export get_mtype, get_solver!, get_phase, get_msglvl, get_nprocs
export set_maxfct!, set_perm!, set_mnum!
export get_maxfct, get_perm, get_mnum
export checkmatrix, checkvec, printstats, pardisoinit, pardiso
export solve, solve!

include("CEnum.jl")

using .CEnum

type PardisoException <: Exception
    info::String
end

type PardisoPosDefException <: Exception
    info::String
end

Base.showerror(io::IO, e::Union{PardisoException, PardisoPosDefException}) = print(io, e.info);


typealias PardisoTypes Union{Float64, Complex128}

abstract AbstractPardisoSolver

function __init__()
    if !(MKL_PARDISO_LIB_FOUND || PARDISO_LIB_FOUND)
        warn("No Pardiso library managed to load.")
    end

    if PARDISO_LIB_FOUND
        try
            global const libpardiso = Libdl.dlopen(PARDISO_PATH, Libdl.RTLD_GLOBAL)
            global const libblas = Libdl.dlopen("libblas", Libdl.RTLD_GLOBAL)
            global const liblapack = Libdl.dlopen("liblapack", Libdl.RTLD_GLOBAL)
            global const libgfortran = Libdl.dlopen("libgfortran", Libdl.RTLD_GLOBAL)
            global const libgomp = Libdl.dlopen("libgomp", Libdl.RTLD_GLOBAL)
            global const init = Libdl.dlsym(libpardiso, "pardisoinit")
            global const pardiso_f = Libdl.dlsym(libpardiso, "pardiso")
            global const pardiso_chkmatrix = Libdl.dlsym(libpardiso, "pardiso_chkmatrix")
            global const pardiso_chkmatrix_z = Libdl.dlsym(libpardiso, "pardiso_chkmatrix_z")
            global const pardiso_printstats = Libdl.dlsym(libpardiso, "pardiso_printstats")
            global const pardiso_printstats_z = Libdl.dlsym(libpardiso, "pardiso_printstats_z")
            global const pardiso_chkvec = Libdl.dlsym(libpardiso, "pardiso_chkvec")
            global const pardiso_chkvec_z = Libdl.dlsym(libpardiso, "pardiso_chkvec_z")
            global const PARDISO_LOADED = true
        catch e
            println("Info: Pardiso did not load because: $e")
            global const PARDISO_LOADED = false
        end
    else
        global const PARDISO_LOADED = false
    end
end

include("../deps/deps.jl")
include("enums.jl")
include("project_pardiso.jl")
include("mkl_pardiso.jl")

# Getters and setters
set_matrixtype!(ps::AbstractPardisoSolver, v::Int) = set_matrixtype!(ps, MatrixType[v][1])
function set_matrixtype!(ps::AbstractPardisoSolver, v::MatrixType)
    ps.mtype = v
end

get_mtype(ps::AbstractPardisoSolver) = ps.mtype
get_iparm(ps::AbstractPardisoSolver, i::Integer) = ps.iparm[i]
get_iparms(ps::AbstractPardisoSolver) = ps.iparm
set_iparm!(ps::AbstractPardisoSolver, i::Integer, v::Integer) = ps.iparm[i] = v

get_mnum(ps::AbstractPardisoSolver) = ps.mnum
set_mnum!(ps::AbstractPardisoSolver, mnum::Integer) = ps.mnum = mnum

get_maxfct(ps::AbstractPardisoSolver) = ps.maxfct
set_maxfct!(ps::AbstractPardisoSolver, maxfct::Integer) = ps.maxfct = maxfct

get_perm(ps::AbstractPardisoSolver) = ps.perm
set_perm!{T <: Integer}(ps::PardisoTypes, perm::Vector{T}) = ps.perm = convert(Vector{Int32}, perm)

get_phase(ps::AbstractPardisoSolver) = ps.phase

set_phase!(ps::AbstractPardisoSolver, v::Int) = set_phase!(ps, Phase[v][1])
function set_phase!(ps::AbstractPardisoSolver, v::Phase)
    ps.phase = v
end

get_msglvl(ps::AbstractPardisoSolver) = ps.msglvl

set_msglvl!(ps::AbstractPardisoSolver, v::Integer) = set_msglvl!(ps, MessageLevel[v][1])
function set_msglvl!(ps::AbstractPardisoSolver, v::MessageLevel)
    ps.msglvl = v
end

function pardisoinit(ps::AbstractPardisoSolver)
    ccall_pardisoinit(ps)
    return
end


function solve{Ti, Tv <: PardisoTypes}(ps::AbstractPardisoSolver, A::SparseMatrixCSC{Tv, Ti},
                                       B::VecOrMat{Tv}, T::Symbol=:N)
  X = copy(B)
  solve!(ps, X, A, B, T)
  return X
end

function solve!{Ti, Tv <: PardisoTypes}(ps::AbstractPardisoSolver, X::VecOrMat{Tv},
                                        A::SparseMatrixCSC{Tv, Ti}, B::VecOrMat{Tv},
                                        T::Symbol=:N)

    pardisoinit(ps)

    # We need to set the transpose flag in PARDISO when we DON'T want
    # a transpose in Julia because we are passing a CSC formatted
    # matrix to PARDISO which expects a CSR matrix.
    if T == :N
        if typeof(ps) == PardisoSolver
            set_iparm!(ps, 12, 1)
        else
            set_iparm!(ps, 12, 2)
        end
    elseif T == :C || T == :T
        set_iparm!(ps, 12, 0)
    else
        throw(ArgumentError("only :T, :N  and :C, are valid transpose symbols"))
    end

    # If we want the matrix to only be transposed and not conjugated
    # we have to conjugate it before sening it to Pardiso due to CSC CSR
    # mismatch.

    # This is the heuristics for choosing what matrix type to use
    ##################################################################
    # - If hermitian try to solve with symmetroc positive definite.
    #   - On pos def exception, solve instead with symmetric indefinite.
    # - If complex and symmetric, solve with symmetric complex solver
    # - Else solve as unsymmetric.
     if ishermitian(A)
        eltype(A) == Float64 ? set_matrixtype!(ps, REAL_SYM_POSDEF) : set_matrixtype!(ps, COMPLEX_HERM_POSDEF)
        try
            if typeof(ps) == PardisoSolver
                   pardiso(ps, X, get_matrix(ps, A, T), B)
            else
                # Workaround for #3
                if eltype(A) == Complex128
                    throw(PardisoPosDefException(""))
                end
                pardiso(ps, X, get_matrix(ps, A, T), B)
            end
        catch e
            isa(e, PardisoPosDefException) || rethrow(e)
            eltype(A) == Float64 ? set_matrixtype!(ps, REAL_SYM_INDEF) : set_matrixtype!(ps, COMPLEX_HERM_INDEF )
            pardiso(ps, X, get_matrix(ps, A, T), B)
        end
    elseif issymmetric(A)
        set_matrixtype!(ps, COMPLEX_SYM)
        pardiso(ps, X, get_matrix(ps, A, T), B)
    else
        eltype(A) == Float64 ? set_matrixtype!(ps, REAL_NONSYM) : set_matrixtype!(ps, COMPLEX_NONSYM)
        pardiso(ps, X, get_matrix(ps, A, T), B)
    end
    original_phase = get_phase(ps)

    # Release memory, TODO: We are running the convert on IA and JA here
    # again which is unnecessary.
    set_phase!(ps, RELEASE_ALL)
    pardiso(ps, X, A, B)
    set_phase!(ps, original_phase)
    return X
end

function pardiso{Ti, Tv <: PardisoTypes}(ps::AbstractPardisoSolver, X::VecOrMat{Tv},
                                         A::SparseMatrixCSC{Tv, Ti}, B::VecOrMat{Tv})

    dim_check(X, A, B)

    if Tv <: Complex && isreal(get_mtype(ps))
        throw(ErrorException(string("input matrix is complex while PardisoSolver ",
                                    "has a real matrix type set")))
    end

    if Tv <: Real && !isreal(get_mtype(ps))
        throw(ErrorException(string("input matrix is real while PardisoSolver ",
                                    "has a complex matrix type set")))
    end

    N = Int32(size(A, 2))

    AA = A.nzval
    IA = convert(Vector{Int32}, A.colptr)
    JA = convert(Vector{Int32}, A.rowval)

    NRHS = Int32(size(B, 2))

    ccall_pardiso(ps, N, AA, IA, JA, NRHS, B, X)
end

function dim_check(X, A, B)
    size(X) == size(B) || throw(DimensionMismatch(string(
                                 "Solution has $(size(X)), ",
                                 "RHS has size as $(size(B)).")))
    size(A,1) == size(B,1) || throw(DimensionMismatch(string(
                                    "Matrix has $(size(A,1)) ",
                                    "rows, RHS has $(size(B,1)) rows.")))
end


end # module
