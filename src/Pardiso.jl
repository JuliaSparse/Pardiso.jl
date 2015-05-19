module Pardiso

using Compat

using Base.LinAlg
using Base.SparseMatrix

import Base.show

export PardisoSolver, MKLPardisoSolver
export set_iparm, set_dparm, set_mtype, set_solver, set_phase, set_msglvl
export get_iparm, get_iparms, get_dparm, get_dparms
export get_mtype, get_solver, get_phase, get_msglvl, get_nprocs
export set_maxfct, set_perm, set_mnum
export get_maxfct, get_perm, get_mnum
export checkmatrix, checkvec, printstats, pardisoinit, pardiso
export solve, solve!

const VALID_MTYPES = [1, 2, -2, 3, 4, -4, 6, 11, 13]
const REAL_MTYPES = [1, 2, -2, 11]
const COMPLEX_MTYPES = [3, 4, -4, 6, 13]
const VALID_MSGLVLS = [0, 1]

type PardisoException <: Exception
    info::ASCIIString
end

type PardisoPosDefException <: Exception
    info::ASCIIString
end

Base.showerror(io::IO, e::Union(PardisoException, PardisoPosDefException)) = print(io, e.info);


@compat const MTYPES = Dict{Int, ASCIIString}(
  1 => "Real structurally symmetric",
  2 => "Real symmetric positive definite",
 -2 => "Real symmetric indefinite",
  3 => "Complex structurally symmetric",
  4 => "Complex Hermitian postive definite",
 -4 => "Complex Hermitian indefinite",
  6 => "Complex symmetric",
 11 => "Real nonsymmetric",
 13 => "Complex nonsymmetric")


typealias PardisoTypes Union(Float64, Complex128)

abstract AbstractPardisoSolver

include("pardiso.jl")
include("mkl_pardiso.jl")

# Getters and setters
function set_mtype(ps::AbstractPardisoSolver, v::Integer)
    v in VALID_MTYPES || throw(ArgumentError(string(
                                    "invalid matrix type, valid matrix ",
                                    "types are $VALID_MTYPES")))
    ps.mtype = v
end
get_mtype(ps::AbstractPardisoSolver) = ps.mtype

get_iparm(ps::AbstractPardisoSolver, i::Integer) = ps.iparm[i]
get_iparms(ps::AbstractPardisoSolver) = ps.iparm
set_iparm(ps::AbstractPardisoSolver, i::Integer, v::Integer) = ps.iparm[i] = v

get_mnum(ps::AbstractPardisoSolver) = ps.mnum
set_mnum(ps::AbstractPardisoSolver, mnum::Integer) = ps.mnum = mnum

get_maxfct(ps::AbstractPardisoSolver) = ps.maxfct
set_maxfct(ps::AbstractPardisoSolver, maxfct::Integer) = ps.maxfct = maxfct

get_perm(ps::AbstractPardisoSolver) = ps.perm
set_perm{T <: Integer}(ps::PardisoTypes, perm::Vector{T}) = ps.perm = convert(Vector{Int32}, perm)

get_phase(ps::AbstractPardisoSolver) = ps.phase

function set_phase(ps::AbstractPardisoSolver, v::Integer)
    v in valid_phases(ps) || throw(ArgumentError(string(
                                    "invalid phase, valid phases ",
                                    "are \n $(valid_phases(ps))")))
    ps.phase = v
end

get_msglvl(ps::AbstractPardisoSolver) = ps.msglvl
function set_msglvl(ps::AbstractPardisoSolver, v::Integer)
    v in VALID_MSGLVLS || throw(ArgumentError(string(
                                "invalid message level, valid message levels ",
                                "are $VALID_MSGLVLS")))
    ps.msglvl = v
end


function pardisoinit(ps::AbstractPardisoSolver)
    ccall_pardisoinit(ps)
    return
end


function solve{Ti, Tv <: PardisoTypes}(ps::PardisoSolver, A::SparseMatrixCSC{Tv, Ti},
                                       B::VecOrMat{Tv}, T::Symbol=:N)
  X = copy(B)
  solve!(ps, X, A, B, T)
  return X
end

function solve!{Ti, Tv <: PardisoTypes}(ps::AbstractPardisoSolver, X::VecOrMat{Tv},
                                        A::SparseMatrixCSC{Tv, Ti}, B::VecOrMat{Tv},
                                        T::Symbol=:N)


    pardisoinit(ps)

    if (T != :N && T != :T)
        throw(ArgumentError("only :T and :N are valid transpose symbols"))
    end

    # We need to set the transpose flag in PARDISO when we DON'T want
    # a transpose in Julia because we are passing a CSC formatted
    # matrix to PARDISO which expects a CSR matrix.
    if T == :N
        set_transposed(ps, true)
    else
        set_transposed(ps, false)
    end

    original_phase = get_phase(ps)

    # If hermitian, try solve pos def, on error, solve with normal symm
    # else solve with unsymm
    if ishermitian(A)
        eltype(A) == Float64 ? set_mtype(ps, 2) : set_mtype(ps, 4)
        try
            pardiso(ps, X, A, B)
        catch e
            isa(e, PardisoPosDefException) || rethrow(e)
            eltype(A) == Float64 ? set_mtype(ps, -2) : set_mtype(ps, -4)
            pardiso(ps, X, A, B)
        end
    else
        eltype(A) == Float64 ? set_mtype(ps, 11) : set_mtype(ps, 13)
        pardiso(ps, X, A, B)
    end

    # Release memory, TODO: We are running the convert on IA and JA here
    # again which is unnecessary.
    set_phase(ps, -1)
    pardiso(ps, X, A, B)
    set_phase(ps, original_phase)
    return X
end

function pardiso{Ti, Tv <: PardisoTypes}(ps::AbstractPardisoSolver, X::VecOrMat{Tv},
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

