__precompile__()

module Pardiso

if !isfile(joinpath(@__DIR__, "..", "deps", "deps.jl"))
    error("""please run Pkg.build("Pardiso") before loading the package""")
end

include("../deps/deps.jl")

function show_build_log()
    logfile = joinpath(@__DIR__, "..", "deps", "build.log")
    if !isfile(logfile)
        error("no log file found")
    else
        println(read(logfile, String))
    end
end

using Libdl
using SparseArrays
using LinearAlgebra
import Base.show

if !LOCAL_MKL_FOUND
    import MKL_jll
end

MKL_LOAD_FAILED = false

mkl_is_available() = (LOCAL_MKL_FOUND || MKL_jll.is_available()) && !MKL_LOAD_FAILED

if LinearAlgebra.BLAS.vendor() === :mkl && LinearAlgebra.BlasInt == Int64
    const MklInt = Int64
    const PARDISO_FUNC = :pardiso_64
else
    const MklInt = Int32
    const PARDISO_FUNC = :pardiso
end

const libmkl_rt = Ref{String}("")

export PardisoSolver, MKLPardisoSolver
export set_iparm!, set_dparm!, set_matrixtype!, set_solver!, set_phase!, set_msglvl!, set_nprocs!
export get_iparm, get_iparms, get_dparm, get_dparms
export get_matrixtype, get_solver, get_phase, get_msglvl, get_nprocs
export set_maxfct!, set_perm!, set_mnum!
export get_maxfct, get_perm, get_mnum
export checkmatrix, checkvec, printstats, pardisoinit, pardiso
export solve, solve!
export get_matrix
export schur_complement, pardisogetschur
export fix_iparm!
export  mkl_is_available, panua_is_available

struct PardisoException <: Exception
    info::String
end

struct PardisoPosDefException <: Exception
    info::String
end

Base.showerror(io::IO, e::Union{PardisoException,PardisoPosDefException}) = print(io, e.info);


const PardisoNumTypes = Union{Float32, ComplexF32, Float64, ComplexF64}

abstract type AbstractPardisoSolver end

function collect_gfortran_lib_candidates(main_ver)
    candidates = String[]
    if Sys.isapple()
        homebrew_gcc = "/usr/local/Cellar/gcc/"
        isdir(homebrew_gcc) || return String[]
        vers = readdir(homebrew_gcc)
        filter!(x -> startswith(x, "$main_ver."), vers)
        for v in vers
            push!(candidates, joinpath(homebrew_gcc, v, "lib/gcc/$main_ver/"))
        end
    elseif Sys.islinux()
        gcc_path = "/usr/lib/gcc/x86_64-linux-gnu/"
        isdir(gcc_path) || return String[]
        vers = readdir(gcc_path)
        filter!(x -> startswith(x, "$main_ver.") || isequal(x, "$main_ver"), vers)
        for v in vers
            push!(candidates, joinpath(gcc_path, v))
        end
    end
    return candidates
end

load_lib_fortran(lib::String, v::Int) = load_lib_fortran(lib, [v])
function load_lib_fortran(lib::String, vs::Vector{Int})
    candidates = String[]
    for v in vs
        append!(candidates, collect_gfortran_lib_candidates(v))
    end
    path = Libdl.find_library(lib, candidates)
    isempty(path) && (path = lib)
    Libdl.dlopen(path * "." * Libdl.dlext, Libdl.RTLD_GLOBAL)
end

# Pardiso
const init = Ref{Ptr}()
const pardiso_f = Ref{Ptr}()
const pardiso_chkmatrix = Ref{Ptr}()
const pardiso_chkmatrix_z = Ref{Ptr}()
const pardiso_printstats = Ref{Ptr}()
const pardiso_printstats_z = Ref{Ptr}()
const pardiso_chkvec = Ref{Ptr}()
const pardiso_chkvec_z = Ref{Ptr}()
const pardiso_get_schur_f = Ref{Ptr}()
const PARDISO_LOADED = Ref(false)

panua_is_available() = PARDISO_LOADED[]

function __init__()
    global MKL_LOAD_FAILED
    if LOCAL_MKL_FOUND
        if Sys.iswindows()
            libmkl_rt[] = "mkl_rt"
        elseif Sys.isapple()
            libmkl_rt[] = "@rpath/libmkl_rt.dylib"
        else
            libmkl_rt[] = "libmkl_rt"
        end
    elseif MKL_jll.is_available()
        libmkl_rt[] = MKL_jll.libmkl_rt_path
    end

    if !haskey(ENV, "PARDISOLICMESSAGE")
        ENV["PARDISOLICMESSAGE"] = 1
    end

    if LOCAL_MKL_FOUND && !haskey(ENV, "MKLROOT")
        @warn "MKLROOT not set, MKL Pardiso solver will not be functional"
    end

    if mkl_is_available()
        try
            libmklpardiso = Libdl.dlopen(libmkl_rt[])
            mklpardiso_f = Libdl.dlsym(libmklpardiso, "pardiso")
        catch e
            @error("MKL Pardiso did not manage to load, error thrown was: $(sprint(showerror, e))")
            MKL_LOAD_FAILED = true
        end
    end

    # This is apparently needed for MKL to not get stuck on 1 thread when
    # libpardiso is loaded in the block below...
    if libmkl_rt[] !== ""
        get_nprocs_mkl()
    end

    if PARDISO_LIB_FOUND
        try
            libpardiso = Libdl.dlopen(PARDISO_PATH)
            init[] = Libdl.dlsym(libpardiso, "pardisoinit")
            pardiso_f[] = Libdl.dlsym(libpardiso, "pardiso")
            pardiso_chkmatrix[] = Libdl.dlsym(libpardiso, "pardiso_chkmatrix")
            pardiso_chkmatrix_z[] = Libdl.dlsym(libpardiso, "pardiso_chkmatrix_z")
            pardiso_printstats[] = Libdl.dlsym(libpardiso, "pardiso_printstats")
            pardiso_printstats_z[] = Libdl.dlsym(libpardiso, "pardiso_printstats_z")
            pardiso_chkvec[] = Libdl.dlsym(libpardiso, "pardiso_chkvec")
            pardiso_chkvec_z[] = Libdl.dlsym(libpardiso, "pardiso_chkvec_z")
            pardiso_get_schur_f[] = Libdl.dlsym(libpardiso, "pardiso_get_schur")

            PARDISO_LOADED[] = true
        catch e
            @error("Pardiso did not manage to load, error thrown was: $(sprint(showerror, e))")
        end
    end
end
include("enums.jl")
include("panua_pardiso.jl")
include("mkl_pardiso.jl")

# Getters and setters
set_matrixtype!(ps::AbstractPardisoSolver, v::Int) = set_matrixtype!(ps, MatrixType(v))
function set_matrixtype!(ps::AbstractPardisoSolver, v::MatrixType)
    ps.mtype = v
end

get_matrixtype(ps::AbstractPardisoSolver) = ps.mtype
get_iparm(ps::AbstractPardisoSolver, i::Integer) = ps.iparm[i]
get_iparms(ps::AbstractPardisoSolver) = ps.iparm
set_iparm!(ps::AbstractPardisoSolver, i::Integer, v::Integer) = ps.iparm[i] = v

get_mnum(ps::AbstractPardisoSolver) = ps.mnum
set_mnum!(ps::AbstractPardisoSolver, mnum::Integer) = ps.mnum = mnum

get_maxfct(ps::AbstractPardisoSolver) = ps.maxfct
set_maxfct!(ps::AbstractPardisoSolver, maxfct::Integer) = ps.maxfct = maxfct

get_perm(ps::AbstractPardisoSolver) = ps.perm
set_perm!(ps::AbstractPardisoSolver, perm::Vector{T}) where {T <: Integer} = ps.perm = convert(Vector{Int32}, perm)

get_phase(ps::AbstractPardisoSolver) = ps.phase

set_phase!(ps::AbstractPardisoSolver, v::Int) = set_phase!(ps, Phase(v))
function set_phase!(ps::AbstractPardisoSolver, v::Phase)
    ps.phase = v
end

get_msglvl(ps::AbstractPardisoSolver) = ps.msglvl

set_msglvl!(ps::AbstractPardisoSolver, v::Integer) = set_msglvl!(ps, MessageLevel(v))
function set_msglvl!(ps::AbstractPardisoSolver, v::MessageLevel)
    ps.msglvl = v
end

function pardisoinit(ps::AbstractPardisoSolver)
    ccall_pardisoinit(ps)
    return
end


function solve(ps::AbstractPardisoSolver, A::SparseMatrixCSC{Tv,Ti},
               B::StridedVecOrMat{Tv}, T::Symbol=:N) where {Ti, Tv <: PardisoNumTypes}
    X = copy(B)
    solve!(ps, X, A, B, T)
    return X
end

function fix_iparm!(ps::AbstractPardisoSolver, T::Symbol)
    # We need to set the transpose flag in PARDISO when we DON'T want
    # a transpose in Julia because we are passing a CSC formatted
    # matrix to PARDISO which expects a CSR matrix.
    if T === :N
        if isa(ps, PardisoSolver)
            set_iparm!(ps, 12, 1)
        else
            # iparm[12] = 1 is complex conjugate in MKL
            set_iparm!(ps, 12, 2)
        end
    elseif T === :C || T === :T
        set_iparm!(ps, 12, 0)
    else
        throw(ArgumentError("only :T, :N and :C, are valid transpose symbols"))
    end
end

# Copied from SparseArrays.jl but with a tweak
# to check symmetry of the sparsity pattern
function _is_hermsym(A::SparseMatrixCSC, check::Function)
    m, n = size(A)
    if m != n; return false; end

    colptr = SparseArrays.getcolptr(A)
    rowval = rowvals(A)
    nzval = nonzeros(A)
    tracker = copy(SparseArrays.getcolptr(A))
    for col in axes(A,2)
        # `tracker` is updated such that, for symmetric matrices,
        # the loop below starts from an element at or below the
        # diagonal element of column `col`"
        for p = tracker[col]:colptr[col+1]-1
            val = nzval[p]
            row = rowval[p]

            # Ignore stored zeros
            if iszero(val)
                continue
            end

            # If the matrix was symmetric we should have updated
            # the tracker to start at the diagonal or below. Here
            # we are above the diagonal so the matrix can't be symmetric.
            if row < col
                return false
            end

            # Diagonal element
            if row == col
                if !check(val, val)
                    return false
                end
            else
                offset = tracker[row]

                # If the matrix is unsymmetric, there might not exist
                # a rowval[offset]
                if offset > length(rowval)
                    return false
                end

                row2 = rowval[offset]

                # row2 can be less than col if the tracker didn't
                # get updated due to stored zeros in previous elements.
                # We therefore "catch up" here while making sure that
                # the elements are actually zero.
                while row2 < col
                    if _isnotzero(nzval[offset])
                        return false
                    end
                    offset += 1
                    row2 = rowval[offset]
                    tracker[row] += 1
                end

                # Non zero A[i,j] exists but A[j,i] does not exist
                if row2 > col
                    return false
                end

                # A[i,j] and A[j,i] exists
                if row2 == col
                    if !check(val, nzval[offset])
                        return false
                    end
                    tracker[row] += 1
                end
            end
        end
    end
    return true
end

isstructurallysymmetric(A::SparseMatrixCSC) = _is_hermsym(A, (x,y) -> true)

function solve!(ps::AbstractPardisoSolver, X::StridedVecOrMat{Tv},
                A::SparseMatrixCSC{Tv,Ti}, B::StridedVecOrMat{Tv},
                T::Symbol=:N) where {Ti, Tv <: PardisoNumTypes}
    LinearAlgebra.checksquare(A)
    set_phase!(ps, ANALYSIS_NUM_FACT_SOLVE_REFINE)

    # This is the heuristics for choosing what matrix type to use
    ##################################################################
    # - If hermitian try to solve with symmetric positive definite.
    #   - On pos def exception, solve instead with symmetric indefinite.
    # - If complex and symmetric, solve with symmetric complex solver
    # - Else solve as unsymmetric.
    if Tv <: Union{Float32, Float64}
        if issymmetric(A)
            set_matrixtype!(ps, REAL_SYM_POSDEF)
            pardisoinit(ps)
            fix_iparm!(ps, T)
            Tv == Float32 && set_iparm!(ps, 28, 1)
            try
                pardiso(ps, X, get_matrix(ps, A, T), B)
            catch e
                set_phase!(ps, RELEASE_ALL)
                pardiso(ps, X, A, B)
                set_phase!(ps, ANALYSIS_NUM_FACT_SOLVE_REFINE)
                if !isa(e, PardisoPosDefException)
                    rethrow()
                end
                set_matrixtype!(ps, REAL_SYM_INDEF)
                pardisoinit(ps)
                fix_iparm!(ps, T)
                Tv == Float32 && set_iparm!(ps, 28, 1)
                pardiso(ps, X, get_matrix(ps, A, T), B)
            end
        else
            if isstructurallysymmetric(A)
                set_matrixtype!(ps, REAL_SYM)
            else
                set_matrixtype!(ps, REAL_NONSYM)
            end
            pardisoinit(ps)
            fix_iparm!(ps, T)
            Tv == Float32 && set_iparm!(ps, 28, 1)
            pardiso(ps, X, get_matrix(ps, A, T), B)
        end
    else # Tv <: Union{ComplexF64, ComplexF32}
        if ishermitian(A)
            set_matrixtype!(ps, COMPLEX_HERM_POSDEF)
            pardisoinit(ps)
            fix_iparm!(ps, T)
            Tv == ComplexF32 && set_iparm!(ps, 28, 1)
            try
                pardiso(ps, X, get_matrix(ps, A, T), B)
            catch e
                set_phase!(ps, RELEASE_ALL)
                pardiso(ps, X, A, B)
                set_phase!(ps, ANALYSIS_NUM_FACT_SOLVE_REFINE)
                if !isa(e, PardisoPosDefException)
                    rethrow()
                end
                set_matrixtype!(ps, COMPLEX_HERM_INDEF)
                pardisoinit(ps)
                fix_iparm!(ps, T)
                Tv == ComplexF32 && set_iparm!(ps, 28, 1)
                pardiso(ps, X, get_matrix(ps, A, T), B)
            end
        else
            if issymmetric(A)
                set_matrixtype!(ps, COMPLEX_SYM)
            elseif isstructurallysymmetric(A)
                set_matrixtype!(ps, COMPLEX_STRUCT_SYM)
            else
                set_matrixtype!(ps, COMPLEX_NONSYM)
            end
            pardisoinit(ps)
            fix_iparm!(ps, T)
            Tv == ComplexF32 && set_iparm!(ps, 28, 1)
            pardiso(ps, X, get_matrix(ps, A, T), B)
        end
    end

    # Release memory, TODO: We are running the convert on IA and JA here
    # again which is unnecessary.
    set_phase!(ps, RELEASE_ALL)
    pardiso(ps, X, A, B)

    set_phase!(ps, ANALYSIS_NUM_FACT_SOLVE_REFINE)
    return X
end

function get_matrix(ps::AbstractPardisoSolver, A, T)
    mtype = get_matrixtype(ps)

    if isposornegdef(mtype)
        if ps isa MKLPardisoSolver
            T == :C && return conj(tril(A))
            return tril(A)
        elseif ps isa PardisoSolver
            T == :T && return tril(A)
            return conj(tril(A))
        end
    end

    if !issymmetric(mtype)
        T == :C && return conj(A)
        return A
    end

    if mtype == COMPLEX_SYM
        T == :C && return conj(tril(A))
        return tril(A)
    end

    error("Unhandled matrix type")
end

function pardiso(ps::AbstractPardisoSolver, X::StridedVecOrMat{Tv}, A::SparseMatrixCSC{Tv,Ti},
                 B::StridedVecOrMat{Tv}) where {Ti, Tv <: PardisoNumTypes}
    if length(X) != 0
        dim_check(X, A, B)
    end

    if Tv <: Complex && isreal(get_matrixtype(ps))
        throw(ErrorException(string("input matrix is complex while PardisoSolver ",
                                    "has a real matrix type set: $(get_matrixtype(ps))")))
    end

    if Tv <: Real && !isreal(get_matrixtype(ps))
        throw(ErrorException(string("input matrix is real while PardisoSolver ",
                                    "has a complex matrix type set: $(get_matrixtype(ps))")))
    end

    if Tv <: Union{Float32, ComplexF32} && typeof(ps) <: MKLPardisoSolver && ps.iparm[28] != 1
        throw(ErrorException(string("input matrix is Float32/ComplexF32 while MKLPardisoSolver ",
                                    "have iparm[28]=$(ps.iparm[28]) rather than 1.")))
    end

    N = size(A, 2)

    resize!(ps.perm, size(B, 1))

    NRHS = size(B, 2)

    ccall_pardiso(ps, N, A.nzval, A.colptr, A.rowval, NRHS, B, X)
end

pardiso(ps::AbstractPardisoSolver) = ccall_pardiso(ps, Int32(0), Float64[], Int32[], Int32[], Int32(0), Float64[], Float64[])
function pardiso(ps::AbstractPardisoSolver, A::SparseMatrixCSC{Tv,Ti}, B::StridedVecOrMat{Tv}) where {Ti, Tv <: PardisoNumTypes}
    pardiso(ps, Tv[], A, B)
end

# populated rows of S determine schur complment block
"""
    schur_complement(ps,A,x) -> S

Schur complement `S` of the submatrix defined by the nonzero entries of `x` in matrix `A`.
If `n=nnz(x)`, then `S` is `n`-by-`n`.

WARNING: for complex `M`, seems to be unstable, made worse as number of nonzero elements in `M` decreases
"""
schur_complement(ps::AbstractPardisoSolver,A,x::SparseVector,T::Symbol=:N) = _schur_complement_permuted(ps,A,x.nzind,T)
schur_complement(ps::AbstractPardisoSolver,A,x::SparseMatrixCSC,T::Symbol=:N) = _schur_complement_permuted(ps,A,unique!(sort!(x.rowval)),T)

# permute A and then compute complement of lower right-hand `n`-by-`n` block
function _schur_complement_permuted(ps,A,rows,T::Symbol)
    P = sparse(vcat(setdiff(1:A.n,rows),rows),1:A.n,1,size(A)...)
    schur_complement(ps,P'*A*P,length(rows),T)
end
# or integer gives last n rows and columns as schur complement block
"""
    schur_complement(ps,M,n) -> S

Schur complement `S` of upper-left block in `M`, where `n` is the size of lower-right block (and therefore also of Schur complement)

WARNING: for complex `M`, seems to be unstable, made worse as number of nonzero elements in `M` decreases
"""
function schur_complement(ps::AbstractPardisoSolver,A::SparseMatrixCSC{Tv},n::Integer,T::Symbol=:N) where Tv <: PardisoNumTypes

    n ≥ size(A,1) ? throw(ErrorException("complement block size n=$n≥A.m=$(A.m)")) : nothing
    # Tv<:Complex ? (@warn "unstable for complex types, unknown why") : nothing

    pardisoinit(ps)
    original_phase = get_phase(ps)
    original_iparms = get_iparms(ps)
    set_iparm!(ps,1,1) # use custom IPARM
    set_iparm!(ps,38,n) # set Schur complement block size to n
    set_phase!(ps,12) # analyze and factorize
    B = Matrix{Tv}(undef,size(A,1),0) # dummy array to feed to pardiso

    # transpose via IPARM(12) doesn't work at factorize step (only on entry to solve step)
    if T==:N
        M = permutedims(A)
        set_iparm!(ps, 12, 1)
    elseif T == :C
        M = conj(permutedims(A))
        set_iparm!(ps, 12, 0)
    elseif T == :T
        M = A
        set_iparm!(ps, 12, 0)
    else
        throw(ArgumentError("only :T, :N and :C, are valid transpose symbols"))
    end

    pardiso(ps,B,M,B)
    S = pardisogetschur(ps) # get schur complement matrix

    set_phase!(ps, RELEASE_ALL)
    pardiso(ps, B, M, B)
    set_phase!(ps, original_phase) # reset phase to user setting
    for i ∈ eachindex(original_iparms)
        set_iparm!(ps,i,original_iparms[i])
    end

    return S
end

"""
    pardisogetschur(ps) -> S

retrieve schur complement from PardisoSolver `ps`.
"""
function pardisogetschur(ps::AbstractPardisoSolver)
    nnzschur = get_iparm(ps, 39)
    nschur = get_iparm(ps,38)
    T = isreal(get_matrixtype(ps)) ? Float64 : ComplexF64
    if nnzschur==0
        return spzeros(T,nschur,nschur)
    else
        S = Vector{T}(undef,nnzschur)
        IS = Vector{Int32}(undef,nschur)
        JS = Vector{Int32}(undef,nnzschur)
        ccall_pardiso_get_schur(ps,S,IS,JS)
        IS = pushfirst!(IS,Int32(1)) # some issue with IS (nschur+1 doesn't seem to work)
        S = permutedims(SparseMatrixCSC(nschur,nschur,IS,JS,S)) # really constructing CSR and then transposing
        return S
    end
end

function dim_check(X, A, B)
    LinearAlgebra.checksquare(A)
    size(X) == size(B) || throw(DimensionMismatch(string("solution has $(size(X)), ",
                                                         "RHS has size as $(size(B)).")))
    size(A, 1) == size(B, 1) || throw(DimensionMismatch(string("matrix has $(size(A,1)) ",
                                                               "rows, RHS has $(size(B,1)) rows.")))
    size(B, 1) == stride(B, 2) || throw(DimensionMismatch(
                                            string("Only memory-contiguous RHS supported")))
end

end # module
