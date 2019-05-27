mutable struct PardisoSolver <: AbstractPardisoSolver
    pt::Vector{Int}
    iparm::Vector{Int32}
    dparm::Vector{Float64}
    mtype::MatrixType
    solver::Solver
    phase::Phase
    msglvl::MessageLevel
    maxfct::Int32
    mnum::Int32
    perm::Vector{Int32}
end

function PardisoSolver()
    if !PARDISO_LOADED[]
      error("pardiso library was not loaded")
    end

    pt = zeros(Int, 64)
    iparm = zeros(Int32, 64)
    dparm = zeros(Float64, 64)
    mtype = REAL_NONSYM
    solver = DIRECT_SOLVER
    phase = ANALYSIS_NUM_FACT_SOLVE_REFINE
    msglvl = MESSAGE_LEVEL_OFF
    # Set numper of processors to CPU_CORES unless "OMP_NUM_THREADS" is set
    if haskey(ENV, "OMP_NUM_THREADS")
        iparm[3] = parse(Int, ENV["OMP_NUM_THREADS"])
    else
        iparm[3] = 1
    end

    mnum = 1
    maxfct = 1
    perm = Int32[]

    ps = PardisoSolver(pt, iparm, dparm, mtype, solver,
                  phase, msglvl, maxfct, mnum, perm)

    finalizer(ps) do ps
        set_phase!(ps, Pardiso.RELEASE_ALL)
        pardiso(ps)
    end

    return ps
end


show(io::IO, ps::PardisoSolver) = print(io, string("$PardisoSolver:\n",
                                  "\tSolver: $(SOLVER_STRING[get_solver(ps)])\n",
                                  "\tMatrix type: $(MATRIX_STRING[get_matrixtype(ps)])\n",
                                  "\tPhase: $(PHASE_STRING[get_phase(ps)])\n",
                                  "\tNum processors: $(get_nprocs(ps))"))



phases(ps::PardisoSolver) = PHASES

set_transposed(ps::PardisoSolver, t::Bool) = t ? set_iparm(ps, 12, 1) : set_iparm(ps, 12, 0)

get_dparm(ps::PardisoSolver, i::Integer) = ps.dparm[i]
get_dparms(ps::PardisoSolver) = ps.dparm
set_dparm!(ps::PardisoSolver, i::Integer, v::AbstractFloat) = ps.dparm[i] = v
get_nprocs(ps::PardisoSolver) = ps.iparm[3]

set_solver!(ps::PardisoSolver, v::Int) = set_solver!(ps, Solver(v))
function set_solver!(ps::PardisoSolver, v::Solver)
    ps.solver = v
end
get_solver(ps::PardisoSolver) = ps.solver


function get_matrix(ps::PardisoSolver, A, T)
    mtype = get_matrixtype(ps)

    if isposornegdef(mtype)
        T == :T && return tril(A)
        return conj(tril(A))
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

@inline function ccall_pardisoinit(ps::PardisoSolver)
    ERR = Ref{Int32}(0)
    ccall(init[], Cvoid,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Float64}, Ptr{Int32}),
          ps.pt, Ref(Int32(ps.mtype)), Ref(Int32(ps.solver)), ps.iparm, ps.dparm, ERR)
     check_error(ps, ERR[])
end


@inline function ccall_pardiso(ps::PardisoSolver, N::Int32, AA::Vector{Tv},
                                   IA, JA, NRHS::Int32, B::StridedVecOrMat{Tv}, X::StridedVecOrMat{Tv}) where {Tv}
    ERR = Ref{Int32}(0)
    ccall(pardiso_f[], Cvoid,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Tv}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Tv},
           Ptr{Int32}, Ptr{Float64}),
          ps.pt, Ref(ps.maxfct), Ref(Int32(ps.mnum)), Ref(Int32(ps.mtype)), Ref(Int32(ps.phase)),
          Ref(N), AA, IA, JA, ps.perm,
          Ref(NRHS), ps.iparm, Ref(Int32(ps.msglvl)), B, X,
          ERR, ps.dparm)
    check_error(ps, ERR[])
end



# Different checks
function printstats(ps::PardisoSolver, A::SparseMatrixCSC{Tv, Ti},
                    B::StridedVecOrMat{Tv}) where {Ti,Tv <: PardisoNumTypes}
    N = Int32(size(A, 2))
    AA = A.nzval
    IA = convert(Vector{Int32}, A.colptr)
    JA = convert(Vector{Int32}, A.rowval)
    NRHS = Int32(size(B, 2))
    ERR = Ref{Int32}(0)
    if Tv <: Complex
        f = pardiso_printstats_z[]
    else
        f = pardiso_printstats[]
    end
    ccall(f, Cvoid,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}, Ptr{Tv},
           Ptr{Int32}),
          Ref(Int32(ps.mtype)), Ref(N), AA, IA, JA, Ref(NRHS), B, ERR)

    check_error(ps, ERR[])
    return
end

function checkmatrix(ps::PardisoSolver, A::SparseMatrixCSC{Tv, Ti}) where {Ti,Tv <: PardisoNumTypes}
    N = Int32(size(A, 1))
    AA = A.nzval
    IA = convert(Vector{Int32}, A.colptr)
    JA = convert(Vector{Int32}, A.rowval)
    ERR = Ref{Int32}(0)

    if Tv <: Complex
        f = pardiso_chkmatrix_z[]
    else
        f = pardiso_chkmatrix[]
    end

    ccall(f, Cvoid,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}),
          Ref(Int32(ps.mtype)), Ref(N), AA, IA,
          JA, ERR)

    check_error(ps, ERR[])
    return
end

function checkvec(ps, B::StridedVecOrMat{Tv}) where {Tv <: PardisoNumTypes}
    N = Int32(size(B, 1))
    NRHS = Int32(size(B, 2))
    ERR = Int32[0]

    if Tv <: Complex
        f = pardiso_chkvec_z[]
    else
        f = pardiso_chkvec[]
    end
    ccall(f, Cvoid,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Int32}),
          Ref(N), Ref(NRHS), B, ERR)

    check_error(ps, ERR[])
    return
end

function check_error(ps::PardisoSolver, err::Integer)
    err != -1  || throw(PardisoException("Input inconsistent."))
    err != -2  || throw(PardisoException("Not enough memory."))
    err != -3  || throw(PardisoException("Reordering problem."))
    err != -4  || throw(PardisoPosDefException("Zero pivot, numerical fact. or iterative refinement problem."))
    err != -5  || throw(PardisoException("Unclassified (internal) error."))
    err != -6  || throw(PardisoException("Preordering failed (matrix types 11, 13 only)."))
    err != -7  || throw(PardisoException("Diagonal matrix problem."))
    err != -8  || throw(PardisoException("32-bit integer overflow problem."))
    err != -10 || throw(PardisoException("No license file pardiso.lic found."))
    err != -11 || throw(PardisoException("License is expired."))
    err != -12 || throw(PardisoException("Wrong username or hostname."))
    err != -100|| throw(PardisoException("Reached maximum number of Krylov-subspace iteration in iterative solver."))
    err != -101|| throw(PardisoException("No sufficient convergence in Krylov-subspace iteration within 25 iterations."))
    err != -102|| throw(PardisoException("Error in Krylov-subspace iteration."))
    err != -103|| throw(PardisoException("Break-Down in Krylov-subspace iteration."))
    return
end
