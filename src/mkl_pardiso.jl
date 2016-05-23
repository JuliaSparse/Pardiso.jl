const MKL_DOMAIN_PARDISO = @compat Int32(4)

type MKLPardisoSolver <: AbstractPardisoSolver
    pt::Vector{Int}
    iparm::Vector{Int32}
    mtype::MatrixType
    solver::Solver
    phase::Phase
    msglvl::MessageLevel
    maxfct::Int32
    mnum::Int32
    perm::Vector{Int32}
end

function MKLPardisoSolver()
    if !MKL_PARDISO_LOADED
        error("mkl library was not loaded, unable to create solver")
    end

    pt = zeros(Int, 64)
    iparm = zeros(Int32, 64)
    mtype = REAL_NONSYM
    solver = DIRECT_SOLVER
    phase = ANALYSIS_NUM_FACT_SOLVE_REFINE
    msglvl = MESSAGE_LEVEL_OFF
    mnum = 1
    maxfct = 1
    perm = Int32[]

    MKLPardisoSolver(pt, iparm, mtype, solver,
                      phase, msglvl, maxfct, mnum, perm)
end


show(io::IO, ps::MKLPardisoSolver) = print(io, string("$MKLPardisoSolver:\n",
                                  "\tMatrix type: $(MTYPES[get_mtype(ps)])\n",
                                  "\tPhase: $(PHASES[get_phase(ps)])\n"))

set_nprocs!(ps::MKLPardisoSolver, n::Integer) = ccall(set_nthreads, Void, (Ptr{Int32},Ptr{Int32}), &(@compat Int32(n)), &MKL_DOMAIN_PARDISO)
get_nprocs(ps::MKLPardisoSolver) = ccall(get_nthreads, Int32, (Ptr{Int32},), &MKL_DOMAIN_PARDISO)

valid_phases(ps::MKLPardisoSolver) = keys(MKL_PHASES)
phases(ps::MKLPardisoSolver) = MKL_PHASES

function get_matrix(ps::MKLPardisoSolver, A, T)
    mtype = get_mtype(ps)

    if isposornegdef(mtype)
        T == :C && return conj(tril(A))
        return tril(A)
    end

    if !symmetric(mtype)
        T == :C && return conj(A)
        return A
    end

    if mtype == COMPLEX_SYM
        T == :C && return conj(tril(A))
        return tril(A)
    end

    error("Unhandled matrix type")
end

function ccall_pardisoinit(ps::MKLPardisoSolver)
    ERR = Ref{Int32}(0)
    ccall(mkl_init, Void,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32}),
          ps.pt, &ps.mtype, ps.iparm)
    check_error(ps, ERR)
end


function ccall_pardiso{Tv}(ps::MKLPardisoSolver, N, AA::Vector{Tv}, IA, JA,
                                   NRHS, B::VecOrMat{Tv}, X::VecOrMat{Tv})
    ERR = Ref{Int32}(0)
    ccall(mkl_pardiso_f, Void,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Tv}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Tv},
           Ptr{Int32}),
          ps.pt, &ps.maxfct, &ps.mnum, &ps.mtype, &ps.phase,
          &N, AA, IA, JA, ps.perm,
          &NRHS, ps.iparm, &ps.msglvl, B, X,
          ERR)
    check_error(ps, ERR[])
end

function check_error(ps::MKLPardisoSolver, err::Integer)
    err != -1  || throw(PardisoException("Input inconsistent."))
    err != -2  || throw(PardisoException("Not enough memory."))
    err != -3  || throw(PardisoException("Reordering problem."))
    err != -4  || throw(PardisoPosDefException("Zero pivot, numerical fact. or iterative refinement problem."))
    err != -5  || throw(PardisoException("Unclassified (internal) error."))
    err != -6  || throw(PardisoException("Preordering failed (matrix types 11, 13 only)."))
    err != -7  || throw(PardisoException("Diagonal matrix is singular"))
    err != -8  || throw(PardisoException("32-bit integer overflow problem."))
    err != -9  || throw(PardisoException("Not enough memory for OOC."))
    err != -10 || throw(PardisoException("Error opening OOC files."))
    err != -11 || throw(PardisoException("Read/write error with OOC files."))
    err != -12 || throw(PardisoException("pardiso_64 called from 32-bit library"))
    return
end


