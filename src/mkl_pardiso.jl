const MKL_DOMAIN_PARDISO = Int32(4)

mutable struct MKLPardisoSolver <: AbstractPardisoSolver
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
    pt = zeros(Int, 64)
    iparm = zeros(Int32, 64)
    mtype = REAL_NONSYM
    solver = DIRECT_SOLVER
    phase = ANALYSIS_NUM_FACT_SOLVE_REFINE
    msglvl = MESSAGE_LEVEL_OFF
    mnum = 1
    maxfct = 1
    perm = Int32[]

    ps = MKLPardisoSolver(pt, iparm, mtype, solver,
                      phase, msglvl, maxfct, mnum, perm)

    return ps
end


show(io::IO, ps::MKLPardisoSolver) = print(io, string("$MKLPardisoSolver:\n",
                                  "\tMatrix type: $(MATRIX_STRING[get_matrixtype(ps)])\n",
                                  "\tPhase: $(PHASE_STRING[get_phase(ps)])"))

set_nprocs!(ps::MKLPardisoSolver, n::Integer) =
    ccall((:mkl_domain_set_num_threads, libmkl_rt), Cvoid, (Ptr{Int32}, Ptr{Int32}), Ref((Int32(n))), Ref(MKL_DOMAIN_PARDISO))
get_nprocs(ps::MKLPardisoSolver) =
    ccall((:mkl_domain_get_max_threads, libmkl_rt), Int32, (Ptr{Int32},), Ref(MKL_DOMAIN_PARDISO))

valid_phases(ps::MKLPardisoSolver) = keys(MKL_PHASES)
phases(ps::MKLPardisoSolver) = MKL_PHASES

function ccall_pardisoinit(ps::MKLPardisoSolver)
    ERR = Ref{Int32}(0)
    ccall((:pardisoinit, libmkl_rt), Cvoid,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32}),
          ps.pt, Ref(Int32(ps.mtype)), ps.iparm)
    check_error(ps, ERR[])
end


function ccall_pardiso(ps::MKLPardisoSolver, N, AA::Vector{Tv}, IA, JA,
                       NRHS, B::StridedVecOrMat{Tv}, X::StridedVecOrMat{Tv}) where {Tv}
    ERR = Ref{Int32}(0)
    ccall((:pardiso, libmkl_rt), Cvoid,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Tv}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Tv},
           Ptr{Int32}),
          ps.pt, Ref(ps.maxfct), Ref(ps.mnum), Ref(Int32(ps.mtype)), Ref(Int32(ps.phase)),
          Ref(N), AA, IA, JA, ps.perm,
          Ref(NRHS), ps.iparm, Ref(Int32(ps.msglvl)), B, X,
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

