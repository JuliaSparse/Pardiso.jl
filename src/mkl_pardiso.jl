# MKL Pardiso functions
if MKL_PARDISO_LOADED
    global const mkl_init = Libdl.dlsym(libmkl_core, "pardisoinit")
    global const mkl_pardiso_f = Libdl.dlsym(libmkl_core, "pardiso")
end

const MKL_PHASES = Dict{Int, ASCIIString}(
 11 => "Analysis",
 12 => "Analysis, numerical factorization",
 13 => "Analysis, numerical factorization, solve, iterative refinement",
 22 => "Numerical factorization",
-22 => "Selected Inversion",
 23 => "Numerical factorization, solve, iterative refinement",
 33 => "Solve, iterative refinement",
331 => "like phase=33, but only forward substitution",
332 => "like phase=33, but only diagonal substitution (if available)",
333 => "like phase=33, but only backward substitution",
  0 => "Release internal memory for L and U matrix number MNUM",
 -1 => "Release all internal memory for all matrices")


type MKLPardisoSolver <: AbstractPardisoSolver
    pt::Vector{Int}
    iparm:: Vector{Int32}
    mtype::Int32
    solver::Int32
    phase::Int32
    msglvl::Int32
    maxfct::Int32
    mnum::Int32
    perm::Vector{Int32}
end

function MKLPardisoSolver()
    if !MKL_PARDISO_LOADED
        error("MKL Pardiso library could not be loaded")
    end

    pt = zeros(Int, 64)
    iparm = zeros(Int32, 64)
    mtype = 11 # Default to real unsymmetric matrices
    solver = 0 # Default to direct solver
    phase = 13 # Default to analysis + fact + solve + refine
    msglvl = 0
    mnum = 1
    maxfct = 1
    perm = Int32[]

    MKLPardisoSolver(pt, iparm, mtype, solver,
                      phase, msglvl, maxfct, mnum, perm)
end


show(io::IO, ps::MKLPardisoSolver) = print(io, string("$MKLPardisoSolver:\n",
                                  "\tMatrix type: $(MTYPES[get_mtype(ps)])\n",
                                  "\tPhase: $(PHASES[get_phase(ps)])\n"))


valid_phases(ps::MKLPardisoSolver) = keys(MKL_PHASES)
phases(ps::MKLPardisoSolver) = MKL_PHASES

@inline function ccall_pardisoinit(ps::MKLPardisoSolver)
    ccall(mkl_init, Void,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32}),
          ps.pt, &ps.mtype, ps.iparm)
    check_error(ERR)
end


@inline function ccall_pardiso(ps::MKLPardisoSolver, N, AA, IA, JA, NRHS, B, X)
    ERR = Int32[0]
    ccall(mkl_pardiso_f, Void,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Tv}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Tv},
           Ptr{Int32}, Ptr{Float64}),
          ps.pt, &ps.maxfct, &ps.mnum, &ps.mtype, &ps.phase,
          &N, AA, IA, JA, ps.perm,
          &NRHS, ps.iparm, &ps.msglvl, B, X,
          ERR)
    check_error(ps, ERR)
end

function check_error(ps::MKLPardisoSolver, err::Vector{Int32})
    err = err[1]
    err != -1  || throw(ErrorException("Input inconsistent."))
    err != -2  || throw(ErrorException("Not enough memory."))
    err != -3  || throw(ErrorException("Reordering problem."))
    err != -4  || throw(ErrorException("Zero pivot, numerical fact. or iterative refinement problem."))
    err != -5  || throw(ErrorException("Unclassified (internal) error."))
    err != -6  || throw(ErrorException("Preordering failed (matrix types 11, 13 only)."))
    err != -7  || throw(ErrorException("Diagonal matrix is singular"))
    err != -8  || throw(ErrorException("32-bit integer overflow problem."))
    err != -9  || throw(ErrorException("Not enough memory for OOC."))
    err != -10 || throw(ErrorException("Error opening OOC files."))
    err != -11 || throw(ErrorException("Read/write error with OOC files."))
    err != -12 || throw(ErrorException("pardiso_64 called from 32-bit library"))
    return
end

