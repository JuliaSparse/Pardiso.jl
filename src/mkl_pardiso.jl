try
    const MKLROOT = ENV["MKLROOT"]
    if Int === Int64
        global const libmkl_gd = Libdl.dlopen(string(MKLROOT, "/lib/intel64/libmkl_gf_lp64"), Libdl.RTLD_GLOBAL)
        global const libmkl_seqential = Libdl.dlopen(string(MKLROOT, "/lib/intel64/libmkl_sequential"), Libdl.RTLD_GLOBAL)
        global const libmkl_core = Libdl.dlopen(string(MKLROOT, "/lib/intel64/libmkl_core"), Libdl.RTLD_GLOBAL)
    else
        # Untested!!
        global const libmkl_gd = Libdl.dlopen(string(MKLROOT, "/lib/ia32/libmkl_gf"), Libdl.RTLD_GLOBAL)
        global const libmkl_seqential = Libdl.dlopen(string(MKLROOT, "/lib/ia32/libmkl_sequential"), Libdl.RTLD_GLOBAL)
        global const libmkl_core = Libdl.dlopen(string(MKLROOT, "/lib/ia32/libmkl_core"), Libdl.RTLD_GLOBAL)
    end
    global const mkl_init = Libdl.dlsym(libmkl_gd, "pardisoinit")
    global const mkl_pardiso_f = Libdl.dlsym(libmkl_gd, "pardiso")
    global const set_nthreads = Libdl.dlsym(libmkl_gd, "mkl_domain_set_num_threads")
    global const get_nthreads = Libdl.dlsym(libmkl_gd, "mkl_domain_get_max_threads")
    global const MKL_PARDISO_LOADED = true
catch
    global const MKL_PARDISO_LOADED = false
end

const MKL_DOMAIN_PARDISO = Int32(4)


@compat const MKL_PHASES = Dict{Int, ASCIIString}(
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
        error("mkl library was not loaded")
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

set_nprocs(ps::MKLPardisoSolver, n::Integer) = ccall(set_nthreads, Void, (Ptr{Int32},Ptr{Int32}), &Int32(n), &MKL_DOMAIN_PARDISO)
get_nprocs(ps::MKLPardisoSolver) = ccall(get_nthreads, Int32, (Ptr{Int32},), &MKL_DOMAIN_PARDISO)

valid_phases(ps::MKLPardisoSolver) = keys(MKL_PHASES)
phases(ps::MKLPardisoSolver) = MKL_PHASES

set_transposed(ps::MKLPardisoSolver, t::Bool) = t ? set_iparm(ps, 12, 2) : set_iparm(ps, 12, 0)

@inline function ccall_pardisoinit(ps::MKLPardisoSolver)
    ERR = Int32[1]
    ccall(mkl_init, Void,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32}),
          ps.pt, &ps.mtype, ps.iparm)
    check_error(ps, ERR)
end


@inline function ccall_pardiso(ps::MKLPardisoSolver, N, AA, IA, JA, NRHS, B, X)
    ERR = Int32[0]
    ccall(mkl_pardiso_f, Void,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Tv}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Tv}, Ptr{Tv},
           Ptr{Int32}),
          ps.pt, &ps.maxfct, &ps.mnum, &ps.mtype, &ps.phase,
          &N, AA, IA, JA, ps.perm,
          &NRHS, ps.iparm, &ps.msglvl, B, X,
          ERR)
    check_error(ps, ERR)
end

function check_error(ps::MKLPardisoSolver, errv::Vector{Int32})
    err = errv[1]
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


