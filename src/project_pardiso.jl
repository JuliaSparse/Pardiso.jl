# Pardiso functions
try
    global const libpardiso = Libdl.dlopen("libpardiso", Libdl.RTLD_GLOBAL)
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

const VALID_SOLVERS = [0, 1]

const SOLVERS = Dict{Int, ASCIIString}(
0 => "Direct",
1 => "Iterative")


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


type PardisoSolver <: AbstractPardisoSolver
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
    if !PARDISO_LOADED
      error("pardiso library was not loaded")
    end

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
        throw(ErrorException("OMP_NUM_THREADS not set"))
    end

    mnum = 1
    maxfct = 1
    perm = Int32[]

    PardisoSolver(pt, iparm, dparm, mtype, solver,
                  phase, msglvl, maxfct, mnum, perm)
end


show(io::IO, ps::PardisoSolver) = print(io, string("$PardisoSolver:\n",
                                  "\tSolver: $(SOLVERS[get_solver(ps)])\n",
                                  "\tMatrix type: $(MTYPES[get_mtype(ps)])\n",
                                  "\tPhase: $(PHASES[get_phase(ps)])\n",
                                  "\tNum processors: $(get_nprocs(ps))"))


valid_phases(ps::PardisoSolver) = keys(PHASES)
phases(ps::PardisoSolver) = PHASES

set_transposed(ps::PardisoSolver, t::Bool) = t ? set_iparm(ps, 12, 1) : set_iparm(ps, 12, 0)

get_dparm(ps::PardisoSolver, i::Integer) = ps.dparm[i]
get_dparms(ps::PardisoSolver) = ps.dparm
set_dparm(ps::PardisoSolver, i::Integer, v::AbstractFloat) = ps.dparm[i] = v
get_nprocs(ps::PardisoSolver) = ps.iparm[3]
function set_solver(ps::PardisoSolver, v::Integer)
    v in keys(SOLVERS) || throw(ArgumentError(string("invalid solver, valid solvers are 0 for",
                        " sparse direct solver, 1 for multi-recursive iterative solver")))
    ps.solver = v
end
get_solver(ps::PardisoSolver) = ps.solver


function get_matrix(ps::PardisoSolver, A, T)
    mtype = get_mtype(ps)

    if mtype in [2, 4, -2, -4]
        T == :T && return tril(A)
        return conj(tril(A))
    end

     if mtype in [11, 13]
        T == :C && return conj(A)
        return A
    end

    if mtype == 6
        T == :C && return conj(tril(A))
        return tril(A)
    end
end

@inline function ccall_pardisoinit(ps::PardisoSolver)
   ERR = Int32[0]
    ccall(init, Void,
          (Ptr{Int}, Ptr{Int32}, Ptr{Int32},
           Ptr{Int32}, Ptr{Float64}, Ptr{Int32}),
          ps.pt, &ps.mtype, &ps.solver, ps.iparm, ps.dparm, ERR)
    check_error(ps, ERR)
end


@inline function ccall_pardiso{Tv}(ps::PardisoSolver, N, AA::Vector{Tv},
                                   IA, JA, NRHS, B::VecOrMat{Tv}, X::VecOrMat{Tv})
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
    check_error(ps, ERR)
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

    check_error(ps, ERR)
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

    check_error(ps, ERR)
    return
end

function checkvec{Tv <: PardisoTypes}(ps, B::VecOrMat{Tv})
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

    check_error(ps, ERR)
    return
end

function check_error(ps::PardisoSolver, errv::Vector{Int32})
    err = errv[1]
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
