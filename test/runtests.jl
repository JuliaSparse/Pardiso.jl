using Pardiso
using Base.Test
using Base.SparseMatrix

srand(1234)
ENV["OMP_NUM_THREADS"] = 1

psolvers = DataType[]
Pardiso.MKL_PARDISO_LOADED && push!(psolvers, MKLPardisoSolver)
Pardiso.PARDISO_LOADED && push!(psolvers, PardisoSolver)

if length(psolvers) == 0
    error("No Pardiso library managed to load. Unable to run tests.")
end

# Test solver + checkers for real matrices
let
for pardiso_type in psolvers
    for data_type in [Float64, Complex128]
        ps = pardiso_type()
        pardisoinit(ps)

        if data_type == Float64
            set_mtype(ps, 11)
        else
            set_mtype(ps, 13)
        end

        A1 = sparse(rand(data_type, 10,10))
        B = rand(data_type, 10, 2)
        X = similar(B)

        # Test unsymmetric, symmetric indef and symmetric posdef
        for A in SparseMatrixCSC[A1, A1 + A1', A1'A1, A1 + A1.']

            solve!(ps, X, A, B)
            @test_approx_eq X A\B
            fill!(X, 0.0)

            X = solve(ps, A, B)
            @test_approx_eq X A\B
            fill!(X, 0.0)

            solve!(ps, X, A, B, :C)
            @test_approx_eq X A'\B
            fill!(X, 0.0)

            X = solve(ps, A, B, :C)
            @test_approx_eq X A'\B
            fill!(X, 0.0)

            solve!(ps, X, A, B, :T)
            @test_approx_eq X A.'\B
            fill!(X, 0.0)

            X = solve(ps, A, B, :T)
            @test_approx_eq X A.'\B
            fill!(X, 0.0)
        end
    end
end
end

# Test some errors
let
for pardiso_type in psolvers

    ps = pardiso_type()

    A = sparse(rand(10,10))
    B = rand(10, 2)
    X = rand(10, 2)

    if pardiso_type == PardisoSolver
        printstats(ps, A, B)
        checkmatrix(ps, A)
        checkvec(ps, B)
    end


    set_mtype(ps, 13)
    @test_throws ErrorException pardiso(ps, X, A, B)
    @test_throws ArgumentError solve(ps, A, B, :P)
    @test_throws ArgumentError solve!(ps, X, A, B, :P)

    set_mtype(ps, 11)
    X = zeros(12, 2)
    @test_throws DimensionMismatch solve!(ps,X, A, B)

    B = rand(12, 2)
    @test_throws DimensionMismatch solve(ps, A, B)
end
end


let
for pardiso_type in psolvers
    ps = PardisoSolver()
    set_iparm(ps, 1, 0)
    pardisoinit(ps)
    @test get_iparm(ps, 1) == 1

    @test_throws ArgumentError set_solver(ps, 2)
    @test_throws ArgumentError set_phase(ps, 5)
    @test_throws ArgumentError set_msglvl(ps, 2)
    @test_throws ArgumentError set_mtype(ps, 15)

    set_solver(ps, 1)
    @test get_solver(ps) == 1

    set_dparm(ps, 5, 13.37)
    @test get_dparm(ps, 5) == 13.37

    set_iparm(ps, 13, 100)
    @test get_iparm(ps, 13) == 100

    set_mtype(ps, 1)
    @test get_mtype(ps) == 1

    set_phase(ps, 12)
    @test get_phase(ps) == 12

    set_msglvl(ps, 1)
    @test get_msglvl(ps) == 1
end
end
