ENV["OMP_NUM_THREADS"] = 2

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

using Pardiso

srand(1234)

psolvers = DataType[]

Pardiso.MKL_PARDISO_LOADED && push!(psolvers, MKLPardisoSolver)
Pardiso.PARDISO_LOADED && push!(psolvers, PardisoSolver)

println("Testing ", psolvers)

if length(psolvers) == 0
    error("No Pardiso library managed to load. Unable to run tests.")
end

# Test solver + for real and complex data
@testset "solving" begin
for pardiso_type in psolvers
    for T in (Float64, Complex128)
        ps = pardiso_type()
        pardisoinit(ps)

        if T == Float64
            set_matrixtype!(ps, 11)
        else
            set_matrixtype!(ps, 13)
        end

        A1 = sparse(rand(T, 10,10))
        B = rand(T, 10, 2)
        X = similar(B)

        # Test unsymmetric, herm indef, herm posdef and symmetric
        for A in SparseMatrixCSC[A1, A1 + A1', A1'A1, A1 + A1.']

            solve!(ps, X, A, B)
            @test X ≈ A\B

            X = solve(ps, A, B)
            @test X ≈ A\B

            solve!(ps, X, A, B, :C)
            @test X ≈ A'\B

            X = solve(ps, A, B, :C)
            @test X ≈ A'\B

            solve!(ps, X, A, B, :T)
            @test X ≈ A.'\B

            X = solve(ps, A, B, :T)
            @test X ≈ A.'\B
        end
    end
end
end #testset

@testset "error checks" begin
for pardiso_type in psolvers

    ps = pardiso_type()

    A = sparse(rand(10,10))
    B = rand(10, 2)
    X = rand(10, 2)

    if typeof(pardiso_type) == PardisoSolver
        printstats(ps, A, B)
        checkmatrix(ps, A)
        checkvec(ps, B)
    end


    set_matrixtype!(ps, 13)
    @test_throws ErrorException pardiso(ps, X, A, B)
    @test_throws ArgumentError solve(ps, A, B, :P)
    @test_throws ArgumentError solve!(ps, X, A, B, :P)

    set_matrixtype!(ps, 11)
    X = zeros(12, 2)
    @test_throws DimensionMismatch solve!(ps,X, A, B)

    B = rand(12, 2)
    @test_throws DimensionMismatch solve(ps, A, B)
end
end # testset


@testset "getters and setters" begin
for pardiso_type in psolvers
    ps = pardiso_type()
    set_iparm!(ps, 1, 0)
    pardisoinit(ps)
    @test get_iparm(ps, 1) == 1

    @test_throws ArgumentError set_phase!(ps, 5)
    @test_throws ArgumentError set_msglvl!(ps, 2)
    @test_throws ArgumentError set_matrixtype!(ps, 15)

    if typeof(pardiso_type) == PardisoSolver
        @test_throws ArgumentError set_solver(ps, 2)

        set_dparm(ps, 5, 13.37)
        @test get_dparm(ps, 5) == 13.37

        set_solver(ps, 1)
        @test get_solver(ps) == 1
    end

    set_iparm!(ps, 13, 100)
    @test get_iparm(ps, 13) == 100

    set_matrixtype!(ps, Pardiso.REAL_SYM)
    @test get_mtype(ps) == Pardiso.REAL_SYM

    set_phase!(ps, Pardiso.ANALYSIS_NUM_FACT)
    @test get_phase(ps) == Pardiso.ANALYSIS_NUM_FACT

    set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
    @test get_msglvl(ps) == Pardiso.MESSAGE_LEVEL_ON
end
end # testset
