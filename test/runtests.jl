ENV["OMP_NUM_THREADS"] = 2


using Pkg
if Sys.isapple()
    Pkg.add(name="MKL_jll"; version = "2023")
end


using Test
using Pardiso
using Random
using SparseArrays
using LinearAlgebra

Random.seed!(1234)

available_solvers = empty([Pardiso.AbstractPardisoSolver])
if Pardiso.mkl_is_available()
    push!(available_solvers, MKLPardisoSolver)
else
    @warn "Not testing MKL Pardiso solver"
end
if Pardiso.PARDISO_LOADED[]
    push!(available_solvers, PardisoSolver)
else
    @warn "Not testing project Pardiso solver"
end

@show Pardiso.MklInt

println("Testing ", available_solvers)

supported_eltypes(ps::PardisoSolver) = (Float64, ComplexF64)
supported_eltypes(ps::MKLPardisoSolver) = (Float32, ComplexF32, Float64, ComplexF64)

# Test solver + for real and complex data
@testset "solving" begin
for pardiso_type in available_solvers
    ps = pardiso_type()
    for T in supported_eltypes(ps)
        A1 = sparse(rand(T, 10,10))
        for B in (rand(T, 10, 2), view(rand(T, 10, 4), 1:10, 2:3))
            X = similar(B)
            # Test unsymmetric, herm indef, herm posdef and symmetric
            for A in SparseMatrixCSC[A1, A1 + A1', A1'A1, transpose(A1) + A1]
                solve!(ps, X, A, B)
                @test X ≈ A\Matrix(B)

                X = solve(ps, A, B)
                @test X ≈ A\Matrix(B)

                solve!(ps, X, A, B, :C)
                @test X ≈ A'\Matrix(B)

                X = solve(ps, A, B, :C)
                @test X ≈ A'\Matrix(B)

                solve!(ps, X, A, B, :T)
                @test X ≈ copy(transpose(A))\Matrix(B)

                X = solve(ps, A, B, :T)
                @test X ≈ copy(transpose(A))\Matrix(B)
            end
        end
    end
end
end #testset

include("../examples/examplesym.jl")
include("../examples/exampleunsym.jl")
include("../examples/exampleherm.jl")
for solver in available_solvers
    example_symmetric(solver)
    example_unsymmetric(solver)
    example_hermitian_psd(solver)
end

if Pardiso.mkl_is_available()
    if Sys.CPU_THREADS >= 4
        @testset "procs" begin
            ps = MKLPardisoSolver()
            np = get_nprocs(ps)
            set_nprocs!(ps, 2)
            @test get_nprocs(ps) == 2
            set_nprocs!(ps, np)
            @test get_nprocs(ps) == np
        end
    end
end

if Pardiso.PARDISO_LOADED[]
@testset "schur" begin
    # reproduce example from Pardiso website
    include("schur_matrix_def.jl")
    @test norm(real(D) - real(C)*rA⁻¹*real(B) - s) < 1e-10*(8)^2
    # @test norm(D - C*A⁻¹*B - S) < 1e-10*(8)^2

    # try some random matrices
    m = 50; n = 15; p = .1
    ps = PardisoSolver()
    for T in (Float64, )#ComplexF64)
        ps = PardisoSolver()
        pardisoinit(ps)
        if T == Float64
            set_matrixtype!(ps, 11)
        else
            set_matrixtype!(ps, 13)
        end
        for j ∈ 1:100
            A = 5I + sprand(T,m,m,p)
            A⁻¹ = inv(Matrix(A))
            B = sprand(T,m,n,p)
            C = sprand(T,n,m,p)
            D = 5I + sprand(T,n,n,p)
            M = [A B; C D]

            # test integer block specification
            S = schur_complement(ps, M, n);
            @test norm(D - C*A⁻¹*B - S) < 1e-10*(m+n)^2

            # test sparse vector block specification
            x = spzeros(T,m+n)
            x[(m+1):(m+n)] .= 1
            S = schur_complement(ps, M, x);
            @test norm(D - C*A⁻¹*B - S) < 1e-10*(m+n)^2

            # test sparse matrix block specification
            x = spzeros(T,m+n,2)
            x[(m+1):(m+n-1),1] .= 1
            x[end,2] = 1
            S = schur_complement(ps, M, x);
            @test norm(D - C*A⁻¹*B - S) < 1e-10*(m+n)^2
        end
    end
end # testset
end

@testset "error checks" begin
for pardiso_type in available_solvers

    ps = pardiso_type()

    A = sparse(rand(10,10))
    B = rand(10, 2)
    X = rand(10, 2)

    if pardiso_type == PardisoSolver
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
for pardiso_type in available_solvers
    ps = pardiso_type()
    set_iparm!(ps, 1, 0)
    pardisoinit(ps)
    @test get_iparm(ps, 1) == 1

    @test_throws ArgumentError set_phase!(ps, 5)
    @test_throws ArgumentError set_msglvl!(ps, 2)
    @test_throws ArgumentError set_matrixtype!(ps, 15)

    if pardiso_type == PardisoSolver
        @test_throws ArgumentError set_solver!(ps, 2)

        set_dparm!(ps, 5, 13.37)
        @test get_dparm(ps, 5) == 13.37

        set_solver!(ps, 1)
        @test Int(get_solver(ps)) == 1
    end

    set_iparm!(ps, 13, 100)
    @test get_iparm(ps, 13) == 100

    set_matrixtype!(ps, Pardiso.REAL_SYM)
    @test get_matrixtype(ps) == Pardiso.REAL_SYM

    set_phase!(ps, Pardiso.ANALYSIS_NUM_FACT)
    @test get_phase(ps) == Pardiso.ANALYSIS_NUM_FACT

    set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
    @test get_msglvl(ps) == Pardiso.MESSAGE_LEVEL_ON
end

@testset "pardiso" begin
    for pardiso_type in available_solvers
        A = sparse(rand(2,2) + im * rand(2,2))
        b = rand(2)          + im * rand(2)
        ps = pardiso_type()
        set_matrixtype!(ps, Pardiso.COMPLEX_NONSYM)
        x = Pardiso.solve(ps, A, b);
        set_phase!(ps, Pardiso.RELEASE_ALL)
        pardiso(ps)
    end
end
end # testset
