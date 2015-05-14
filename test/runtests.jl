using Pardiso
using Base.Test

srand(1234)

# Test solver + checkers for real matrices
let
    ps = PardisoSolver()
    pardisoinit(ps)
    set_mtype(ps, 11)
    set_solver(ps, 0)

    A = sparse(rand(10,10))
    B = rand(10, 2)
    X = zeros(10, 2)

    printstats(ps, A, B)
    checkmatrix(ps, A, B)
    checkvec(B)

    solve!(ps, X, A, B)
    @test_approx_eq X A\B
    fill!(X, 0.0)

    X = solve(ps, A, B)
    @test_approx_eq X A\B
    fill!(X, 0.0)

    solve!(ps, X, A, B, :T)
    @test_approx_eq X A'\B
    fill!(X, 0.0)

    X = solve(ps, A, B, :T)
    @test_approx_eq X A'\B
    fill!(X, 0.0)
end

# Test some errors
let
    ps = PardisoSolver()

    A = sparse(rand(10,10))
    B = rand(10, 2)
    X = rand(10, 2)

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

# Test solver + checkers for complex matrices
let
    ps = PardisoSolver()
    pardisoinit(ps)
    set_mtype(ps, 13)
    set_solver(ps, 0)

    A = sparse(rand(Complex128, 10, 10))
    B = rand(Complex128, 10, 2)
    X = zeros(Complex128, 10, 2)

    printstats(ps, A, B)
    checkmatrix(ps, A, B)
    checkvec(B)

    solve!(ps, X, A, B)
    @test_approx_eq X A\B
    fill!(X, zero(Complex128))

    X = solve(ps, A, B)
    @test_approx_eq X A\B
    fill!(X, zero(Complex128))

    solve!(ps, X, A, B, :T)
    @test_approx_eq X A.'\B
    fill!(X, zero(Complex128))

    X = solve(ps, A, B, :T)
    @test_approx_eq X A.'\B
    fill!(X, zero(Complex128))

    set_mtype(ps, 11)
    @test_throws ErrorException pardiso(ps, X, A, B)
end

let
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