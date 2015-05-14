using Pardiso
using Base.Test

srand(1234)

#let
    ps = PardisoSolver()
    set_mtype(ps, 11)
    set_solver(ps, 0)

    Af = rand(10, 10)
    A = sparse(Af)
    b = rand(10, 2)
    x = zeros(10, 2)

    solve!(ps, x, A, b)
    @test_approx_eq x A\b
    fill!(x, 0.0)

    x = solve(ps, A, b)
    @test_approx_eq x A\b
    fill!(x, 0.0)

    solve!(ps, x, A, b, :T)
    @test_approx_eq x A'\b
    fill!(x, 0.0)

    x = solve(ps, A, b, :T)
    @test_approx_eq x A'\b
    fill!(x, 0.0)

    set_mtype(ps, 13)
    @test_throws ErrorException pardiso(ps, x, A, b)
#end

let
    ps = PardisoSolver()
    set_mtype(ps, 13)
    set_solver(ps, 0)

    Af = rand(Complex128, 10, 10)
    A = sparse(Af)
    b = rand(Complex128, 10, 2)
    x = zeros(Complex128, 10, 2)

     solve!(ps, x, A, b)
    @test_approx_eq x A\b
    fill!(x, zero(Complex128))

    x = solve(ps, A, b)
    @test_approx_eq x A\b
    fill!(x, zero(Complex128))

    solve!(ps, x, A, b, :T)
    @test_approx_eq x A.'\b
    fill!(x, zero(Complex128))

    x = solve(ps, A, b, :T)
    @test_approx_eq x A.'\b
    fill!(x, zero(Complex128))

    set_mtype(ps, 11)
    @test_throws ErrorException pardiso(ps, x, A, b)
end
