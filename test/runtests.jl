using Pardiso
using Base.Test

srand(1234)

let
    set_mtype(11)
    set_phase(13)
    set_solver(0)

    Af = rand(10, 10)
    A = sparse(Af)
    b = rand(10, 2)
    x = zeros(10, 2)

    check_matrix(A, b)
    check_vec(b)
    print_stats(x, A, b)

    pA_ldiv_B!(x, A, b)
    @test_approx_eq x A\b
    fill!(x, 0.0)

    x = pA_ldiv_B(A, b)
    @test_approx_eq x A\b
    fill!(x, 0.0)

    pAt_ldiv_B!(x, A, b)
    @test_approx_eq x A'\b
    fill!(x, 0.0)

    x = pAt_ldiv_B(A, b)
    @test_approx_eq x A'\b
    fill!(x, 0.0)

    set_mtype(13)
    @test_throws ErrorException pardiso(x, A, b)
end

let
    set_mtype(13)
    set_phase(13)
    set_solver(0)
    Af = rand(Complex128, 10, 10)
    A = sparse(Af)
    b = rand(Complex128, 10, 2)
    x = zeros(Complex128, 10, 2)

    check_matrix(A, b)
    check_vec(b)
    print_stats(x, A, b)

    pA_ldiv_B!(x, A, b)
    @test_approx_eq x A\b
    fill!(x, zero(Complex128))

    x = pA_ldiv_B(A, b)
    @test_approx_eq x A\b
    fill!(x, zero(Complex128))

    pAt_ldiv_B!(x, A, b)
    @test_approx_eq x A.'\b
    fill!(x, zero(Complex128))

    x = pAt_ldiv_B(A, b)
    @test_approx_eq x A.'\b
    fill!(x, zero(Complex128))

    set_mtype(11)
    @test_throws ErrorException pardiso(x, A, b)
end