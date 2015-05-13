using Pardiso
using Base.Test


set_mtype(11)
set_phase(13)
set_solver(0)

let
    srand(1234)
    Af = rand(10, 10)
    A = sparse(Af)
    b = rand(10)
    x = zeros(10)

    check_matrix(A, b)
    check_vec(b)
    print_stats(x, A, b)

    pA_ldiv_B!(x, A, b)
    @test_approx_eq x A\b

    x = pA_ldiv_B(A, b)
    @test_approx_eq x A\b

    pAt_ldiv_B!(x, A, b)
    @test_approx_eq x A'\b

    x = pAt_ldiv_B(A, b)
    @test_approx_eq x A'\b
end

