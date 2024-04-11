using Pkg;
Pkg.activate(@__DIR__);
using OptQit, Yao, MosekTools, Convex
using LinearAlgebra
using Test

function fj(j::Int, n::Int)
    (n - 1) == j && return 1
    return 0
end

function Uj(j::Int)
    res = zeros(ComplexF64, 8, 8)
    for n in 1:4
        y = fj(j, n)
        res[(n - 1) * 2 + y + 1, (n - 1) * 2 + 1] = 1
        res[(n - 1) * 2 + 1, (n - 1) * 2 + y + 1] = 1
    end
    for n in 1:8
        if all(res[n, :] .== 0)
            res[n, n] = 1
        end
    end
    return res
end

function deutsch_joza(
    Pjs::AbstractVector{MT},
    U_fs::AbstractVector{MT};
    optimizer=Mosek.Optimizer,
    silent=true,
) where {MT<:AbstractMatrix}
    n = 2
    res = 0.0

    # order is Aq1, Aq2, Aq3, Bq1, Bq2, Bq3
    Jd = ComplexVariable((2^(n + 1))^2, (2^(n + 1))^2)

    constraints = [
        Jd ⪰ 0,
        partialtrace(Jd, [4, 5, 6], [2, 2, 2, 2, 2, 2]) == LinearAlgebra.I(2^(n + 1)),
    ]

    for j in 0:3
        psi = zero_state(n + 1)
        apply!(psi, repeat(n + 1, H, collect(1:(n + 1))))
        apply!(psi, put(n + 1, (1:(n + 1)) => matblock(U_fs[j + 1])))
        # apply U_r
        # applya U_f again, the second query

        ρ = state(density_matrix(psi))
        fin_state = partialtrace(
            Jd * kron(transpose(ρ), LinearAlgebra.I(2^(n + 1))),
            [1, 2, 3],
            [2, 2, 2, 2, 2, 2],
        )
        res += real(tr(fin_state * Pjs[j + 1]) / 4.0)
    end

    obj = maximize(res, constraints)
    solve!(obj, optimizer; silent_solver=silent)
    return evaluate(Jd), obj.optval
end

function deutsch_joza2(
    Pjs::AbstractVector{MT},
    U_fs::AbstractVector{MT};
    optimizer=Mosek.Optimizer,
    silent=true,
) where {MT<:AbstractMatrix}
    n = 2
    res = 0.0

    # order is Aq1i, Aq2i, Aq3i,  Aq1o, Aq2o, Aq3o, Bq1i, Bq2i, Bq3i, Bq1o, Bq2o, Bq3o
    Jd = ComplexVariable((2^(n + 1))^4, (2^(n + 1))^4)

    constraints = [
        Jd ⪰ 0,
        partialtrace(Jd, [4, 5, 6], [2, 2, 2, 2, 2, 2]) == LinearAlgebra.I(2^(n + 1)),
    ]

    for j in 0:3
        psi = zero_state(n + 1)
        apply!(psi, repeat(n + 1, H, collect(1:(n + 1))))
        apply!(psi, put(n + 1, (1:(n + 1)) => matblock(U_fs[j + 1])))
        apply!(psi, put(n + 1, (1:(n + 1)) => matblock(U_rand)))
        apply!(psi, put(n + 1, (1:(n + 1)) => matblock(U_fs[j + 1])))

        ρ = state(density_matrix(psi))
        fin_state = partialtrace(
            Jd * kron(transpose(ρ), LinearAlgebra.I(2^(n + 1))),
            [1, 2, 3],
            [2, 2, 2, 2, 2, 2],
        )
        res += real(tr(fin_state * Pjs[j + 1]) / 4.0)
    end

    obj = maximize(res, constraints)
    solve!(obj, optimizer; silent_solver=silent)
    return evaluate(Jd), obj.optval
end

@testset "Deutsch Joza" begin
    Pjs = [kron(mat(projector(product_state(2, i - 1))), [1.0 0.0; 0.0 1.0]) for i in 1:4]

    U_fs = [Uj(j) for j in 1:4]

    D_channel, op_val = deutsch_joza(Pjs, U_fs; optimizer=Mosek.Optimizer, silent=true)

    op_val
    @test isapprox(op_val, 0.25, atol=1e-6)
end

@testset "Deutsch Joza2" begin
    Pjs = [kron(mat(projector(product_state(2, i - 1))), [1.0 0.0; 0.0 1.0]) for i in 1:4]

    U_fs = [Uj(j) for j in 1:4]

    D_channel, op_val = deutsch_joza2(Pjs, U_fs; optimizer=Mosek.Optimizer, silent=true)

    op_val
    @test isapprox(op_val, 0.25, atol=1e-6)
end

function O_f(n::Int)
    res = diagm(ones(ComplexF64, 4))
    res[n, n] = -1
    return res
end

function grover(
    Pjs::AbstractVector{P}, O_fs::AbstractVector{MT}; optimizer=Mosek.Optimizer, silent=true
) where {MT<:AbstractMatrix,P<:PrimitiveBlock}
    n = 2
    # we are modeling two qubit circuit here
    # order is Aq1, Aq2, Bq1, Bq2 
    Jd = ComplexVariable((2^n)^2, (2^n)^2)

    constraints = [Jd ⪰ 0, partialtrace(Jd, [3, 4], [2, 2, 2, 2]) == LinearAlgebra.I(2^n)]

    res = 0.0
    for j in 0:3
        psi = zero_state(n)
        apply!(psi, repeat(n, H, collect(1:n)))
        apply!(psi, put(n, (1:n) => matblock(O_fs[j + 1])))
        ρ = state(density_matrix(psi))
        fin_state = partialtrace(
            Jd * kron(transpose(ρ), LinearAlgebra.I(2^n)), [1, 2], [2, 2, 2, 2]
        )
        res += real(tr(fin_state * mat(Pjs[j + 1])) / 4.0)
    end
    obj = maximize(res, constraints)
    solve!(obj, optimizer; silent_solver=silent)
    return evaluate(Jd), obj.optval
end

@testset "Grover" begin
    Pjs = [projector(product_state(2, i - 1)) for i in 1:4]

    O_fs = [O_f(i) for i in 1:4]

    D_channel, op_val = grover(Pjs, O_fs; optimizer=Mosek.Optimizer, silent=true)

    D_channel
    @test isapprox(op_val, 1.0, atol=1e-6)
end
