using Pkg;
Pkg.activate(dirname(@__FILE__));
using OptQit, Yao, MosekTools, Convex
using LinearAlgebra
using Test

function fj(j::Int, n::Int)
    (n - 1) == j && return 1
    return 0
end

function Uj(j::Int)
    res = zeros(8, 8)
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

function deutsch_joza_yao(n::Int, U_f, D)
    res = ComplexF64(0.0)
    Pjs = Matrix{ComplexF64}[]

    for k in 1:4
        pj = zeros(ComplexF64, 4)
        pj[k] = 1
        push!(Pjs, kron(pj * pj', mat(I2)))
    end

    for j in 0:3
        psi = zero_state(n + 1)
        apply!(psi, repeat(n + 1, H, collect(1:(n + 1))))
        apply!(psi, put(n + 1, (1:3) => matblock(U_f)))
        apply!(psi, put(n + 1, (1:3) => matblock(D)))
        res += tr(state(density_matrix(psi)) * Pjs[j + 1]) / 4.0
    end
    return res
end

deutsch_joza_yao(2, Uj(0), rand_unitary(8))

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
