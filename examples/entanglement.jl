# %% [markdown]
# # Quantum Game
# [Paper](https://arxiv.org/pdf/1303.2849.pdf)

# %%
using Pkg;
Pkg.activate(dirname(@__FILE__));
using OptQit, Yao, MosekTools, Convex
using LinearAlgebra

function cabxy(a, b, x, y)
    return (-1)^(a + b + x * y)
end

# function cabxy(a, b, x, y)
# 	return (-1)^(a * x + b * y)
# end

# function cabxy(a, b, x, y)
# 	return (0.5)^(a + b + x * y) + (-0.5)^(x + y + a * b)
# end

# function cabxy(a, b, x, y)
# 	return exp(a)^b^x^y 
# end

function target1(As, Bs, ρ, score::Function)
    return real(
        sum([
            score(a, b, x, y) * tr(kron(As[(a) * 2 + x + 1], Bs[(b) * 2 + y + 1]) * ρ) for
            a in 0:1, b in 0:1, x in 0:1, y in 0:1
        ]),
    )
end

function chsh_op_state(
    target::Function,
    As::AbstractVector{MT},
    Bs::AbstractVector{MT},
    restriction::Expr,
    nlevels::Pair{Int,Int};
    optimizer=Mosek.Optimizer,
    silent=true,
) where {MT<:AbstractMatrix}
    dim = reduce(*, nlevels)
    ρ = ComplexVariable(dim, dim)

    constriants = [ρ in :SDP, tr(ρ) == 1.0]

    if restriction == :PPT
        push!(constriants, partialtranspose(ρ, 1, nlevels) in :SDP)
    end

    p = maximize(target(As, Bs, ρ), constraints)

    solve!(p, optimizer; silent_solver=silent)
    p.status == Convex.MathOptInterface.OPTIMAL && return evaluate(ρ), p.optval
    return error("Failed to find the optimal State.")
end

function chsh_op_povm(
    target::String,
    As,
    Bs,
    ρ,
    nlevels::Pair{Int,Int};
    optimizer=Mosek.Optimizer,
    silent=true,
)
    v0 = rand_state(1; nlevel=n)
    v1 = rand_state(1; nlevel=n)

    if isnothing(Bs)
        Bs = Matrix{ComplexF64}[]
        B00 = v0.state * v0.state'
        B10 = I - B00
        B01 = v1.state * v1.state'
        B11 = I - B01
        Bs = [B00, B01, B10, B11]
    end

    if isnothing(As)
        target == "ρ" && error("whoa")
        As = [@variable(model, [1:n, 1:n] in HermitianPSDCone()) for _ in 1:4]
        @constraint(model, As[1] + As[3] == LinearAlgebra.I)
        @constraint(model, As[2] + As[4] == LinearAlgebra.I)
    end

    @objective(model, Max,)

    solve!(p, optimizer; silent_solver=silent)

    @assert p.status == Convex.MathOptInterface.OPTIMAL

    if target == "ρ"
        return objective_value(model), value.(ρ)
    else
        return objective_value(model), [value.(As[ii]) for ii in 1:4]
    end
end

@testset "entanglement" begin
    dit = 2
    # ρ_init = density_matrix(rand_state(2; nlevel = dit)).state
    ρ_init = [
        0.3943 0.0348 0.0348 0.3282
        0.0348 0.1057 0.0303 -0.0348
        0.0348 0.0303 0.1057 -0.0348
        0.3282 -0.0348 -0.0348 0.3943
    ]

    As = nothing
    Bs = nothing
    for jj in 1:100
        obj_val1, As = chsh_op_povm("As", nothing, nothing, ρ_init, dit)
        obj_val2, Bs = chsh_op_povm("As", nothing, As, ρ_init, dit)
        # obj_val3, ρ = chsh_op_povm("ρ", As, Bs, nothing, dit)
        println(obj_val1, " ", obj_val2)
        # println(obj_val1, " ", obj_val2, " ", obj_val3)
    end

    @test isapprox(obj_val2, 2 * sqrt(2), atol=1e-4)
end

@testset "entanglement" begin
    dit = 3
    ρ_init = density_matrix(rand_state(2; nlevel=dit)).state
    for _ in 1:10
        As = nothing
        Bs = nothing
        obj_val1, As = chsh_op_povm("As", nothing, nothing, ρ_init, dit)
        obj_val2, Bs = chsh_op_povm("As", nothing, As, ρ_init, dit)
        obj_val3, ρ = chsh_op_povm("ρ", As, Bs, nothing, dit)
        println(obj_val1, " ", obj_val2, " ", obj_val3)
    end

    @test isapprox(obj_val2, 2 * sqrt(2), atol=1e-4)
end

@testset "entanglement" begin
    dit = 5
    ρ_init = density_matrix(rand_state(2; nlevel=dit)).state
    for _ in 1:10
        As = nothing
        Bs = nothing
        obj_val1, As = chsh_op_povm("As", nothing, nothing, ρ_init, dit)
        obj_val2, Bs = chsh_op_povm("As", nothing, As, ρ_init, dit)
        obj_val3, ρ = chsh_op_povm("ρ", As, Bs, nothing, dit)
        println(obj_val1, " ", obj_val2, " ", obj_val3)
    end

    using BitBasis

    @test isapprox(obj_val2, 2 * sqrt(2), atol=1e-4)
end
