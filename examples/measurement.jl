# %% [markdown]
# # Quantum Marginal Problem

# %%
using Pkg;
Pkg.activate(dirname(@__FILE__));
using OptQit, Yao, MosekTools, Convex
using LinearAlgebra

"""
Given two mixed states ρ and σ, find the measurement that maximizes the
measurement result difference between the two states. They should be equal to
the trace distance between the two states.
"""
function sdp_measurement(
    ρ::AbstractMatrix, σ::AbstractMatrix; optimizer=Mosek.Optimizer, silent=true
)
    M = ComplexVariable(size(ρ)...)

    constraints = [Diagonal(ones(size(ρ, 1))) ⪰ M, M in :SDP]

    p = maximize(real(tr(M * (ρ - σ))), constraints)

    solve!(p, optimizer; silent_solver=silent)
    p.status == Convex.MathOptInterface.OPTIMAL && return evaluate(M), p.optval
    return error("Failed to find the optimal measurement.")
end

sdp_measurement(state(rand_density_matrix(2)), state(rand_density_matrix(2)))
