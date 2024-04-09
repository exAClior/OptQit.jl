# %% 
using Pkg;
Pkg.activate(dirname(@__FILE__));
using OptQit, Yao, MosekTools, Convex
using SCS
using LinearAlgebra
using Test
"""
Use Semidefinite Programming to find the Choi representation of a channel.
The channel maps ρ1 to σ1 and ρ2 to σ2.

Following https://shuvomoy.github.io/blogs/posts/Solving_semidefinite_programming_problems_in_Julia/
"""
function sdp_Choi_rep(
    ρ1::AbstractMatrix,
    ρ2::AbstractMatrix,
    σ1::AbstractMatrix,
    σ2::AbstractMatrix;
    optimizer=Mosek.Optimizer,
    silent=true,
)
    N_A = size(ρ1, 1)
    N_B = size(σ1, 1)

    # this is a CP map from A'A to A'B
    J1 = ComplexVariable(N_A * N_B, N_A * N_B)

    constraints = [
        J1 ⪰ 0,
        partialtrace(J1, 2, [N_A, N_B]) == LinearAlgebra.I(N_A),
        partialtrace(J1 * kron(ρ1', LinearAlgebra.I(N_B)), 1, [N_A, N_B]) == σ1,
        partialtrace(J1 * kron(ρ2', LinearAlgebra.I(N_B)), 1, [N_A, N_B]) == σ2,
    ]

    p = satisfy(constraints)

    solve!(p, optimizer; silent_solver=silent)
    # p.status != Convex.MathOptInterface.OPTIMAL && return error("SDP failed to find a solution")
    return evaluate(J1)
end

@testset "SDP Choi trivial" begin
    num_qubits = 1

    ρ1 = rand_density_matrix(num_qubits).state
    ρ2 = rand_density_matrix(num_qubits).state

    J1 = sdp_Choi_rep(ρ1, ρ2, ρ1, ρ2; optimizer=SCS.Optimizer)

    @test isapprox(
        partial_tr(
            DensityMatrix(J1 * kron(ρ1', mat(igate(num_qubits)))),
            tuple((num_qubits + 1):(num_qubits + num_qubits)...),
        ).state,
        ρ1,
        atol=1e-6,
    )

    @test isapprox(
        partial_tr(
            DensityMatrix(J1 * kron(ρ2', mat(igate(num_qubits)))),
            tuple((num_qubits + 1):(num_qubits + num_qubits)...),
        ).state,
        ρ2,
        atol=1e-6,
    )
end

@testset "SDP Choi non-trivial" begin
    num_qubits = 1

    ρ1 = density_matrix(zero_state(num_qubits)).state
    ρ2 = density_matrix(ghz_state(num_qubits)).state

    σ1 = density_matrix(arrayreg(bit"1")).state
    σ2 = density_matrix((arrayreg(bit"0") - arrayreg(bit"1")) / sqrt(2)).state

    J1 = sdp_Choi_rep(ρ1, ρ2, σ1, σ2)

    @test isapprox(
        partial_tr(
            DensityMatrix(J1 * kron(ρ1', mat(igate(num_qubits)))),
            tuple((num_qubits + 1):(num_qubits + num_qubits)...),
        ).state,
        σ1,
        atol=1e-6,
    )

    @test isapprox(
        partial_tr(
            DensityMatrix(J1 * kron(ρ2', mat(igate(num_qubits)))),
            tuple((num_qubits + 1):(num_qubits + num_qubits)...),
        ).state,
        σ2,
        atol=1e-6,
    )
end
