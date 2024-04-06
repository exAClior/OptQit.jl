# %% 
using Pkg;
Pkg.activate(dirname(@__FILE__));
using OptQit, Yao, MosekTools, Convex
using LinearAlgebra
using Test
"""
Use Semidefinite Programming to find the Choi representation of a channel.
The channel maps ρ1 to σ1 and ρ2 to σ2.

Following https://shuvomoy.github.io/blogs/posts/Solving_semidefinite_programming_problems_in_Julia/

# TODO
This is the naive version without considering c1 and c2 which
always gives a solution and the validity of the solution is 
guaranteed when only c1 = 1.0 and c2 = 0.0.
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

    A_qubits = Yao.log2i(N_A)
    B_qubits = Yao.log2i(N_B)

    # this is a CP map from A'A to A'B
    J1 = ComplexVariable(N_A * N_B, N_A * N_B)
    # c1 = Variable()
    # c2 = Variable()

    constraints = [
        J1 ⪰ 0,
        partialtrace(J1, 1, [N_A, N_B]) == LinearAlgebra.I(N_B),
        partialtrace(J1 * kron(ρ1', LinearAlgebra.I(N_B)), 1, [N_A, N_B]) == σ1,
        partialtrace(J1 * kron(ρ2', LinearAlgebra.I(N_B)), 1, [N_A, N_B]) == σ2,
    ]

    # p = minimize(c1+c2, constraints)
    p = satisfy(constraints)

    solve!(p, optimizer; silent_solver=silent)
    @show typeof(p.status), p.status
    p.status != :Optimal || return error("SDP failed to find a solution")
    # (c1 == 1.0 && c2 == 0.0) || return error("Invalid solution")

    # return evaluate(J1), evaluate(c1), evaluate(c2)
    return evaluate(J1)
end

@testset "SDP Choi trivial" begin
    num_qubits = 1

    ρ1 = rand_density_matrix(num_qubits).state
    ρ2 = rand_density_matrix(num_qubits).state

    J1 = sdp_Choi_rep(ρ1, ρ2, ρ1, ρ2)

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
