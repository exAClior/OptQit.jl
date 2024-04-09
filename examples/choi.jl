# %% 
using Pkg;
Pkg.activate(dirname(@__FILE__));
using OptQit, Yao, MosekTools, Convex
using LinearAlgebra
using Test
"""
Use Semidefinite Programming to find the Choi representation of a channel.
The channel maps ρ1 to σ1 and ρ2 to σ2.
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
    J2 = ComplexVariable(N_A * N_B, N_A * N_B)
    Jn = J1 - J2
    c1 = Variable()
    c2 = Variable()

    constraints = [
        J1 ⪰ 0,
        J2 ⪰ 0,
        c1 >= 0,
        c2 >= 0,
        partialtrace(J1, 2, [N_A, N_B]) == c1 * LinearAlgebra.I(N_A),
        partialtrace(J2, 2, [N_A, N_B]) == c2 * LinearAlgebra.I(N_A),
        partialtrace(Jn * kron(ρ1', LinearAlgebra.I(N_B)), 1, [N_A, N_B]) == σ1,
        partialtrace(Jn * kron(ρ2', LinearAlgebra.I(N_B)), 1, [N_A, N_B]) == σ2,
    ]

    obj = minimize(c1 + c2, constraints)

    solve!(obj, optimizer; silent_solver=silent)
    return evaluate(J1), evaluate(J2), evaluate(c1), evaluate(c2)
end

@testset "SDP Choi trivial" begin
    num_qubits = 1

    ρ1 = state(rand_density_matrix(num_qubits))
    ρ2 = state(rand_density_matrix(num_qubits))

    J1, J2, c1, c2 = sdp_Choi_rep(ρ1, ρ2, ρ1, ρ2; optimizer=SCS.Optimizer)

    c1
    @test isapprox(c1, 1.0, atol=1e-3)
    @test isapprox(c2, 0.0, atol=1e-3)

    @test isapprox(
        partial_tr(
            DensityMatrix(J1 * kron(ρ1', mat(igate(num_qubits)))),
            tuple((num_qubits + 1):(num_qubits + num_qubits)...),
        ).state,
        ρ1,
        atol=1e-3,
    )

    @test isapprox(
        partial_tr(
            DensityMatrix(J1 * kron(ρ2', mat(igate(num_qubits)))),
            tuple((num_qubits + 1):(num_qubits + num_qubits)...),
        ).state,
        ρ2,
        atol=1e-3,
    )
end

@testset "SDP Choi non-trivial" begin
    num_qubits = 1

    ρ1 = state(density_matrix(zero_state(num_qubits)))
    ρ2 = state(density_matrix(ghz_state(num_qubits)))

    σ1 = state(density_matrix(arrayreg(bit"1")))
    σ2 = state(density_matrix((arrayreg(bit"0") - arrayreg(bit"1")) / sqrt(2)))

    J1, J2, c1, c2 = sdp_Choi_rep(ρ1, ρ2, σ1, σ2)

    @test c1 ≈ 1.0 atol = 1e-8
    @test c2 ≈ 0.0 atol = 1e-8

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
