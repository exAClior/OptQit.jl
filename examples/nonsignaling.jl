# %%
using Pkg;
Pkg.activate(dirname(@__FILE__));
using OptQit, Yao, MosekTools, Convex
using LinearAlgebra

# operator must be a quantum channel
function choi_of_op(op::AbstractMatrix)
    dim = size(op, 1)
    res = zeros(dim * dim, dim * dim)
    for i in 1:dim, j in 1:dim
        ρij = zeros(dim, dim)
        ρij[i, j] = 1
        res .+= kron(ρij, op * ρij * op')
    end
    return res
end

# %% [markdown]
# We can represent the 4 possible functions with 1 bit input and 1 bit output as Choi matrix representation of classical-quantum channel

# %% 

# %% [markdown]
# We denote whether a function is balanced or not by a single qubit. In the quantum circuit that implements the algorithm to test if a function is balanced or not, the output qubit will be in state $|0\rangle\langle 0|$ if the function is constant and $|1\rangle\langle 1|$ if the function is balanced. Hence, the target function we try to maximize is

# %%
function target(Mis)
    s0 = state(density_matrix(product_state(bit"0")))
    s1 = state(density_matrix(product_state(bit"1")))
    return (
        tr(s0 * partialtrace(Mis[1] * kron(s0, I(2)), 1, [2, 2])) +
        tr(s1 * partialtrace(Mis[2] * kron(s0, I(2)), 1, [2, 2])) +
        tr(s1 * partialtrace(Mis[3] * kron(s0, I(2)), 1, [2, 2])) +
        tr(s0 * partialtrace(Mis[4] * kron(s0, I(2)), 1, [2, 2]))
    ) / 4.0
end

function deutsch(
    n_a::Int, n_b::Int, optimizer=Mosek.Optimizer, silent=true, signaling=false
)
    n_a, n_b = 1, 1
    optimizer = Mosek.Optimizer
    silent = true
    fs = [
        [1 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0],
        [0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 1],
        [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1],
        [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0],
    ]
    # order is A_i, B_i, A_o, B_o
    J_pi = ComplexVariable(2^n_a * 2^n_b * 2^n_a * 2^n_b, 2^n_a * 2^n_b * 2^n_a * 2^n_b)

    constraints = [
        J_pi in :SDP,
        partialtrace(J_pi, [3, 4], [2^n_a, 2^n_b, 2^n_a, 2^n_b]) ==
        LinearAlgebra.I(2^n_a * 2^n_b),
        partialtrace(J_pi, 4, [2^n_a, 2^n_b, 2^n_a, 2^n_b]) ==
        partialtrace(J_pi, 3, [2^n_a, 2^n_b, 2^n_a, 2^n_b]),
        partialtrace(J_pi, 3, [2^n_a, 2^n_b, 2^n_a, 2^n_b]) == kron(
            partialtrace(J_pi, [2, 4], [2^n_a, 2^n_b, 2^n_a, 2^n_b]),
            LinearAlgebra.I(2^n_b) ./ n_b,
        ),
    ]

    Mis = [
        partialtrace(
            J_pi * kron(transpose(fs[i]), LinearAlgebra.I(2^n_a * 2^n_b)),
            [2, 3],
            [2^n_a, 2^n_b, 2^n_a, 2^n_b],
        ) for i in 1:4
    ]

    obj = maximize(target(Mis), constraints)

    solve!(obj, optimizer; silent_solver=silent)

    return value(J_pi), obj.optval
end
