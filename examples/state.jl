# %% [markdown]
# # Introduction
# In this notebook, we convert the contents in book into code.

# %% 
using Pkg;
Pkg.activate(dirname(@__FILE__));
using OptQit, Yao, MosekTools, Convex

# %% [markdown]
# ## Quantum State Estimation

# ### Example 1: State Reconstruction from Measurement result
# Given a set of measurement and measurement results,
# we try to reconstruct the state of the system. 

# %%

function eg1()
    optimizer = Mosek.Optimizer
    n_qubits = 2

    # secret state, we cheat a little to provide it here
    σ = state(rand_density_matrix(n_qubits))

    ρ = ComplexVariable(2^n_qubits, 2^n_qubits)

    num_meas = 5
    meas_ops = [mat(kron(rand([X,Y,Z],n_qubits)...)) for _ in 1:num_meas]
    meas_vals = [tr(σ*meas_ops[i]) for i in 1:num_meas] 

    constrains = [
        ρ in :SDP,
        tr(ρ) == 1,
        [tr(M_i * ρ) == m_i for (M_i, m_i) in zip(meas_ops, meas_vals)]...,
    ]

    p = satisfy(constrains)
    solve!(p, optimizer; silent_solver=true)

    @show p
    res = evaluate(ρ)
    @show res 
    @assert all([tr(res*meas_ops[i]) ≈ meas_vals[i] for i in 1:num_meas])
end

eg1()

# ### Example 2: Trace Distance Estimation

function eg2()
    optimizer = Mosek.Optimizer
    n_qubits = 2

    ρ = ComplexVariable(2^n_qubits, 2^n_qubits)
    σ = state(rand_density_matrix(n_qubits))

    num_meas = 5
    meas_ops = [rand_hermitian(2^n_qubits) for _ in 1:num_meas]
    meas_vals = [tr(σ*meas_ops[i])+rand(0:0.001) for i in 1:num_meas] 

    constrains = [
        ρ ⪰ 0,
        tr(ρ) == 1,
        [tr(M_i * ρ) == m_i for (M_i, m_i) in zip(meas_ops, meas_vals)]...,
    ]


    p = minimize(nuclearnorm(ρ - σ)/2.0, constrains)
    solve!(p, optimizer; silent_solver=true)

    @show p.status
    @show p.optval
    res = evaluate(ρ)
    @show  res 
    @assert res ≈ σ 
end


function eg3()
    optimizer = Mosek.Optimizer
    n_qubits = 2

    ρ = ComplexVariable(2^n_qubits, 2^n_qubits)
    op_X = ComplexVariable(2^n_qubits, 2^n_qubits)
    σ = state(rand_density_matrix(n_qubits))

    num_meas = 5
    meas_ops = [rand_hermitian(2^n_qubits) for _ in 1:num_meas]
    meas_vals = [tr(σ*meas_ops[i])+rand(0:0.001) for i in 1:num_meas] 

    constrains = [
        ρ ⪰ 0,
        tr(ρ) == 1,
        [tr(M_i * ρ) == m_i for (M_i, m_i) in zip(meas_ops, meas_vals)]...,
        ρ - σ ⪰ - op_X,
        ρ - σ ⪯ op_X,
    ]


    p = minimize(real(tr(op_X)/2.0), constrains)
    solve!(p, optimizer; silent_solver=true)

    @show p.status
    @show p.optval
    res = evaluate(ρ)
    @show  res 
    @assert res ≈ σ 
end

# %% [markdown]
# ## Quantum Magrinal Problem

# %% 
function eg4()
    sub_system_dim = 2 
    num_sys = 3
    epsilon = 0.0001
    σ = state(rand_density_matrix(num_sys, nlevel=sub_system_dim))

    ρ = ComplexVariable(sub_system_dim^num_sys, sub_system_dim^num_sys)
    
    ρXY = partialtrace(σ,3, repeat([sub_system_dim], num_sys) )
    ρXZ = partialtrace(σ,2, repeat([sub_system_dim], num_sys) )
    ρYZ = partialtrace(σ,1, repeat([sub_system_dim], num_sys) )


    constrains = [
        ρ ⪰ 0,
        tr(ρ) == 1,
        nuclearnorm(partialtrace(ρ,3,repeat([sub_system_dim],num_sys)) - ρXY)/2.0 <= epsilon ,
        nuclearnorm(partialtrace(ρ,2,repeat([sub_system_dim],num_sys)) - ρXZ)/2.0 <= epsilon,
        nuclearnorm(partialtrace(ρ,1,repeat([sub_system_dim],num_sys)) - ρYZ) /2.0 <= epsilon,
    ]

    p = satisfy(constrains)

    solve!(p, Mosek.Optimizer; silent_solver=true)

    @show p.status
    res = evaluate(ρ)
    @show res
    @show res .- σ
end


eg4()

# %% [markdown]
# # Exercises