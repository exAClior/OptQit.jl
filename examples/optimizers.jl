# %% [markdown]
# CVX Supports a variety of optimization solvers. Some have julia wrapper. We document those here.

# [Decision tree for Optimizer](https://plato.asu.edu/guide.html)

# 1. [SDPT3](https://jump.dev/JuMP.jl/stable/packages/SDPT3/) requires matlab
# 2. [SEDUMI](https://jump.dev/JuMP.jl/stable/packages/SeDuMi/) requires Matlab
# 3. [Gurobi](https://jump.dev/JuMP.jl/stable/packages/Gurobi/) requires [Gurobi](https://www.gurobi.com/academia/academic-program-and-licenses/) license, very annoying
# 4. [Mosek](https://github.com/MOSEK/Mosek.jl) easiest to [install](https://docs.mosek.com/latest/install/installation.html) and [acquire license](https://docs.mosek.com/10.1/licensing/quickstart.html#i-don-t-have-a-license-file-yet), use this!  [Documentation](https://docs.mosek.com/latest/install/installation.html)

# For the full list of supported solvers, see [Jump.jl documentation](https://jump.dev/JuMP.jl/stable/packages/solvers/). To list a few ones that is Plug and Play

# 1. [Ipopt](https://jump.dev/JuMP.jl/stable/packages/Ipopt/) 

# If you are looking to implement your own solver in Julia. You might find this [paper](http://euler.nmt.edu/~brian/csdppaper.pdf) useful. 

using Mosek, MosekTools
using Convex, SCS

A = [0.47213595 0.11469794+0.48586827im; 0.11469794-0.48586827im 0.52786405]
B = ComplexVariable(2, 2)
ρ = kron(A, B)

constraints = [
    partialtrace(ρ, 1, [2; 2]) == [1 0; 0 0]
    tr(ρ) == 1
    ρ in :SDP
]
p = satisfy(constraints)
@time solve!(p, Mosek.Optimizer; silent_solver=true)
@time solve!(p, SCS.Optimizer; silent_solver=true)
p.status
