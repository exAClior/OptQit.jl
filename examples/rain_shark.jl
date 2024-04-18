using Pkg;
Pkg.activate(dirname(@__FILE__));
using OptQit, Yao, MosekTools, Convex

# %% [markdown]
# Let's follow an example.
# Pirate on the sea, enemy ship firing from direction of (1,1)
# Sharks are surrounding you in the water, they trace out a circle of radius 1.
# The sea is quite small, so you can only move to a point (x1, x2) such that $sqrt(1+x1^2) + x2^2 \leq 100$
# What is the best location you can be?

# The first step is to declare variables 
# The order of the arguments is the size, the sign, and then the Convex.VarType (i.e., integer, binary, or continuous), and any may be omitted to use the default.
size = (2,1)
sign = Positive() # or Negative(), NoSign(), ComplexSign()
vartype = IntVar # or BinVar, ContVar

@show y = Variable(size,sign,vartype)

z = Semidefinite(2) # alias for  Variable((m, m), NoSign(), ContVar)
z = ComplexVariable(size) # alias for Variable(size,ComplexSign(),ContVar)

# %% [markdown]
# Ofcourse, it is easy to define your own variables. We won't go over it here. You may find examples [here](https://jump.dev/Convex.jl/stable/advanced/#Custom-Variable-Types)

x1,x2 = Variable(2)


# %% [markdown]
# The next step is to define the constraints and the objective function
# Notice, the special way of building expression from variables 
# Go back to slides now.

arr = Matrix{Float64}[]
push!(arr, rand(2,2))
push!(arr, rand(4,4))

constraints = [
    sqrt(1 + square(x1)) + square(x2) <= 1
    square(x1) + square(x2) <= 1
]

objective = minimize(x1 + x2, constraints)

solve!(objective, Mosek.Optimizer)

objective.status
objective.optval
evaluate(x1)
evaluate(x2)





