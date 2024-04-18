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

using Pkg;
Pkg.activate(dirname(@__FILE__));
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

# %% [markdown]
# # Convex Optimization
# Optimization is useful in life. An optimization problem can be defined  
# with the following components:
# 1. Configuration space for variables
# 2. Objective function
# 3. Constraints on the variables
# The optimization problem can be put into different categories
# when these components follows different constraints.
# For example, the problem is called a linear programming problem
# if the objective function and the constraints are linear.
# Similarily, the problem is called a convex optimization problem if 
# 1. Objective is convex
# 2. Equality constraint is linear (we could replace the equality constraint with two inequality constraints. They differ by a sign which requires the function to be both convex and concave by the next point. The only function that is both convex and concave is the linear function.) 
# 3. Inequality constraint is convex
# If you could expression an optimization problem as a convex optimization problem you are in luck. Because there are efficient algorithms to solve convex optimization problems and the solution is guaranteed to be the global minimum.
# For an optimization problem with constraint, we could always map it to an
# optimization problem in higher dimension without constraint. 
# For example, we just add a penalty function to the objective function. The penalty function is a function that is zero when the constraint is satisfied and infinity when the constraint is violated.(Think Lagrange multiplier)
# More visually, given the problem of finding the spot of no rain given that the rain cloud is moving in the (1,1) direction. And there are sharks in the water. The constraint is that you cannot swim in the water. The penalty function is the shark. The optimization problem is to find the spot of no rain that is closest to the rain cloud and furthest from the shark. We could elevate the problem from 2D to 3D. The third dimension is the objective function modified after the penalty function. The two problems will have the same solution due to the principle of strong duality.(Need more reading)

# Now that we converted constrained optimization to an unconstrained optimization and maintained the convexity of the problem. We could use the following algorithms to solve the problem:
# The dual definition of convexity for the functions mentioned above are what 
# made possible of finding the global minimum of the function from local information. For example, by definition of a convex function, it needs to be above the tangent plane at all points. When you find a tanget plane that is horizontal, you have found the global minimum. You can find such minimum by either solving the problem of $\nabla f = 0$ or by using newton's method. 

# ## Necessity of Displined Convex Programming
# In the section above, we have see the reason why we like convex optimization.
# For an optimization problem to be convex, we need to guarantee that the constraint and the objective function are convex. This is done via DCP. It is basically a derivation rule that tells you from the convexity of component of a function whether the function is convex or not. And there are function that is neither convex or concave i.e sin. 
# ## Reference
# 1. [What is Mathematical Optimization](https://www.youtube.com/watch?v=AM6BY4btj-M) is a set of videos that explains basics of convex optimization
