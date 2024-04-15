using Pkg;
Pkg.activate(dirname(@__FILE__));
using OptQit, Yao, MosekTools, Convex

# %% [markdown]
# The first step is to declare variables 
# The order of the arguments is the size, the sign, and then the Convex.VarType (i.e., integer, binary, or continuous), and any may be omitted to use the default.
size = (2,1)
sign = Positive() # or Negative(), NoSign(), ComplexSign()
vartype = IntVar # or BinVar, ContVar

@show y = Variable(size,sign,vartype)

z = Semidefinite(2) # alias for  Variable((m, m), NoSign(), ContVar)
z = ComplexVariable(size) # alias for Variable(size,ComplexSign(),ContVar)

# %% [markdown]
# Ofcourse, it is easy to define your own variables. For example, let's define a probability amplitude variable

# %%
mutable struct ProbabilityAmpVector <: Convex.AbstractVariable
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int, Int}
    value::Convex.ValueOrNothing
    vexity::Convex.Vexity
    function ProbabilityAmpVector(d)
        this = new(:ProbabilityAmpVector, 0, (d,1), nothing, Convex.AffineVexity())
        this.id_hash = objectid(this)
        this
    end
end

Convex.constraints(p::ProbabilityAmpVector) = [ norm(p,2) == 1 ]
Convex.sign(::ProbabilityAmpVector) = Convex.ComplexSign()
Convex.vartype(::ProbabilityAmpVector) = Convex.ContVar

(p::ProbabilityVector)(x) = dot(p, x)

x1 = Variable(1)
x2 = Variable()


# %% [markdown]
# Now that we have variables, we could form expressions with them




