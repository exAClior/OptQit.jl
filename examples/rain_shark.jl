using Pkg;
Pkg.activate(dirname(@__FILE__));
using OptQit, Yao, MosekTools, Convex

# %% [markdown]
# The first step is to declare variables 
# The order of the arguments is the size, the sign, and then the Convex.VarType (i.e., integer, binary, or continuous), and any may be omitted to use the default.
size = (2,1)
sign = Positive() # or Negative()
vartype = IntVar # or BinVar
@show y = Variable(size,sign,vartype)

z = Semidefinite(2)

cx1 = ComplexVariable()

x1 = Variable(1)
x2 = Variable()



# %% [markdown]
# Now that we have variables, we could form expressions with them




