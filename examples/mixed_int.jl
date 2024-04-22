using Pkg;
Pkg.activate(dirname(@__FILE__));
using OptQit, Yao, MosekTools, Convex

x = Variable(4, :Int)
p = minimize(sum(x), x >= 0.5)
solve!(p, Mosek.Optimizer; silent_solver = true)
evaluate(x)
