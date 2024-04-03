using Convex
using Yao


ρ = ComplexVariable(4,4)

A = Matrix(rand_density_matrix(ComplexF64,2))
B = Matrix(rand_density_matrix(ComplexF64,3))

AB = kron(A,B)

@show B .- partialtrace(AB,1,[2^2,2^3])


@assert isapprox(B,partialtrace(AB,1,[2^2,2^3]))
