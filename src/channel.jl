using Convex, SCS, LinearAlgebra 

function ispsd(A::AbstractMatrix; tol::Float64=eps(Float64)^(3/4))
    size(A, 1) == size(A, 2) || throw(DimensionMismatch("Matrix must be square"))
    A == A' || return false 
    L, = cholesky(A+(tol+eps(Float64))*I)
    return all(L.L .> 0)
end

function ispsd_wit(A::AbstractMatrix; tol::Float64=eps(Float64)^(3/4))
    size(A, 1) == size(A, 2) || throw(DimensionMismatch("Matrix must be square"))
    T = (A + A')/2
    evals,evecs = eigen(T)
    @show minimum(evals), evals[1]
    return (minimum(evals) >= -tol ,evecs[:,1])
end

using Yao
A = state(rand_density_matrix(3))
A = round.(A, digits=3)
A = A .+ rand(ComplexF64)*1e-3

A = rand(8,8)

@time isposdef(A)

@time yn, wit = ispsd_wit(A)
@time ispsd(A)
