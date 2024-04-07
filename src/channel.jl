# is this necessary?
function ispsd(A::AbstractMatrix{T}; tol::Real=eps(real(T))^(3 / 4)) where T    
    return isposdef(A + (tol + eps(tol)) * I)
end


"""
    ispsd_wit(A::AbstractMatrix; tol::Float64=eps(Float64)^(3/4))

Check if a matrix `A` is positive semidefinite (PSD) eigenvalues.

## Arguments
- `A::AbstractMatrix`: The matrix to be checked.
- `tol::Float64`: Tolerance value for eigenvalue comparison. Default is `eps(Float64)^(3/4)`.

## Returns
- `(is_psd::Bool, eigenvector::AbstractVector)`: A tuple containing a boolean value indicating whether `A` is PSD and the eigenvector corresponding to the smallest eigenvalue.

"""
function ispsd_wit(A::AbstractMatrix{T}; tol::Real=eps(real(T))^(3 / 4)) where T
    ishermitian(A) || return false, zeros(T,size(A, 1)) 
    S = (A + A') / 2
    evals, evecs = eigen(S)
    return evals[1] >= -tol, T(evecs[:, 1])
end