
function partialtrace(a, syss::Vector{Int}, dims::Vector{Int})
    syss_sorted = sort(syss; rev=true)
    dims = copy(dims)
    for sys in syss_sorted
        a = partialtrace(a, sys, dims)
        dims = deleteat!(dims, sys)
    end
    return a
end
