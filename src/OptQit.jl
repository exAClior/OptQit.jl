module OptQit

using Convex,SCS, LinearAlgebra
# Write your package code here.

export ispsd, ispsd_wit
include("channel.jl")

include("state.jl")
end
