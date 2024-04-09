module OptQit

using Convex, SCS, LinearAlgebra
# Write your package code here.

export ispsd, ispsd_wit
include("channel.jl")

include("state.jl")

import Convex.partialtrace
export partialtrace
include("utils.jl")
end
