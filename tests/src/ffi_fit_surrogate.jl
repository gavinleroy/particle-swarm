module surrogate

include("TestUtils.jl")

const tu = TestUtils

function lambda(A::Array{Cdouble})

end

# TODO change to @ctype ...
type blackbox_fn = Int32

end # module
