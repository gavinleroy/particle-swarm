module PSO_GE


# FFI stubs for LRU factorization
include("TestUtils.jl")

const tu = TestUtils

function GE_solve(A::Matrix{<:Real}, b::Vector{<:Real})
    @assert ndims(A) == 2
    @assert size(A, 1) == size(A, 2)

    N = size(A, 1)
    Ab = hcat(A, b)
    Ab = tu.column_major_to_row(Ab)
    x = Vector{Cdouble}(undef, N)

    GC.@preserve Ab begin
        retcode = ccall(
            (:gaussian_elimination_solve, :libpso),
            Cint,                                  # return type
            (Csize_t, Ptr{Cdouble}, Ptr{Cdouble}), # parameter types
            N, Ab, x                               # actual arguments
        )
    end

    return x
end

end #module
