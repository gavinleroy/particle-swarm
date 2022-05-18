module PLU

using LinearAlgebra

# FFI stubs for LRU factorization
include("TestUtils.jl")

const tu = TestUtils

# PLU Factorization struct

struct PLU_factorization_C
    L::Ptr{Cdouble}
    U::Ptr{Cdouble}
    p::Ptr{Cint}
end

struct PLU_factorization
    L::Vector{Cdouble}
    U::Vector{Cdouble}
    p::Vector{Cint}
end

# See https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/#Garbage-Collection-Safety
# - Julia objects that are not rooted (i.e not alive) can get GC-ed.
# - T* types corresponds to Ref{T} for (non-array) inout arguments whose memory is managed by Julia, and to Ptr{T} for arrays
#   see https://github.com/JuliaLang/julia/issues/29850).
#   If the input object (of type T) does not have a C-like structure, __do not__ explicitly convert it to a type T'
#   instead define an overload for unsafe_convert(T', x::T) (if unsafe), cconvert (if safe but might change type),
#   or convert (if safe and type stable).
# - Never return a Ptr from a `convert` call ; ccall will call unsafe_convert if necessary (and make sure the
#   object is not GC-ed during the execution of the ccall). 
# - Alternatively: manually convert the object `obj` to pointer but mark the region surrounding the ccall as
#   GC.@preserve obj begin ... end (and more if the C code keeps a reference to the object!)
# QUESTION: relationship between convert/unsafe_convert/cconvert?

# unsafe_convert because it offers no guarentee that the memory in the arrays
# stays available.
function Base.convert(::Type{PLU_factorization_C}, x::PLU_factorization)
    PLU_factorization_C(pointer(x.L), pointer(x.U), pointer(x.p))
end

function alloc_PLU_factorization(N)::PLU_factorization
    L = Vector{Cdouble}(undef, N * N)
    U = Vector{Cdouble}(undef, N * N)
    p = Vector{Cint}(undef, N)
    return PLU_factorization(L, U, p)
end

function PLU_factorize(M::Matrix{Cdouble})

    @assert ndims(M) == 2
    @assert size(M, 1) == size(M, 2)

    N = size(M, 1) # M should be an NxN matrix
    Mp::Vector{Cdouble} = tu.column_major_to_row(M)
    plu::PLU_factorization = alloc_PLU_factorization(N)

    GC.@preserve plu begin
        retcode = ccall(
            (:plu_factorize, :libpso),
            Cint,                                           # return type
            (Cint, Ptr{Cdouble}, Ref{PLU_factorization_C}), # parameter types
            N, Mp, plu                                      # actual arguments
        )

        @assert retcode == 0

        lr = tu.array_to_column_major(plu.L, N, N)
        ur = tu.array_to_column_major(plu.U, N, N)
        # C code is 0-base indexed but Julia is 1-base
        pr = plu.p .+ 1
    end

    return (lr, ur, pr)
end

# function Vector{Float64}(m::Matrix{Float64})
#     tu.column_major_to_array(m)
# end

function PLU_solve(M::Matrix{Cdouble}, b::Vector{Cdouble})::Vector{Cdouble}
    N = size(M, 1)
    L, U, p = lu(M) # Use Julia's LA functions for isolation

    # Convert Matrices
    Mp::Vector{Cdouble} = tu.column_major_to_row(M)
    Lp::Vector{Cdouble} = tu.column_major_to_row(L)
    Up::Vector{Cdouble} = tu.column_major_to_row(U)
    x::Vector{Cdouble} = Vector{Cdouble}(undef, N)
    p = p .- 1 # Julia 1-based indexing at it again

    plu = PLU_factorization(Lp, Up, p)

    GC.@preserve plu begin
        retcode = ccall(
            (:plu_solve, :libpso),
            Cint,
            (Cint, Ref{PLU_factorization_C}, Ptr{Cdouble}, Ptr{Cdouble}),
            N, plu, b, x
        )
    end

    @assert retcode == 0

    return x
end

end # module
