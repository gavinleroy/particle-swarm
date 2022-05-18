module PSO_BLOOM

struct rounding_bloom end

function rounding_bloom_alloc()::Ptr{rounding_bloom}
    ccall(
        (:rounding_bloom_alloc, :libpso),
        Ptr{rounding_bloom},
        ()
      )
end


function rounding_bloom_init(bloom::Ptr{rounding_bloom}, entries::Int, error::Real, epsilon::Real, dim::Int, lower_bound::Vector{T}) where {T<:Real}
    ret = ccall(
        (:rounding_bloom_init, :libpso),
        Cint,
        (Ptr{rounding_bloom}, Cint, Cdouble, Cdouble, Cint, Ptr{Cdouble}),
        bloom, entries, error, epsilon, dim, lower_bound
      )
    ret
end


function rounding_bloom_check_add(bloom::Ptr{rounding_bloom}, x::Vector{T}, add::Bool) where {T<:Real}
    ret = ccall(
        (:rounding_bloom_check_add, :libpso),
        Cint,
        (Ptr{rounding_bloom}, Cint, Ptr{Cdouble}, Cint),
        bloom, length(x), x, add
      )
    ret
end


function rounding_bloom_print(bloom::Ptr{rounding_bloom})
    ret = ccall(
        (:rounding_bloom_print, :libpso),
        Cvoid,
        (Ptr{rounding_bloom},),
        bloom
      )
    ret
end

end #module
