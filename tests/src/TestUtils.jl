module TestUtils

using Printf: @printf, @sprintf
using ProgressBars

export alloc_aligned_vec,
    fill_c_vec,
    column_major_to_row,
    array_to_column_major,
    starting_test,
    sym_n,
    @run_random

# Changes the major order of the array memory layout.
# Row -> Column and vice-versa
# This functions works on *any* sized 2-dimensional array
# and not just square matrices.
change_major_order(X::AbstractArray, sizes...=size(X)...) =
    permutedims(reshape(X, sizes...), length(sizes):-1:1)

column_major_to_row(M::AbstractArray)::AbstractVector =
    collect(reshape(change_major_order(M), length(M)))

array_to_column_major(A::AbstractArray, sizes...)::AbstractArray =
    collect(change_major_order(reshape(A, sizes)))

function alloc_aligned_vec(::Type{T}, dims...) where T
    # @assert isbits(T)
    local n
    align = 32
    n = prod(dims)
    if n % align != 0
        n = ((n ÷ align) + 1) * align
    end
    sz = sizeof(T) * n
    @assert n % align == 0
    pt = ccall(:aligned_alloc, Ptr{Cvoid}, (Base.Csize_t, Base.Csize_t), align, sz)
    # pt = Libc.calloc(prod(dims), sizeof(T))
    @assert pt != C_NULL
    return unsafe_wrap(Array, convert(Ptr{T}, pt), dims; own = true)

end

function fill_c_vec(A::AbstractArray, v::AbstractVector)
    @assert ndims(A) == 2
    @assert length(A) == length(v)
    (rows, cols) = size(A)
    for i in 0:(rows - 1)
        for j in 0:(cols - 1)
            v[(i * cols + j) + 1] = A[i + 1, j + 1]
        end
    end
end

function fill_c_vec(V::AbstractVector, v::AbstractVector)
    @assert length(V) == length(v)
    N = length(V)
    for i in 1:N
        v[i] = V[i]
    end
end

function starting_test(msg)
    @printf "[starting]: %s\n" msg
end

function sym_n(n)
    A = rand(n, n)
    return A+A'
end

# XXX I know this is a little much. Especially the required thunk
# that will be invoked with the iterating N. But ...
# I don't personally like declaring variables that get shared accross
# compilation environments (I.E an unhygienic use of 'n')
macro run_random_N(iters, step, MAX_N, msg, thunk)
    es = eval(step)
    ei = eval(iters)
    MN = eval(MAX_N)
    all = ((MN / es) * ei)
    iters = Iterators.product((1:ei), (es:es:MN))
    total_tests = ProgressBar(iters)
    m = @sprintf "%d random %s instances" all msg
    return quote
        starting_test($m)
        @time for (_, n) in $total_tests
            $(esc(thunk))(n)
        end
    end
end

end # module
