module TestUtils

using Logging
using Printf: @sprintf, @info
using ProgressBars

const ldl = Base.Libc.Libdl

global dlopen_flags = (ldl.RTLD_LAZY | ldl.RTLD_DEEPBIND | ldl.RTLD_LOCAL)
global libpso

export alloc_aligned_vec,
    fill_c_vec,
    column_major_to_row,
    array_to_column_major,
    starting_test,
    sym_n,
    set_lib,
    lookup,
    @run_random

function set_lib(sym)
    sym = Symbol(sym)
    global libpso = ldl.dlopen(
        string(sym),
        dlopen_flags;
        throw_error=true
    )
    return
end

function lookup(library, sym::Symbol)
    return ldl.dlsym(library, Symbol(sym))
end

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
    @info "[starting]: $msg\n"
end

function sym_n(n)
    A = rand(n, n)
    return A+A'
end

# Matrices generated are of form:
# |  Q  | P |
# | P^t | 0 |
#  === AND ===
# | P |  Q  |
# | 0 | P^t |
# Where we assume that the dimensions are: 
# | c x c | c x d | * | d |  =  | d |
# | d x c | d x d |   | c |     | c |
# === AND ===
# | c x d | c x c | * | c |  =  | c |
# | d x d | d x c |   | d |     | d |
function build_matrices(c::Int64, d::Int64)
    P = rand(c, d)
    Q = sym_n(c)

    Z = zeros(d, d)

    R1 = hcat(Q, P)
    R2 = hcat(P', Z)

    R = vcat(R1, R2)

    M1 = hcat(P, Q)
    M2 = hcat(Z, P')

    M = vcat(M1, M2)

    R, M
end

# Matrix generated is of form:
# | P |  Q  |
# | 0 | P^t |
# Where we assume that the dimensions are: 
# | c x d | c x c | * | c |  =  | c |
# | d x d | d x c |   | d |     | d |
function build_block_triangular_matrix(c::Int64, d::Int64)
    P = rand(c, d)
    Q = sym_n(c)

    Z = zeros(d, d)

    M1 = hcat(P, Q)
    M2 = hcat(Z, P')

    M = vcat(M1, M2)

    # display(M)
    M
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
