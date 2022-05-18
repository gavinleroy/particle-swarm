module TestUtils

export column_major_to_row, array_to_column_major

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

end # module
