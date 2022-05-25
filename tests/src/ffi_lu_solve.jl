module LU

using Test: @test, @testset
using Printf: @sprintf
using LinearAlgebra

include("TestUtils.jl")

const tu = TestUtils

function lu_solve(N, A, p, b)
    GC.@preserve A b begin
        retcode = ccall(
            (:lu_solve, :libpso),
            Cint,
            (Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}),
            N, A, p, b
        )
    end
    return retcode
end

function LU_solve(M::Matrix{Cdouble}, b::Vector{Cdouble})::Vector{Cdouble}
    @assert ndims(M) == 2
    @assert size(M, 1) == size(M, 2)
    N = size(M, 1) # M should be an NxN matrix

    Mp::AbstractArray = tu.alloc_aligned_vec(Cdouble, N * N)
    ipiv::AbstractVector = tu.alloc_aligned_vec(Cint, N)
    bp::AbstractVector = tu.alloc_aligned_vec(Cdouble, N)

    tu.fill_c_vec(M, Mp)
    tu.fill_c_vec(b, bp)

    retcode = lu_solve(N, Mp, ipiv, bp)
    @assert retcode == 0

    return bp
end

# --------------------------------------

function solve_tests()
    @testset "LU Solve Tests" begin

        tu.starting_test("LU Solve")

        function test_lu_solve(A::Matrix, b::Vector)
            jx = A\b
            x = LU_solve(A, b)
            @test A * x ≈ b
            @test jx ≈ x
        end

        test_lambda = (n) -> begin
            M = rand(n, n)
            b = rand(n)
            test_lu_solve(M, b)
        end


        @testset "Random LU Solve Small" begin
            # tu.@run_random_N 100 2^5 2^8 "small lu solve" test_lambda
            tu.@run_random_N 2 50 600 "small lu solve" test_lambda
        end


        # @testset "Random LU Solve Large" begin
        #     tu.@run_random_N 1 2^10 2^12 "large lu solve" test_lambda
        # end

    end
end

function perf_tests()
    n = 2^10

    local A
    local b

    A_vec = tu.alloc_aligned_vec(Cdouble, n * n)
    b_vec = tu.alloc_aligned_vec(Cdouble, n)
    ipiv = tu.alloc_aligned_vec(Cint, n)

    while true
        A = tu.sym_n(n)
        b = rand(n)
        tu.fill_c_vec(A, A_vec)
        tu.fill_c_vec(b, b_vec)
        # Generate random matrices until we find one with a solution
        lu_solve(n, A_vec, ipiv, b_vec) != 0 || break
    end

    tu.fill_c_vec(A, A_vec)
    tu.fill_c_vec(b, b_vec)

    tu.starting_test(@sprintf "LU perf comparison with A[%d, %d]x = b[%d]" n n n)

    # NOTE this preserve shouldn't be necessary because
    # a Ptr{Cdouble} Base.unsafe_convert already exists.
    GC.@preserve A_vec b_vec begin
        retcode = ccall(
            (:perf_test_lu_solve, :libpso),
            Cint,                                             # return type
            (Csize_t, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}), # parameter types
            n, A_vec, ipiv, b_vec                             # actual arguments
        )
    end
    @assert retcode == 0
end


end # module
