module PSO_GE

using Test: @test, @testset
using Printf: @sprintf
using TryCatch

# FFI stubs for LRU factorization
include("TestUtils.jl")

const tu = TestUtils

function ge_solve(N, Ab, x)
    # NOTE this preserve shouldn't be necessary because
    # a Ptr{Cdouble} Base.unsafe_convert already exists.
    GC.@preserve Ab x begin
        retcode = ccall(
            (:gaussian_elimination_solve, :libpso),
            Cint,                                  # return type
            (Csize_t, Ptr{Cdouble}, Ptr{Cdouble}), # parameter types
            N, Ab, x                               # actual arguments
        )
    end
    return retcode
end

function GE_solve(A::Matrix{<:Real}, b::Vector{<:Real})
    @assert ndims(A) == 2
    @assert size(A, 1) == size(A, 2)
    N = size(A, 1)
    Ab = hcat(A, b)
    Ab_vec = tu.alloc_aligned_vec(Cdouble, length(Ab))
    tu.fill_c_vec(Ab, Ab_vec)
    x = tu.alloc_aligned_vec(Cdouble, N)
    retcode = ge_solve(N, Ab_vec, x)
    @assert retcode == 0
    return x
end

# --------------------------------------

function solve_tests()
    @testset "GE Solve Tests" begin

        tu.starting_test("GE Solve")

        function test_ge_solve(A::Matrix, b::Vector)
            jx = A\b
            local x
            # NOTE if the system isnt' solvable
            # we will get an assertion error from GE_solve
            @try begin
                x = GE_solve(A, b)
            @catch e::AssertionError
                return
            end

            # @test A * jx ≈ b
            @test A * x ≈ b
            @test jx ≈ x
        end

        test_lambda = (n) -> begin
            M = tu.sym_n(n)
            b = rand(n)
            test_ge_solve(M, b)
        end

        @testset "Random GE Solve Small" begin
            tu.@run_random_N 100 2^5 2^7 "small ge solve" test_lambda
        end


        # @testset "Random GE Solve Large" begin
        #     tu.@run_random_N 1 2^10 2^12 "large ge solve" test_lambda
        # end

    end
end

function perf_tests()
    n = 512 # 2^8
    ab_n = n * (n + 1)
    local Ab
    Ab_vec = tu.alloc_aligned_vec(Cdouble, ab_n)
    x = tu.alloc_aligned_vec(Cdouble, n)

    while true
        A = tu.sym_n(n)
        b = rand(n)
        Ab = hcat(A, b)
        @assert length(Ab) == ab_n
        tu.fill_c_vec(Ab, Ab_vec)
        # Generate random matrices until we find one with a solution
        ge_solve(n, Ab_vec, x) != 0 || break
    end

    tu.fill_c_vec(Ab, Ab_vec)

    tu.starting_test(@sprintf "GE perf comparison with A[%d, %d]x = b[%d]" n n n)

    # NOTE this preserve shouldn't be necessary because
    # a Ptr{Cdouble} Base.unsafe_convert already exists.
    GC.@preserve Ab x begin
        retcode = ccall(
            (:perf_test_ge_solve, :libpso),
            Cint,                                  # return type
            (Csize_t, Ptr{Cdouble}, Ptr{Cdouble}), # parameter types
            n, Ab_vec, x                           # actual arguments
        )
    end
    @assert retcode == 0
end

end #module
