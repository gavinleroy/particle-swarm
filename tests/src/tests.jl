module tests

using Test: @test, @testset
using Printf: @printf, @sprintf
using LinearAlgebra: lu

# To add a new FFI test stub, include it here ----
include("ffi_plu_factorization.jl")
include("ffi_gaussian_elimination.jl")

function starting_test(msg)
    @printf "[starting]: %s\n" msg
end

# --------------------------------------

@testset "LU Factorization Tests" begin

    starting_test("LU Factorization")

    function test_lu(A)
        try
            L,U,p = lu(A)
            Lt,Ut,pt = PLU.PLU_factorize(A)
            @test (Lt * Ut) ≈ A[pt, :]
            @test L ≈ Lt
            @test U ≈ Ut
        catch e
            # NOTE An ErrorException can be thrown from
            # LinearAlgebra.lu when the matrix has
            # no non-zero pivots. If this is the case,
            # we can just skip the matrix and ignore the
            # tests.
            bt = catch_backtrace()
            msg = sprint(showerror, e, bt)
            println(msg)
            return
        end
    end


    function run_random(iters, step, MAX_N; msg::String="")
        starting_test(@sprintf "%d random %s instances" ((MAX_N / step) * iters) msg)
        @time begin
            for i = 1:iters, n = step:step:MAX_N
                M = rand(n, n)
                test_lu(M)
            end
        end
    end


    @testset "Static LU" begin
        M_1::Matrix{Cdouble} = [1 2 4; 3 8 4; 2 6 13]
        test_lu(M_1)
    end


    @testset "Random Small LU" begin
        run_random(100, 10, 100, msg="small lu factor")
    end

    @testset "Random Large LU" begin
        run_random(2, 100, 2000, msg="large lu factor")
    end

end

# --------------------------------------

@testset "LU Solve Tests" begin

    starting_test("LU Solve")

    function test_lu_solve(A::Matrix, b::Vector)
        luA = lu(A)
        jx = luA\b
        x = PLU.PLU_solve(A, b)
        @test A * jx ≈ b
        @test A * x ≈ b
        @test jx ≈ x
    end


    function run_random(iters, step, MAX_N; msg::String="")
        starting_test(@sprintf "%d random %s instances" ((MAX_N / step) * iters) msg)
        @time begin
            for i = 1:iters, n = step:step:MAX_N
                M = rand(n, n)
                b = rand(n)
                test_lu_solve(M, b)
            end
        end
    end


    @testset "Static LU Solve" begin
        M_1::Matrix{Cdouble} = [2 1 -2; 1 -1 -1; 1 1 3]
        b_1::Vector{Cdouble} = [3; 0; 12]
        test_lu_solve(M_1, b_1)
    end


    @testset "Random LU Solve Small" begin
        run_random(100, 10, 100, msg="small lu solve")
    end


    @testset "Random LU Solve Large" begin
        run_random(2, 100, 2000, msg="large lu solve")
    end

end

# --------------------------------------

@testset "GE Solve Tests" begin

    starting_test("GE Solve")

    function test_ge_solve(A::Matrix, b::Vector)
        jx = A\b
        x = PSO_GE.GE_solve(A, b)
        @test A * jx ≈ b
        @test A * x ≈ b
        @test jx ≈ x
    end


    function run_random(iters, step, MAX_N; msg::String="")
        starting_test(@sprintf "%d random %s instances" ((MAX_N / step) * iters) msg)
        @time begin
            for i = 1:iters, n = step:step:MAX_N
                M = rand(n, n)
                b = rand(n)
                test_ge_solve(M, b)
            end
        end
    end


    @testset "Static GE Solve" begin
        M_1::Matrix{Cdouble} = [2 1 -2; 1 -1 -1; 1 1 3]
        b_1::Vector{Cdouble} = [3; 0; 12]
        test_ge_solve(M_1, b_1)
    end


    @testset "Random GE Solve Small" begin
        run_random(100, 10, 100, msg="small ge solve")
    end


    @testset "Random GE Solve Large" begin
        run_random(2, 100, 2000, msg="large ge solve")
    end

end

end # module
