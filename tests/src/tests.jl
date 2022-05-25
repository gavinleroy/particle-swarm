module tests

# To add a new FFI test stub, include it here ----
include("ffi_lu_solve.jl")
include("ffi_gaussian_elimination.jl")

# PSO_GE.solve_tests()
# PSO_GE.perf_tests()

LU.solve_tests()
# LU.perf_tests()

end # module
