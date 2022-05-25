extern "C"
{
#include "perf_lu_solve.h"
};

#include "cstring"

#include "../lu_solve.h"
#include "PerformanceTester.hpp"

template <> class ArgumentRestorer<lu_solve_fun_t>
{
private:
  int N0 = 0;
  double *A0 = nullptr;
  double *b0 = nullptr;

public:
  ArgumentRestorer<lu_solve_fun_t>(int N, double *A, int *ipiv, double *b)
  {
    N0 = N;

    size_t A_s = N * N;

    A0 = (double *)std::malloc(A_s * sizeof(double));
    b0 = (double *)std::malloc(N * sizeof(double));

    std::memcpy(A0, A, A_s * sizeof(double));
    std::memcpy(b0, b, N * sizeof(double));
  }

  void restore_arguments(int N, double *A, int *ipiv, double *b)
  {
    size_t A_s = N * N;
    std::memcpy(A, A0, A_s * sizeof(double));
    std::memcpy(b, b0, N * sizeof(double));
  }

  ~ArgumentRestorer<lu_solve_fun_t>()
  {
    free(b0);
    free(A0);
  }
};

static PerformanceTester<lu_solve_fun_t> perf_tester;

extern "C" void add_function_LU_SOLVE(lu_solve_fun_t f, char *name, int flop)
{
  std::string nm = std::string(name);
  std::cout << "Register " << nm << "\n";
  perf_tester.add_function(f, nm, flop);
}

extern "C" int perf_test_lu_solve(int N, double *A, int *ipiv, double *b)
{
  register_functions_LU_SOLVE();
  return perf_tester.perf_test_all_registered(N, A, ipiv, b);
}
