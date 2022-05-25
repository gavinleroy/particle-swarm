#include "gaussian_elimination_solver.h"
#include "helpers.h"
#include "perf_testers/perf_ge_solve.h"

#include <stdlib.h>

#define DIM 3
int main(int argc, char **argv)
{
  double A[DIM * DIM] = {1., 2., 3., 4., 5., 6, 7, 8, 10};
  double b[DIM] = {5., 6., 7.};

  size_t Ab_s = DIM * (DIM + 1);
  double *Ab = aligned_alloc(32, Ab_s * sizeof(double));
  double *x = malloc(DIM * sizeof(double));

  for (int i = 0; i < DIM; i++)
  {
    for (int j = 0; j < DIM; j++)
    {
      Ab[i * (DIM + 1) + j] = A[i * DIM + j];
    }
    Ab[i * (DIM + 1) + DIM] = b[i];
  }
  for (int i = 0; i < DIM; i++)
  {
    x[i] = 0;
  }

  print_rect_matrixd(Ab, DIM, DIM + 1, "Ab_before");
  print_vectord(x, DIM, "x_before");

  gaussian_elimination_solve(DIM, Ab, x);

  print_rect_matrixd(Ab, DIM, DIM + 1, "Ab_after");
  print_vectord(x, DIM, "x_after");

  for (int i = 0; i < DIM; i++)
  {
    for (int j = 0; j < DIM; j++)
    {
      Ab[i * (DIM + 1) + j] = A[i * DIM + j];
    }
    Ab[i * (DIM + 1) + DIM] = b[i];
  }

  for (int i = 0; i < DIM; i++)
  {
    x[i] = 0;
  }

  print_rect_matrixd(Ab, DIM, DIM + 1, "Ab_before");
  print_vectord(x, DIM, "x_before");

  // perf test
  register_functions_GE_SOLVE();
  perf_test_ge_solve(DIM, Ab, x);

  print_rect_matrixd(Ab, DIM, DIM + 1, "Ab_after");
  print_vectord(x, DIM, "x_after");

  free(x);
  free(Ab);
}
