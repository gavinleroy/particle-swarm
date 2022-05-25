#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "lu_solve.h"
// #include "perf_testers/perf_lu_solve.h"

// #define N 32768L
#define N 527L

static double drand()
{
  double a = 5.0;
  return ((double)rand() / (double)(RAND_MAX)) * a;
}

int main()
{

  printf("Solving\n");

  double *A = (double *)aligned_alloc(32, N * N * sizeof(double));
  double *b = (double *)aligned_alloc(32, N * sizeof(double));
  int *ipiv = (int *)aligned_alloc(32, N * sizeof(int));

  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      A[i * N + j] = drand();
  for (int i = 0; i < N; ++i)
    b[i] = drand();

  int ret = lu_solve(N, A, ipiv, b);

  printf("Done\n");
  // perf_test_lu_solve(N, A, ipiv, b);

  return 0;
}
