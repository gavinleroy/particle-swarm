
#include "helpers.h"

#include <stdio.h>
#include <inttypes.h>


void print_rect_matrixd(double const * A, int n, int m, char const *name)
{
  printf("%s = [\n", name);
  int i, j;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      printf("%.3f ", A[i * m + j]);
    }
    printf("\n");
  }
  printf("]\n");
}

void print_matrixd(double const * A, int N, char const *name)
{
  print_rect_matrixd(A, N, N, name);
}


void print_vectord(double const * M, int N, char const * name)
{
  printf("%s = [", name);
  int i;
  for (i = 0; i < N; i++)
  {
    printf("%.4f", M[i]);
    if (i < N - 1) printf(", ");
  }
  printf("]\n");
}

void print_vectori(int const * M, int N, char const * name)
{
  printf("%s = [", name);
  int i;
  for (i = 0; i < N; i++)
  {
    printf("%d", M[i]);
    if (i < N - 1) printf(", ");
  }
  printf("]\n");
}

void print_vectoru64(uint64_t const * M, int N, char const * name)
{
  printf("%s = [", name);
  int i;
  for (i = 0; i < N; i++)
  {
    printf("%" PRIu64 " ", M[i]);
    if (i < N - 1) printf(", ");
  }
  printf("]\n");
}