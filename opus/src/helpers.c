
#include "helpers.h"

#include "stdlib.h"
#include <inttypes.h>

#include "math.h"

void print_rect_matrixd(double const *A, int n, int m, char const *name)
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

void print_matrixd(double const *A, int N, char const *name)
{
  print_rect_matrixd(A, N, N, name);
}

void print_vectord(double const *M, int N, char const *name)
{
  printf("%s = [", name);
  int i;
  for (i = 0; i < N; i++)
  {
    printf("%.4f", M[i]);
    if (i < N - 1)
      printf(", ");
  }
  printf("]\n");
}

void print_vectori(int const *M, int N, char const *name)
{
  printf("%s = [", name);
  int i;
  for (i = 0; i < N; i++)
  {
    printf("%d", M[i]);
    if (i < N - 1)
      printf(", ");
  }
  printf("]\n");
}

void print_vectoru64(uint64_t const *M, int N, char const *name)
{
  printf("%s = [", name);
  int i;
  for (i = 0; i < N; i++)
  {
    printf("%" PRIu64 " ", M[i]);
    if (i < N - 1)
      printf(", ");
  }
  printf("]\n");
}

double rand_between(double a, double b)
{
  // see http://c-faq.com/lib/randrange.html
  return a + (double)rand() / ((double)RAND_MAX / (b - a));
}

double dist2(size_t dim, double const *x, double const *y)
{
  double s = 0;
  for (size_t i = 0; i < dim; i++)
  {
    double v = x[i] - y[i];
    s += v * v;
  }
  return s;
}

double dist(size_t dim, double const *x, double const *y)
{
  return sqrt(dist2(dim, x, y));
}
