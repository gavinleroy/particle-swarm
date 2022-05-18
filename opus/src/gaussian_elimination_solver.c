
#include "gaussian_elimination_solver.h"

#include <math.h>
#include <stdio.h>

#include "helpers.h"

#define DEBUG_GE_SOLVER 0

static void swapd(double *a, int i, int j)
{
  double t = a[i];
  a[i] = a[j];
  a[j] = t;
}

// Textual transformation assumes matrix Ab and matrix size `N`
// are in scope.
#define MAT_Ab(ROW, COL) (Ab)[(N + 1) * (ROW) + (COL)]
int gaussian_elimination_solve(int N, double *Ab, double *x)
{
  // Note: works inplace on Ab, outputs in x
  // Ab= [ A | b], sie N x (N+1)

  // Gaussian elimination with partial pivoting

#if DEBUG_GE_SOLVER
  printf("Before elimination:");
  print_rect_matrixd(Ab, N, N + 1, "Ab0");
#endif

  for (int k = 0; k < N - 1; k++)
  {
    // Find largest possible pivot in submatrix of A
    double p = 0.;
    int pivot_row_idx = -1;
    for (int i = k; i < N; i++)
    {
      double v = MAT_Ab(i, k);
      if (fabs(v) > fabs(p))
      {
        p = v;
        pivot_row_idx = i;
      }
    }

    if (pivot_row_idx < 0)
    {
      // singular matrix
      fprintf(
          stderr,
          "ERROR: gaussian elimination failed: cannot find non-zero pivot for "
          "sub-matrix %d\n",
          k);
      return -1;
    }

    if (k != pivot_row_idx)
    {
// swap the rows in Ab
#if DEBUG_GE_SOLVER
      printf("Swap rows %d <-> %d\n", k, pivot_row_idx);
#endif
      for (int j = 0; j < N + 1; j++)
      {
        swapd(Ab, k * (N + 1) + j, pivot_row_idx * (N + 1) + j);
      }
    }

    // elimination: on A __and b__
    // note: the first substraction could be skipped / or replaced by =0
    // as it is known to give value 0 (and not used in back substitution)
    // Keep for debugging for now

    for (int i = k + 1; i < N; i++)
    {
      double r = MAT_Ab(i, k) / p;
      for (int j = k; j < N + 1; j++)
      {
        MAT_Ab(i, j) -= r * MAT_Ab(k, j);
      }
    }

#if DEBUG_GE_SOLVER
    printf("Elimination step %d:\n", k);
    print_rect_matrixd(Ab, N, N + 1, "Ab");
#endif
  }

  if (fabs(MAT_Ab(N - 1, N - 1)) < 1e-3)
  {
    // singular matrix
    fprintf(stderr, "ERROR: gaussian elimination failed: last pivot is 0\n");
    return -1;
  }

  // A is now upper triangular

  // Backward substitution
  // U x = y
  for (int i = N - 1; i >= 0; i--)
  {
    // value in b
    double v = MAT_Ab(i, N);
    for (int j = i + 1; j < N; j++)
    {
      v -= MAT_Ab(i, j) * x[j];
    }
    x[i] = v / MAT_Ab(i, i);
  }

  return 0;
}
