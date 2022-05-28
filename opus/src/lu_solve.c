#include "lu_solve.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "helpers.h"

#ifdef TEST_MKL
#include "mkl.h"
#endif

#ifdef TEST_PERF
#include "perf_testers/perf_lu_solve.h"
#include "perf_testers/perf_mmm.h"
#endif

// XXX Assume 'A' is an NxN defined in scope
#define AIX(ROW, COL) (A)[(N) * (ROW) + (COL)]
#define IX(ROW, COL) ((N) * (ROW) + (COL))
// NOTE When working with a matrix inset in another, you must
// index into it using this macro. Explicitly specifying the
// leading dimension which *must* be available.
#define MIX(MAT, LDIM, ROW, COL) (MAT)[((LDIM) * (ROW) + (COL))]

#define ONE 1.E0
#define ERR_THRESHOLD 1.0E-6 // FIXME is this small / big enough

#define APPROX_EQUAL(l, r) (fabs((l) - (r)) <= ERR_THRESHOLD)
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

int lu_solve_0(int N, double *A, int *ipiv, double *b);
int lu_solve_1(int N, double *A, int *ipiv, double *b);
int lu_solve_2(int N, double *A, int *ipiv, double *b);
#ifdef TEST_MKL
int lu_solve_3(int N, double *A, int *ipiv, double *b);
int lu_solve_4(int N, double *A, int *ipiv, double *b);
#endif

/** @brief Entry function to solve system A * x = b
 *         After exit b is overwritten with solution vector x.
 *
 * @param N Number of columns and rows of A.
 * @param A Real valued NxN matrix A.
 * @param ipiv Buffer for internal usage when pivoting.
 * @param b Real valued Nx1 vector b.
 */
int lu_solve(int N, double *A, int *ipiv, double *b)
{
  return lu_solve_2(N, A, ipiv, b);
}

// -----------------
// LAPACK routines
#if 0
// NOTE these functions have no declaration. However, a
// corresponding <name>_XXX  definition is made for each
// iteration XXX. This way performance testing is consistent.
// Forward declaration are here merely to have useful doc
// comments in one place.
// And YES, they are incorrectly using an (s) for single
// precision floats when they should be using a (d) for
// double precision (but I didn't notice until too late).

/** @brief Finds the index of the first element having maximum
 *         absolute value.
 *
 * @param N Number of elements in A.
 * @param A Real valued vector.
 * @param stride Space between elements
 * @return Index of the first element with the maximum absolute value.
 */
int isamax(int N, double *A, int stride);

/** @brief Swap interchanges two vectors.
 *
 * @param N The number of elements in vectors X and Y to swap.
 * @param X The first vector of real valued elements to swap.
 * @param incx The storage stride between elements in X.
 * @param Y The second vector of real valued elements to swap.
 * @param incy The storage stride between elements in Y.
 */
void sswap(int N, double *X, int incx, double *Y, int incy);

/** @brief Perform a series of row interchanges, one for each row
 *         k1 - k2 of A.
 *
 * @param N The number of column in MxN matrix A.
 * @param A MxN real valued matrix.
 * @param LDA Leading dimension of A.
 * @param k1 The first element of ipiv to interchange rows.
 * @param k2 (k2 - k1) elements of ipiv to do row interchanges.
 * @param ipiv The array of vector pivot indices such that i <-> ipiv[i].
 * @param incx The increment between indices in ipiv.
 */
void slaswp(int N, double *A, int LDA, int k1, int k2, int *ipiv, int incx);

/** @brief Solves matrix equation A * X = B
 *         A is assumed to be Non-transposed (L)ower
 *         triangular with a unit diagonal
 *         After exit, the matrix B is overwritten with
 *         the solution matrix X.
 *
 *  @param M Number of rows (height) of B.
 *  @param N Number of cols (width) of B.
 *  @param A Real valued matrix A.
 *  @param LDA Leading dimension of A.
 *  @param B Real valued matrix B.
 *  @param LDB Leading dimension of B.
 */
void strsm_L(int M, int N, double *A, int LDA, double *B, int LDB);

/** @brief equivalent to strsm_L except A is an (U)pper triangular matrix
 *         with a *Non-unit* diagonal. It is still assumed to be
non-transposed.
 */
void strsm_U(int M, int N, double *A, int LDA, double *B, int LDB);

// Assume No-transpose for both matrices
/** @brief compute C := alpha * A * B + beta * C
 *
 * @param M Number of rows (height) of A and C.
 * @param N Number of cols (width) of B and C.
 * @param K Number of cols (width) of A and rows (height) of B.
 * @param alpha Scalar alpha.
 * @param A Real valued MxK matrix A.
 * @param LDA Leading dimension of A.
 * @param B Real valued KxN matrix B.
 * @param LDB Leading dimension of B.
 * @param beta Scalar beta.
 * @param C Real valued MxN matrix C.
 * @param Leading dimension of C.
 */
void sgemm(int M, int N, int K, double alpha, double *A, int LDA, double *B,
            int LDB, double beta, double *C, int LDC);

/** @brief Factor A = P * L * U in place using BLAS1 / BLAS2 functions
 *
 * @param M The number of rows (height) of A.
 * @param N The number of columns (width) of A.
 * @param A Real valued matrix in which to factor [L\U].
 * @param LDA The leading dimension of the matrix in memory A.
 * @param ipiv Pivot indices for A.
 */
int sgetf2(int M, int N, double *A, int LDA, int *ipiv);

// Equivalent to SGETRS (no-transpose).
/** @brief Solves system of linear equations A * x = b.
 *         After exit b is overwritten with solution vector x.
 *
 * @param N Number of rows and columns in A.
 * @param A Real valued NxN matrix A which has been factored into [L\U].
 * @param ipiv Pivot indices used when factoring A.
 * @param b Real valued vector with N elements.
 */
int sgetrs(int N, double *A, int *ipiv, double *b);

#endif // FALSE to remove the above documentation.

// -----------------
// Utilities

static void swapd(double *a, int i, int j)
{
  double t = a[i];
  a[i] = a[j];
  a[j] = t;
}

static int ideal_block(int M, int N) { return 64; }

// ------------------------------------------------------------------
// Implementation start

// NOTE All functions preceded by a _ are internal and unexposed by the API.
// They are prefixed because they often have the same name with a different
// signature than those exposed.

/** ------------------------------------------------------------------
 * Base implementation
 */
int lu_solve_0(int N, double *A, int *ipiv, double *b)
{
  int i, j, k;

  for (i = 0; i < N; ++i)
  {

    // == Partial Pivoting ==

    double p_t;
    double p_v;
    int p_i;

    p_v = fabs(AIX(i, i));
    p_i = i;

    for (k = i + 1; k < N; ++k)
    {
      p_t = fabs(AIX(k, i));
      if (p_t > p_v)
      {
        p_v = p_t;
        p_i = k;
      }
    }

    if (APPROX_EQUAL(p_v, 0.))
    {
      fprintf(stderr, "ERROR: LU Solve singular matrix\n");
      return -1;
    }

    if (i != p_i)
    {

#if DEBUG_LU_SOLVER
      printf("p %f pivot row %d\n", p_v, p_i);
      printf("swap rows %d <-> %d\n", i, p_i);
#endif

      // Swap immediately the b
      swapd(b, i, p_i);

      // Swap rows k and p_i
      for (j = 0; j < N; ++j)
      {
        swapd(A, IX(i, j), IX(p_i, j));
      }
    }

    // === ===

    // BLAS 1 Scale vector
    for (j = i + 1; j < N; ++j)
    {
      AIX(j, i) = AIX(j, i) / AIX(i, i);
    }

    // BLAS 2 rank-1 update
    for (j = i + 1; j < N; ++j)
    {
      for (int k = i + 1; k < N; ++k)
      {
        AIX(j, k) = AIX(j, k) - AIX(j, i) * AIX(i, k);
      }
    }
  }

  // A now contains L (below diagonal)
  //                U (above diagonal)

  // Forward substitution
  for (k = 0; k < N; ++k)
  {
    if (!APPROX_EQUAL(b[k], 0.))
    {
      for (i = k + 1; i < N; ++i)
      {
        b[i] = b[i] - b[k] * AIX(i, k);
      }
    }
  }

  // Backward substitution
  for (k = N - 1; k >= 0; --k)
  {
    if (!APPROX_EQUAL(b[k], 0.))
    {
      b[k] = b[k] / AIX(k, k); // Non-unit diagonal
      for (i = 0; i < k; ++i)
      {
        b[i] = b[i] - b[k] * AIX(i, k);
      }
    }
  }

  return 0;
}

/** ------------------------------------------------------------------
 * Base implementation and LAPACK base impls.
 *
 * Using blocked (outer function) and delayed updates.
 */
static int isamax_1(int N, double *A, int stride)
{
  // NOTE N can be 0 if we iterate to the end of the rows in _lu_solve_XXX
  // this saves us from initializing the ipiv array. We could also do this
  // by just setting the end index to itself.

  // assert(0 < N);
  if (0 == N)
    return 0;

  assert(0 < stride);

  // TODO speacialize the case when `STRIDE` == 1.

  double p_v, p_t;
  int i, p_i, ix;

  p_v = fabs(A[0]);
  p_i = 0;
  ix = 0 + stride;

  for (i = 1; i < N; ++i, ix += stride)
  {
    assert(ix == i * stride);

    p_t = fabs(A[ix]);
    if (p_t > p_v)
    {
      p_v = p_t;
      p_i = i;
    }
  }
  return p_i;
}

static void sswap_1(int N, double *X, int incx, double *Y, int incy)
{
  // TODO special case when incx == incy == 1
  // ^^^ This is our case actually (row-major iteration).

  assert(0 < N);
  assert(0 < incx);
  assert(0 < incy);

  double t;
  int i, ix, iy;

  for (i = 0, ix = 0, iy = 0; i < N; ++i, ix += incx, iy += incy)
  {
    t = X[ix];
    X[ix] = Y[iy];
    Y[iy] = t;
  }

  return;
}

static void slaswp_1(int N, double *A, int LDA, int k1, int k2, int *ipiv,
                     int incx)
{
  // interchange row i with row IPIV[k1 + (i - k1) * abs(incx)]

  assert(0 < incx); // Do not cover the negative case

  int ix, ix0;
  int i, k, p_i;
  double tmp;
  ix0 = k1;
  ix = ix0;

  for (i = k1; i < k2; ++i, ix += incx)
  {
    p_i = ipiv[ix];
    if (p_i != i)
    {
      for (k = 0; k < N; ++k)
      {
        tmp = MIX(A, LDA, i, k);
        MIX(A, LDA, i, k) = MIX(A, LDA, p_i, k);
        MIX(A, LDA, p_i, k) = tmp;
      }
    }
  }
}

static void strsm_L_1(int M, int N, double *A, int LDA, double *B, int LDB)
{
  int i, j, k;

  for (j = 0; j < N; ++j)
  {
    for (k = 0; k < M; ++k)
    {
      if (!APPROX_EQUAL(MIX(B, LDB, k, j), 0.))
      {
        for (i = k + 1; i < M; ++i)
        {
          MIX(B, LDB, i, j) =
              MIX(B, LDB, i, j) - MIX(B, LDB, k, j) * MIX(A, LDA, i, k);
        }
      }
    }
  }
}

static void strsm_U_1(int M, int N, double *A, int LDA, double *B, int LDB)
{
  int i, j, k;

  for (j = 0; j < N; ++j)
  {
    for (k = M - 1; k >= 0; --k)
    {
      if (!APPROX_EQUAL(MIX(B, LDB, k, j), 0.))
      {
        MIX(B, LDB, k, j) = MIX(B, LDB, k, j) / MIX(A, LDA, k, k);
        for (i = 0; i < k; ++i)
        {
          MIX(B, LDB, i, j) =
              MIX(B, LDB, i, j) - MIX(B, LDB, k, j) * MIX(A, LDA, i, k);
        }
      }
    }
  }
}

static void sgemm_1(int M, int N, int K, double alpha, double *A, int LDA,
                    double *B, int LDB, double beta, double *C, int LDC)
{

  // NOTE as written below, we specialize to alpha = -1 beta = 1
  assert(APPROX_EQUAL(beta, ONE));
  assert(APPROX_EQUAL(alpha, -ONE));

  int i, j, k;
  double tmp;

  for (j = 0; j < N; ++j)
  {
    // BETA = 1
    for (k = 0; k < K; ++k)
    {
      tmp = alpha * MIX(B, LDB, k, j);
      for (i = 0; i < M; ++i)
      {
        MIX(C, LDC, i, j) = MIX(C, LDC, i, j) + tmp * MIX(A, LDA, i, k);
      }
    }
  }
}

void sgemm_intel(int M, int N, int K, double alpha, double *A, int LDA,
                 double *B, int LDB, double beta, double *C, int LDC)
{
  // Update trailing submatrix
  return cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha,
                     A, LDA, B, LDB, beta, C, LDC);
}

int sgetf2_1(int M, int N, double *A, int LDA, int *ipiv)
{

  int i, j;

  // NOTE I have increased this to iterate *until* N however you could
  // stop at < N - 1 if you cover the bounds case.
  for (i = 0; i < MIN(M, N); ++i)
  {
    // --------
    // ISAMAX
    double p_v;
    int p_i;

    p_i = i + isamax_1(M - i, &MIX(A, LDA, i, i), LDA);
    p_v = MIX(A, LDA, p_i, i);

    if (APPROX_EQUAL(p_v, 0.))
    {
      fprintf(stderr, "ERROR: LU Solve singular matrix\n");
      return -1;
    }

    ipiv[i] = p_i;

    if (i != p_i)
    {
      sswap_1(N, &MIX(A, LDA, i, 0), 1, &MIX(A, LDA, p_i, 0), 1);
    }

    // TODO compute the minimum machine safe integer such that
    // 1 / fmin doesn't overflow. This can be inlind because
    // we are targeting SKYLAKE.

    // BLAS 1 Scale vector
    for (j = i + 1; j < M; ++j)
    {
      MIX(A, LDA, j, i) = MIX(A, LDA, j, i) / MIX(A, LDA, i, i);
    }

    if (i < MIN(M, N))
    {

      // BLAS 2 rank-1 update
      for (j = i + 1; j < M; ++j)
      {
        for (int k = i + 1; k < N; ++k)
        {
          MIX(A, LDA, j, k) =
              MIX(A, LDA, j, k) - MIX(A, LDA, j, i) * MIX(A, LDA, i, k);
        }
      }
    }
  }

  return 0;
}

// Equivalent to SGETRS (no-transpose).
static int sgetrs_1(int N, double *A, int *ipiv, double *b)
{
  // A now contains L (below diagonal)
  //                U (above diagonal)
  // Swap pivot rows in b
  slaswp_1(1, b, 1, 0, N, ipiv, 1);

#if DEBUG_LU_SOLVER
  printf("B: ");
  for (int i = 0; i < N; ++i)
    printf("%.4f  ", b[i]);
  printf("\n");
#endif

  // Forward substitution
  strsm_L_1(N, 1, A, N, b, 1);
  // Backward substitution
  strsm_U_1(N, 1, A, N, b, 1);
  return 0;
}

int lu_solve_1(int N, double *A, int *ipiv, double *b)
{
  /**
   * FIXME for simplicity we will choose a block size of 32.
   * In the real system we can choose a dynamic block (better), or we can
   * always choose a fixed block and handle clean-up cases afterward (worse).
   */

  int NB = 32, retcode;

  int M = N;   // Square matrix
  int LDA = N; // N is the leading dimension

  int ib, IB, k;

  // BLocked factor A into [L \ U]

  if (N < NB)
  {
    // Use unblocked code
    retcode = sgetf2_1(M, N, A, LDA, ipiv);

#if DEBUG_LU_SOLVER
    printf("IPIV: ");
    for (i = 0; i < N; ++i)
      printf("%d  ", ipiv[i]);
    printf("\n");
#endif

    if (retcode != 0)
      return retcode;
  }
  else
  {

    for (ib = 0; ib < MIN(M, N); ib += NB)
    {
      IB = MIN(MIN(M, N) - ib, NB);

      retcode = sgetf2_1(M - ib, IB, &AIX(ib, ib), LDA, ipiv + ib);

      if (retcode != 0)
        return retcode;

      // Update the pivot indices
      for (k = ib; k < MIN(M, ib + IB); ++k)
      {
        ipiv[k] = ipiv[k] + ib;
      }

      // Apply interchanges to columns 0 : ib
      slaswp_1(ib, A, LDA, ib, ib + IB, ipiv, 1);

      if (ib + IB < N)
      {
        // Apply interchanges to columns ib + IB : N
        slaswp_1(N - ib - IB, &AIX(0, ib + IB), LDA, ib, ib + IB, ipiv, 1);

        // Compute the block row of U
        strsm_L_1(IB, N - ib - IB, &AIX(ib, ib), LDA, &AIX(ib, ib + IB), LDA);

        if (ib + IB < M) // NOTE for a square matrix this will always be true.
        {
          // Update trailing submatrix
          sgemm_1(M - ib - IB, N - ib - IB, IB, -ONE, //
                  &AIX(ib + IB, ib), LDA,             //
                  &AIX(ib, ib + IB), LDA,             //
                  ONE,                                //
                  &AIX(ib + IB, ib + IB), LDA         //
          );
        }
      }
    }
  }

  // Solve the system with A
  retcode = sgetrs_1(N, A, ipiv, b);

  return retcode;
}

/** ------------------------------------------------------------------
 * Base implementation and LAPACK base impls.
 * Using C optimizations, but NO vectorization (only setting it up).
 *
 * Using blocked (outer function) and delayed updates.
 *
 */
// FIXME XXX we can play with NB
static int isamax_2(int N, double *A, int stride)
{
  // NOTE N can be 0 if we iterate to the end of the rows in _lu_solve_XXX
  // this saves us from initializing the ipiv array. We could also do this
  // by just setting the end index to itself.

  assert(0 <= N);

  // TODO FIXME
  if (N < 1 || stride == 0)
    return 0;

  assert(0 < stride);

  // TODO speacialize the case when `STRIDE` == 1.

  double p_v, p_t;
  int i, p_i, ix;

  p_v = fabs(A[0]);
  p_i = 0;
  ix = 0 + stride;

  for (i = 1; i < N; ++i, ix += stride)
  {
    assert(ix == i * stride);

    p_t = fabs(A[ix]);
    if (p_t > p_v)
    {
      p_v = p_t;
      p_i = i;
    }
  }
  return p_i;
}

static void sswap_2(int N, double *X, int incx, double *Y, int incy)
{
  // TODO special case when incx == incy == 1
  // ^^^ This is our case actually (row-major iteration).

  assert(0 < N);
  assert(0 < incx);
  assert(0 < incy);

  double   //
      t,   //
      t_0, //
      t_1, //
      t_2, //
      t_3, //
      t_4, //
      t_5, //
      t_6, //
      t_7, //

      x_i_0, //
      x_i_1, //
      x_i_2, //
      x_i_3, //
      x_i_4, //
      x_i_5, //
      x_i_6, //
      x_i_7, //

      y_i_0, //
      y_i_1, //
      y_i_2, //
      y_i_3, //
      y_i_4, //
      y_i_5, //
      y_i_6, //
      y_i_7  //
      ;

  int i, ix, iy;

  if (incx == 1 && incy == 1)
  {
    for (i = 0; i < N - 7; i += 8)
    {

      x_i_0 = X[i + 0];
      x_i_1 = X[i + 1];
      x_i_2 = X[i + 2];
      x_i_3 = X[i + 3];
      x_i_4 = X[i + 4];
      x_i_5 = X[i + 5];
      x_i_6 = X[i + 6];
      x_i_7 = X[i + 7];

      y_i_0 = Y[i + 0];
      y_i_1 = Y[i + 1];
      y_i_2 = Y[i + 2];
      y_i_3 = Y[i + 3];
      y_i_4 = Y[i + 4];
      y_i_5 = Y[i + 5];
      y_i_6 = Y[i + 6];
      y_i_7 = Y[i + 7];

      X[i + 0] = y_i_0;
      X[i + 1] = y_i_1;
      X[i + 2] = y_i_2;
      X[i + 3] = y_i_3;
      X[i + 4] = y_i_4;
      X[i + 5] = y_i_5;
      X[i + 6] = y_i_6;
      X[i + 7] = y_i_7;

      Y[i + 0] = x_i_0;
      Y[i + 1] = x_i_1;
      Y[i + 2] = x_i_2;
      Y[i + 3] = x_i_3;
      Y[i + 4] = x_i_4;
      Y[i + 5] = x_i_5;
      Y[i + 6] = x_i_6;
      Y[i + 7] = x_i_7;
    }

    for (; i < N; ++i)
    {
      t_0 = X[i + 0];
      X[i + 0] = Y[i + 0];
      Y[i + 0] = t_0;
    }
  }
  else
  {
    for (i = 0, ix = 0, iy = 0; i < N; ++i, ix += incx, iy += incy)
    {
      t = X[ix];
      X[ix] = Y[iy];
      Y[iy] = t;
    }
  }
}

static void slaswp_2(int N, double *A, int LDA, int k1, int k2, int *ipiv,
                     int incx)
{
  // NOTE ipiv is layed out sequentially in memory so we can special case it
  assert(0 < incx);  // Do not cover the negative case
  assert(1 == incx); // Do not cover the negative case

  int n32;
  int i, j, k, p_i;

  double tmp;
  double        //
      a__i_k_0, //
      a__i_k_1, //
      a__i_k_2, //
      a__i_k_3, //
      a__i_k_4, //
      a__i_k_5, //
      a__i_k_6, //
      a__i_k_7, //

      a_pi_k_0, //
      a_pi_k_1, //
      a_pi_k_2, //
      a_pi_k_3, //
      a_pi_k_4, //
      a_pi_k_5, //
      a_pi_k_6, //
      a_pi_k_7  //
      ;

  n32 = (N / 32) * 32;

  if (0 < n32)
  {
    for (j = 0; j < n32; j += 32)
    {
      for (i = k1; i < k2; ++i)
      {
        p_i = ipiv[i];
        if (p_i != i)
        {
          for (k = j; k < j + 32 - 7; k += 8)
          {
            a__i_k_0 = MIX(A, LDA, i, k + 0);
            a__i_k_1 = MIX(A, LDA, i, k + 1);
            a__i_k_2 = MIX(A, LDA, i, k + 2);
            a__i_k_3 = MIX(A, LDA, i, k + 3);
            a__i_k_4 = MIX(A, LDA, i, k + 4);
            a__i_k_5 = MIX(A, LDA, i, k + 5);
            a__i_k_6 = MIX(A, LDA, i, k + 6);
            a__i_k_7 = MIX(A, LDA, i, k + 7);

            a_pi_k_0 = MIX(A, LDA, p_i, k + 0);
            a_pi_k_1 = MIX(A, LDA, p_i, k + 1);
            a_pi_k_2 = MIX(A, LDA, p_i, k + 2);
            a_pi_k_3 = MIX(A, LDA, p_i, k + 3);
            a_pi_k_4 = MIX(A, LDA, p_i, k + 4);
            a_pi_k_5 = MIX(A, LDA, p_i, k + 5);
            a_pi_k_6 = MIX(A, LDA, p_i, k + 6);
            a_pi_k_7 = MIX(A, LDA, p_i, k + 7);

            MIX(A, LDA, i, k + 0) = a_pi_k_0;
            MIX(A, LDA, i, k + 1) = a_pi_k_1;
            MIX(A, LDA, i, k + 2) = a_pi_k_2;
            MIX(A, LDA, i, k + 3) = a_pi_k_3;
            MIX(A, LDA, i, k + 4) = a_pi_k_4;
            MIX(A, LDA, i, k + 5) = a_pi_k_5;
            MIX(A, LDA, i, k + 6) = a_pi_k_6;
            MIX(A, LDA, i, k + 7) = a_pi_k_7;

            MIX(A, LDA, p_i, k + 0) = a__i_k_0;
            MIX(A, LDA, p_i, k + 1) = a__i_k_1;
            MIX(A, LDA, p_i, k + 2) = a__i_k_2;
            MIX(A, LDA, p_i, k + 3) = a__i_k_3;
            MIX(A, LDA, p_i, k + 4) = a__i_k_4;
            MIX(A, LDA, p_i, k + 5) = a__i_k_5;
            MIX(A, LDA, p_i, k + 6) = a__i_k_6;
            MIX(A, LDA, p_i, k + 7) = a__i_k_7;
          }
        }
      }
    }
  }

  // Leftover cases from blocking by 32
  if (n32 != N)
  {
    for (i = k1; i < k2; ++i)
    {
      p_i = ipiv[i];
      if (i != p_i)
      {
        for (k = n32; k < N - 7; k += 8)
        {
          a__i_k_0 = MIX(A, LDA, i, k + 0);
          a__i_k_1 = MIX(A, LDA, i, k + 1);
          a__i_k_2 = MIX(A, LDA, i, k + 2);
          a__i_k_3 = MIX(A, LDA, i, k + 3);
          a__i_k_4 = MIX(A, LDA, i, k + 4);
          a__i_k_5 = MIX(A, LDA, i, k + 5);
          a__i_k_6 = MIX(A, LDA, i, k + 6);
          a__i_k_7 = MIX(A, LDA, i, k + 7);

          a_pi_k_0 = MIX(A, LDA, p_i, k + 0);
          a_pi_k_1 = MIX(A, LDA, p_i, k + 1);
          a_pi_k_2 = MIX(A, LDA, p_i, k + 2);
          a_pi_k_3 = MIX(A, LDA, p_i, k + 3);
          a_pi_k_4 = MIX(A, LDA, p_i, k + 4);
          a_pi_k_5 = MIX(A, LDA, p_i, k + 5);
          a_pi_k_6 = MIX(A, LDA, p_i, k + 6);
          a_pi_k_7 = MIX(A, LDA, p_i, k + 7);

          MIX(A, LDA, i, k + 0) = a_pi_k_0;
          MIX(A, LDA, i, k + 1) = a_pi_k_1;
          MIX(A, LDA, i, k + 2) = a_pi_k_2;
          MIX(A, LDA, i, k + 3) = a_pi_k_3;
          MIX(A, LDA, i, k + 4) = a_pi_k_4;
          MIX(A, LDA, i, k + 5) = a_pi_k_5;
          MIX(A, LDA, i, k + 6) = a_pi_k_6;
          MIX(A, LDA, i, k + 7) = a_pi_k_7;

          MIX(A, LDA, p_i, k + 0) = a__i_k_0;
          MIX(A, LDA, p_i, k + 1) = a__i_k_1;
          MIX(A, LDA, p_i, k + 2) = a__i_k_2;
          MIX(A, LDA, p_i, k + 3) = a__i_k_3;
          MIX(A, LDA, p_i, k + 4) = a__i_k_4;
          MIX(A, LDA, p_i, k + 5) = a__i_k_5;
          MIX(A, LDA, p_i, k + 6) = a__i_k_6;
          MIX(A, LDA, p_i, k + 7) = a__i_k_7;
        }

        for (; k < N; ++k)
        {
          a__i_k_0 = MIX(A, LDA, i, k + 0);
          a_pi_k_0 = MIX(A, LDA, p_i, k + 0);
          MIX(A, LDA, p_i, k + 0) = a__i_k_0;
          MIX(A, LDA, i, k + 0) = a_pi_k_0;
        }
      }
    }
  }
}

static void strsm_L_2(int M, int N, double *A, int LDA, double *B, int LDB)
{
  int i, j, k;

  for (j = 0; j < N; ++j)
  {
    for (k = 0; k < M; ++k)
    {
      for (i = k + 1; i < M; ++i)
      {
        MIX(B, LDB, i, j) =
            MIX(B, LDB, i, j) - MIX(B, LDB, k, j) * MIX(A, LDA, i, k);
      }
    }
  }
}

static void strsm_U_2(int M, int N, double *A, int LDA, double *B, int LDB)
{
  int i, j, k;

  for (j = 0; j < N; ++j)
  {
    for (k = M - 1; k >= 0; --k)
    {
      MIX(B, LDB, k, j) = MIX(B, LDB, k, j) / MIX(A, LDA, k, k);
      for (i = 0; i < k; ++i)
      {
        MIX(B, LDB, i, j) =
            MIX(B, LDB, i, j) - MIX(B, LDB, k, j) * MIX(A, LDA, i, k);
      }
    }
  }
}

static void sgemm_2(int M, int N, int K, double alpha, double *restrict A,
                    int LDA, double *restrict B, int LDB, double beta,
                    double *restrict C, int LDC)
{

  /**
   * A[ M x K ] * B[ K x N ] = C[ M x N ]
   */

  // NOTE as written below, we specialize to alpha = -1 beta = 1
  assert(APPROX_EQUAL(beta, ONE));
  assert(APPROX_EQUAL(alpha, -ONE));

  // ---------------
  // https://en.wikichip.org/wiki/intel/microarchitectures/skylake_(client)
  //
  // NOTE Skylake has a 32000 B L1 cache
  // With cache lines of 64 B

#define CACHE_BLOCK 56
  const int NB = CACHE_BLOCK;
  const int DBL_LMT = 512 * 2;

  double __attribute__((aligned(32))) AL[CACHE_BLOCK * CACHE_BLOCK];
  double __attribute__((aligned(32))) BL[CACHE_BLOCK * CACHE_BLOCK];
  double __attribute__((aligned(32))) CL[CACHE_BLOCK * CACHE_BLOCK];

  // NOTE choose MU + NU + MU * NU <= 16
  const int MU = 4;
  const int NU = 2;
  const int KU = 8;

  int i, j, k,      //
      ii, jj, kk,   //
      iii, jjj, kkk //
      ;

  i = 0;
  ii = 0;
  iii = 0;

  j = 0;
  jj = 0;
  jjj = 0;

  k = 0;
  kk = 0;
  kkk = 0;

  double             //
      A_iii_0_kkk_0, //
      A_iii_0_kkk_1, //
      A_iii_0_kkk_2, //
      A_iii_0_kkk_3, //
      A_iii_0_kkk_4, //
      A_iii_0_kkk_5, //
      A_iii_0_kkk_6, //
      A_iii_0_kkk_7, //

      A_iii_1_kkk_0, //
      A_iii_1_kkk_1, //
      A_iii_1_kkk_2, //
      A_iii_1_kkk_3, //
      A_iii_1_kkk_4, //
      A_iii_1_kkk_5, //
      A_iii_1_kkk_6, //
      A_iii_1_kkk_7, //

      A_iii_2_kkk_0, //
      A_iii_2_kkk_1, //
      A_iii_2_kkk_2, //
      A_iii_2_kkk_3, //
      A_iii_2_kkk_4, //
      A_iii_2_kkk_5, //
      A_iii_2_kkk_6, //
      A_iii_2_kkk_7, //

      A_iii_3_kkk_0, //
      A_iii_3_kkk_1, //
      A_iii_3_kkk_2, //
      A_iii_3_kkk_3, //
      A_iii_3_kkk_4, //
      A_iii_3_kkk_5, //
      A_iii_3_kkk_6, //
      A_iii_3_kkk_7, //
      A_iii_3_kkk_8, //

      C_iii_0_jjj_0, //
      C_iii_0_jjj_1, //
      C_iii_1_jjj_0, //
      C_iii_1_jjj_1, //
      C_iii_2_jjj_0, //
      C_iii_2_jjj_1, //
      C_iii_3_jjj_0, //
      C_iii_3_jjj_1, //

      B_kkk_0_jjj_0, //
      B_kkk_0_jjj_1, //
      B_kkk_1_jjj_0, //
      B_kkk_1_jjj_1, //
      B_kkk_2_jjj_0, //
      B_kkk_2_jjj_1, //
      B_kkk_3_jjj_0, //
      B_kkk_3_jjj_1, //
      B_kkk_4_jjj_0, //
      B_kkk_4_jjj_1, //
      B_kkk_5_jjj_0, //
      B_kkk_5_jjj_1, //
      B_kkk_6_jjj_0, //
      B_kkk_6_jjj_1, //
      B_kkk_7_jjj_0, //
      B_kkk_7_jjj_1  //

      ;

  // -----------------------------------------------------------------
  // Basic unrolled MMM with no register blocking

  // NOTE this should be our case as the incoming matrix A is tall and skinny.
  if (1 || M > N) // FIXME remove automatic true
  {

    // i blocked
    for (i = 0; i <= M - NB; i += NB)
    {
      // j blocked
      for (j = 0; j <= N - NB; j += NB)
      {

        // k blocked
        for (k = 0; k <= K - NB; k += NB)

          // Cache blocking
          for (ii = i; ii <= (i + NB) - MU; ii += MU)
            for (jj = j; jj <= (j + NB) - NU; jj += NU)
              for (kk = k; kk <= (k + NB) - KU; kk += KU)

              /* // Register blocking */
              /* for (kkk = kk; kkk < kk + KU; ++kkk) */
              /*   for (iii = ii; iii < ii + MU; ++iii) */
              /*     for (jjj = jj; jjj < jj + NU; ++jjj) */
              {

                // FIXME reorder stores and loads

                C_iii_0_jjj_0 = MIX(C, LDC, ii + 0, jj + 0);
                C_iii_0_jjj_1 = MIX(C, LDC, ii + 0, jj + 1);
                C_iii_1_jjj_0 = MIX(C, LDC, ii + 1, jj + 0);
                C_iii_1_jjj_1 = MIX(C, LDC, ii + 1, jj + 1);
                C_iii_2_jjj_0 = MIX(C, LDC, ii + 2, jj + 0);
                C_iii_2_jjj_1 = MIX(C, LDC, ii + 2, jj + 1);
                C_iii_3_jjj_0 = MIX(C, LDC, ii + 3, jj + 0);
                C_iii_3_jjj_1 = MIX(C, LDC, ii + 3, jj + 1);

                B_kkk_0_jjj_0 = MIX(B, LDB, kk + 0, jj + 0);
                B_kkk_1_jjj_0 = MIX(B, LDB, kk + 1, jj + 0);
                B_kkk_2_jjj_0 = MIX(B, LDB, kk + 2, jj + 0);
                B_kkk_3_jjj_0 = MIX(B, LDB, kk + 3, jj + 0);
                B_kkk_4_jjj_0 = MIX(B, LDB, kk + 4, jj + 0);
                B_kkk_5_jjj_0 = MIX(B, LDB, kk + 5, jj + 0);
                B_kkk_6_jjj_0 = MIX(B, LDB, kk + 6, jj + 0);
                B_kkk_7_jjj_0 = MIX(B, LDB, kk + 7, jj + 0);

                B_kkk_0_jjj_1 = MIX(B, LDB, kk + 0, jj + 1);
                B_kkk_1_jjj_1 = MIX(B, LDB, kk + 1, jj + 1);
                B_kkk_2_jjj_1 = MIX(B, LDB, kk + 2, jj + 1);
                B_kkk_3_jjj_1 = MIX(B, LDB, kk + 3, jj + 1);
                B_kkk_4_jjj_1 = MIX(B, LDB, kk + 4, jj + 1);
                B_kkk_5_jjj_1 = MIX(B, LDB, kk + 5, jj + 1);
                B_kkk_6_jjj_1 = MIX(B, LDB, kk + 6, jj + 1);
                B_kkk_7_jjj_1 = MIX(B, LDB, kk + 7, jj + 1);

                A_iii_0_kkk_0 = MIX(A, LDA, ii + 0, kk + 0);
                A_iii_0_kkk_1 = MIX(A, LDA, ii + 0, kk + 1);
                A_iii_0_kkk_2 = MIX(A, LDA, ii + 0, kk + 2);
                A_iii_0_kkk_3 = MIX(A, LDA, ii + 0, kk + 3);
                A_iii_0_kkk_4 = MIX(A, LDA, ii + 0, kk + 4);
                A_iii_0_kkk_5 = MIX(A, LDA, ii + 0, kk + 5);
                A_iii_0_kkk_6 = MIX(A, LDA, ii + 0, kk + 6);
                A_iii_0_kkk_7 = MIX(A, LDA, ii + 0, kk + 7);

                A_iii_1_kkk_0 = MIX(A, LDA, ii + 1, kk + 0);
                A_iii_1_kkk_1 = MIX(A, LDA, ii + 1, kk + 1);
                A_iii_1_kkk_2 = MIX(A, LDA, ii + 1, kk + 2);
                A_iii_1_kkk_3 = MIX(A, LDA, ii + 1, kk + 3);
                A_iii_1_kkk_4 = MIX(A, LDA, ii + 1, kk + 4);
                A_iii_1_kkk_5 = MIX(A, LDA, ii + 1, kk + 5);
                A_iii_1_kkk_6 = MIX(A, LDA, ii + 1, kk + 6);
                A_iii_1_kkk_7 = MIX(A, LDA, ii + 1, kk + 7);

                A_iii_2_kkk_0 = MIX(A, LDA, ii + 2, kk + 0);
                A_iii_2_kkk_1 = MIX(A, LDA, ii + 2, kk + 1);
                A_iii_2_kkk_2 = MIX(A, LDA, ii + 2, kk + 2);
                A_iii_2_kkk_3 = MIX(A, LDA, ii + 2, kk + 3);
                A_iii_2_kkk_4 = MIX(A, LDA, ii + 2, kk + 4);
                A_iii_2_kkk_5 = MIX(A, LDA, ii + 2, kk + 5);
                A_iii_2_kkk_6 = MIX(A, LDA, ii + 2, kk + 6);
                A_iii_2_kkk_7 = MIX(A, LDA, ii + 2, kk + 7);

                A_iii_3_kkk_0 = MIX(A, LDA, ii + 3, kk + 0);
                A_iii_3_kkk_1 = MIX(A, LDA, ii + 3, kk + 1);
                A_iii_3_kkk_2 = MIX(A, LDA, ii + 3, kk + 2);
                A_iii_3_kkk_3 = MIX(A, LDA, ii + 3, kk + 3);
                A_iii_3_kkk_4 = MIX(A, LDA, ii + 3, kk + 4);
                A_iii_3_kkk_5 = MIX(A, LDA, ii + 3, kk + 5);
                A_iii_3_kkk_6 = MIX(A, LDA, ii + 3, kk + 6);
                A_iii_3_kkk_7 = MIX(A, LDA, ii + 3, kk + 7);
                A_iii_3_kkk_8 = MIX(A, LDA, ii + 3, kk + 8);

                // ------

                C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
                C_iii_0_jjj_1 -= B_kkk_0_jjj_1 * A_iii_0_kkk_0;
                C_iii_1_jjj_0 -= B_kkk_0_jjj_0 * A_iii_1_kkk_0;
                C_iii_1_jjj_1 -= B_kkk_0_jjj_1 * A_iii_1_kkk_0;
                C_iii_2_jjj_0 -= B_kkk_0_jjj_0 * A_iii_2_kkk_0;
                C_iii_2_jjj_1 -= B_kkk_0_jjj_1 * A_iii_2_kkk_0;
                C_iii_3_jjj_0 -= B_kkk_0_jjj_0 * A_iii_3_kkk_0;
                C_iii_3_jjj_1 -= B_kkk_0_jjj_1 * A_iii_3_kkk_0;

                C_iii_0_jjj_0 -= B_kkk_1_jjj_0 * A_iii_0_kkk_1;
                C_iii_0_jjj_1 -= B_kkk_1_jjj_1 * A_iii_0_kkk_1;
                C_iii_1_jjj_0 -= B_kkk_1_jjj_0 * A_iii_1_kkk_1;
                C_iii_1_jjj_1 -= B_kkk_1_jjj_1 * A_iii_1_kkk_1;
                C_iii_2_jjj_0 -= B_kkk_1_jjj_0 * A_iii_2_kkk_1;
                C_iii_2_jjj_1 -= B_kkk_1_jjj_1 * A_iii_2_kkk_1;
                C_iii_3_jjj_0 -= B_kkk_1_jjj_0 * A_iii_3_kkk_1;
                C_iii_3_jjj_1 -= B_kkk_1_jjj_1 * A_iii_3_kkk_1;

                C_iii_0_jjj_0 -= B_kkk_2_jjj_0 * A_iii_0_kkk_2;
                C_iii_0_jjj_1 -= B_kkk_2_jjj_1 * A_iii_0_kkk_2;
                C_iii_1_jjj_0 -= B_kkk_2_jjj_0 * A_iii_1_kkk_2;
                C_iii_1_jjj_1 -= B_kkk_2_jjj_1 * A_iii_1_kkk_2;
                C_iii_2_jjj_0 -= B_kkk_2_jjj_0 * A_iii_2_kkk_2;
                C_iii_2_jjj_1 -= B_kkk_2_jjj_1 * A_iii_2_kkk_2;
                C_iii_3_jjj_0 -= B_kkk_2_jjj_0 * A_iii_3_kkk_2;
                C_iii_3_jjj_1 -= B_kkk_2_jjj_1 * A_iii_3_kkk_2;

                C_iii_0_jjj_0 -= B_kkk_3_jjj_0 * A_iii_0_kkk_3;
                C_iii_0_jjj_1 -= B_kkk_3_jjj_1 * A_iii_0_kkk_3;
                C_iii_1_jjj_0 -= B_kkk_3_jjj_0 * A_iii_1_kkk_3;
                C_iii_1_jjj_1 -= B_kkk_3_jjj_1 * A_iii_1_kkk_3;
                C_iii_2_jjj_0 -= B_kkk_3_jjj_0 * A_iii_2_kkk_3;
                C_iii_2_jjj_1 -= B_kkk_3_jjj_1 * A_iii_2_kkk_3;
                C_iii_3_jjj_0 -= B_kkk_3_jjj_0 * A_iii_3_kkk_3;
                C_iii_3_jjj_1 -= B_kkk_3_jjj_1 * A_iii_3_kkk_3;

                C_iii_0_jjj_0 -= B_kkk_4_jjj_0 * A_iii_0_kkk_4;
                C_iii_0_jjj_1 -= B_kkk_4_jjj_1 * A_iii_0_kkk_4;
                C_iii_1_jjj_0 -= B_kkk_4_jjj_0 * A_iii_1_kkk_4;
                C_iii_1_jjj_1 -= B_kkk_4_jjj_1 * A_iii_1_kkk_4;
                C_iii_2_jjj_0 -= B_kkk_4_jjj_0 * A_iii_2_kkk_4;
                C_iii_2_jjj_1 -= B_kkk_4_jjj_1 * A_iii_2_kkk_4;
                C_iii_3_jjj_0 -= B_kkk_4_jjj_0 * A_iii_3_kkk_4;
                C_iii_3_jjj_1 -= B_kkk_4_jjj_1 * A_iii_3_kkk_4;

                C_iii_0_jjj_0 -= B_kkk_5_jjj_0 * A_iii_0_kkk_5;
                C_iii_0_jjj_1 -= B_kkk_5_jjj_1 * A_iii_0_kkk_5;
                C_iii_1_jjj_0 -= B_kkk_5_jjj_0 * A_iii_1_kkk_5;
                C_iii_1_jjj_1 -= B_kkk_5_jjj_1 * A_iii_1_kkk_5;
                C_iii_2_jjj_0 -= B_kkk_5_jjj_0 * A_iii_2_kkk_5;
                C_iii_2_jjj_1 -= B_kkk_5_jjj_1 * A_iii_2_kkk_5;
                C_iii_3_jjj_0 -= B_kkk_5_jjj_0 * A_iii_3_kkk_5;
                C_iii_3_jjj_1 -= B_kkk_5_jjj_1 * A_iii_3_kkk_5;

                C_iii_0_jjj_0 -= B_kkk_6_jjj_0 * A_iii_0_kkk_6;
                C_iii_0_jjj_1 -= B_kkk_6_jjj_1 * A_iii_0_kkk_6;
                C_iii_1_jjj_0 -= B_kkk_6_jjj_0 * A_iii_1_kkk_6;
                C_iii_1_jjj_1 -= B_kkk_6_jjj_1 * A_iii_1_kkk_6;
                C_iii_2_jjj_0 -= B_kkk_6_jjj_0 * A_iii_2_kkk_6;
                C_iii_2_jjj_1 -= B_kkk_6_jjj_1 * A_iii_2_kkk_6;
                C_iii_3_jjj_0 -= B_kkk_6_jjj_0 * A_iii_3_kkk_6;
                C_iii_3_jjj_1 -= B_kkk_6_jjj_1 * A_iii_3_kkk_6;

                C_iii_0_jjj_1 -= B_kkk_7_jjj_1 * A_iii_0_kkk_7;
                C_iii_0_jjj_0 -= B_kkk_7_jjj_0 * A_iii_0_kkk_7;
                C_iii_1_jjj_1 -= B_kkk_7_jjj_1 * A_iii_1_kkk_7;
                C_iii_1_jjj_0 -= B_kkk_7_jjj_0 * A_iii_1_kkk_7;
                C_iii_2_jjj_0 -= B_kkk_7_jjj_0 * A_iii_2_kkk_7;
                C_iii_2_jjj_1 -= B_kkk_7_jjj_1 * A_iii_2_kkk_7;
                C_iii_3_jjj_0 -= B_kkk_7_jjj_0 * A_iii_3_kkk_7;
                C_iii_3_jjj_1 -= B_kkk_7_jjj_1 * A_iii_3_kkk_7;

                // -----

                MIX(C, LDC, ii + 0, jj + 0) = C_iii_0_jjj_0;
                MIX(C, LDC, ii + 0, jj + 1) = C_iii_0_jjj_1;
                MIX(C, LDC, ii + 1, jj + 0) = C_iii_1_jjj_0;
                MIX(C, LDC, ii + 1, jj + 1) = C_iii_1_jjj_1;
                MIX(C, LDC, ii + 2, jj + 0) = C_iii_2_jjj_0;
                MIX(C, LDC, ii + 2, jj + 1) = C_iii_2_jjj_1;
                MIX(C, LDC, ii + 3, jj + 0) = C_iii_3_jjj_0;
                MIX(C, LDC, ii + 3, jj + 1) = C_iii_3_jjj_1;
              }

        // k overflow
        for (; k < K; ++k)

          // Cache blocking
          for (ii = i; ii <= (i + NB) - MU; ii += MU)
            for (jj = j; jj <= (j + NB) - NU; jj += NU)
            /* // Registier blocking */
            /* for (iii = ii; iii < ii + MU; ++iii) */
            /*   for (jjj = jj; jjj < jj + NU; ++jjj) */
            {
              C_iii_0_jjj_0 = MIX(C, LDC, ii + 0, jj + 0);
              C_iii_0_jjj_1 = MIX(C, LDC, ii + 0, jj + 1);

              C_iii_1_jjj_0 = MIX(C, LDC, ii + 1, jj + 0);
              C_iii_1_jjj_1 = MIX(C, LDC, ii + 1, jj + 1);

              C_iii_2_jjj_0 = MIX(C, LDC, ii + 2, jj + 0);
              C_iii_2_jjj_1 = MIX(C, LDC, ii + 2, jj + 1);

              C_iii_3_jjj_0 = MIX(C, LDC, ii + 3, jj + 0);
              C_iii_3_jjj_1 = MIX(C, LDC, ii + 3, jj + 1);

              B_kkk_0_jjj_0 = MIX(B, LDB, k + 0, jj + 0);
              B_kkk_0_jjj_1 = MIX(B, LDB, k + 0, jj + 1);

              A_iii_0_kkk_0 = MIX(A, LDA, ii + 0, k + 0);
              A_iii_1_kkk_0 = MIX(A, LDA, ii + 1, k + 0);
              A_iii_2_kkk_0 = MIX(A, LDA, ii + 2, k + 0);
              A_iii_3_kkk_0 = MIX(A, LDA, ii + 3, k + 0);

              // ----

              C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
              C_iii_0_jjj_1 -= B_kkk_0_jjj_1 * A_iii_0_kkk_0;
              C_iii_1_jjj_0 -= B_kkk_0_jjj_0 * A_iii_1_kkk_0;
              C_iii_1_jjj_1 -= B_kkk_0_jjj_1 * A_iii_1_kkk_0;
              C_iii_2_jjj_0 -= B_kkk_0_jjj_0 * A_iii_2_kkk_0;
              C_iii_2_jjj_1 -= B_kkk_0_jjj_1 * A_iii_2_kkk_0;
              C_iii_3_jjj_0 -= B_kkk_0_jjj_0 * A_iii_3_kkk_0;
              C_iii_3_jjj_1 -= B_kkk_0_jjj_1 * A_iii_3_kkk_0;

              // ----

              MIX(C, LDC, ii + 0, jj + 0) = C_iii_0_jjj_0;
              MIX(C, LDC, ii + 0, jj + 1) = C_iii_0_jjj_1;
              MIX(C, LDC, ii + 1, jj + 0) = C_iii_1_jjj_0;
              MIX(C, LDC, ii + 1, jj + 1) = C_iii_1_jjj_1;
              MIX(C, LDC, ii + 2, jj + 0) = C_iii_2_jjj_0;
              MIX(C, LDC, ii + 2, jj + 1) = C_iii_2_jjj_1;
              MIX(C, LDC, ii + 3, jj + 0) = C_iii_3_jjj_0;
              MIX(C, LDC, ii + 3, jj + 1) = C_iii_3_jjj_1;
            }
      }

      // j overflow
      for (; j < N; ++j)
      {

        // k blocked
        for (k = 0; k <= K - NB; k += NB)

          // Cache blocking
          for (ii = i; ii <= (i + NB) - MU; ii += MU)
            for (kk = k; kk <= (k + NB) - KU; kk += KU)

            // Registier blocking
            /* for (kkk = kk; kkk < kk + KU; ++kkk) */
            /*   for (iii = ii; iii < ii + MU; ++iii) */
            {

              C_iii_0_jjj_0 = MIX(C, LDC, ii + 0, j + 0);
              C_iii_1_jjj_0 = MIX(C, LDC, ii + 1, j + 0);
              C_iii_2_jjj_0 = MIX(C, LDC, ii + 2, j + 0);
              C_iii_3_jjj_0 = MIX(C, LDC, ii + 3, j + 0);

              B_kkk_0_jjj_0 = MIX(B, LDB, kk + 0, j + 0);
              B_kkk_1_jjj_0 = MIX(B, LDB, kk + 1, j + 0);
              B_kkk_2_jjj_0 = MIX(B, LDB, kk + 2, j + 0);
              B_kkk_3_jjj_0 = MIX(B, LDB, kk + 3, j + 0);
              B_kkk_4_jjj_0 = MIX(B, LDB, kk + 4, j + 0);
              B_kkk_5_jjj_0 = MIX(B, LDB, kk + 5, j + 0);
              B_kkk_6_jjj_0 = MIX(B, LDB, kk + 6, j + 0);
              B_kkk_7_jjj_0 = MIX(B, LDB, kk + 7, j + 0);

              A_iii_0_kkk_0 = MIX(A, LDA, ii + 0, kk + 0);
              A_iii_0_kkk_1 = MIX(A, LDA, ii + 0, kk + 1);
              A_iii_0_kkk_2 = MIX(A, LDA, ii + 0, kk + 2);
              A_iii_0_kkk_3 = MIX(A, LDA, ii + 0, kk + 3);
              A_iii_0_kkk_4 = MIX(A, LDA, ii + 0, kk + 4);
              A_iii_0_kkk_5 = MIX(A, LDA, ii + 0, kk + 5);
              A_iii_0_kkk_6 = MIX(A, LDA, ii + 0, kk + 6);
              A_iii_0_kkk_7 = MIX(A, LDA, ii + 0, kk + 7);

              A_iii_1_kkk_0 = MIX(A, LDA, ii + 1, kk + 0);
              A_iii_1_kkk_1 = MIX(A, LDA, ii + 1, kk + 1);
              A_iii_1_kkk_2 = MIX(A, LDA, ii + 1, kk + 2);
              A_iii_1_kkk_3 = MIX(A, LDA, ii + 1, kk + 3);
              A_iii_1_kkk_4 = MIX(A, LDA, ii + 1, kk + 4);
              A_iii_1_kkk_5 = MIX(A, LDA, ii + 1, kk + 5);
              A_iii_1_kkk_6 = MIX(A, LDA, ii + 1, kk + 6);
              A_iii_1_kkk_7 = MIX(A, LDA, ii + 1, kk + 7);

              A_iii_2_kkk_0 = MIX(A, LDA, ii + 2, kk + 0);
              A_iii_2_kkk_1 = MIX(A, LDA, ii + 2, kk + 1);
              A_iii_2_kkk_2 = MIX(A, LDA, ii + 2, kk + 2);
              A_iii_2_kkk_3 = MIX(A, LDA, ii + 2, kk + 3);
              A_iii_2_kkk_4 = MIX(A, LDA, ii + 2, kk + 4);
              A_iii_2_kkk_5 = MIX(A, LDA, ii + 2, kk + 5);
              A_iii_2_kkk_6 = MIX(A, LDA, ii + 2, kk + 6);
              A_iii_2_kkk_7 = MIX(A, LDA, ii + 2, kk + 7);

              A_iii_3_kkk_0 = MIX(A, LDA, ii + 3, kk + 0);
              A_iii_3_kkk_1 = MIX(A, LDA, ii + 3, kk + 1);
              A_iii_3_kkk_2 = MIX(A, LDA, ii + 3, kk + 2);
              A_iii_3_kkk_3 = MIX(A, LDA, ii + 3, kk + 3);
              A_iii_3_kkk_4 = MIX(A, LDA, ii + 3, kk + 4);
              A_iii_3_kkk_5 = MIX(A, LDA, ii + 3, kk + 5);
              A_iii_3_kkk_6 = MIX(A, LDA, ii + 3, kk + 6);
              A_iii_3_kkk_7 = MIX(A, LDA, ii + 3, kk + 7);
              A_iii_3_kkk_8 = MIX(A, LDA, ii + 3, kk + 8);

              // ------

              C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
              C_iii_1_jjj_0 -= B_kkk_0_jjj_0 * A_iii_1_kkk_0;
              C_iii_2_jjj_0 -= B_kkk_0_jjj_0 * A_iii_2_kkk_0;
              C_iii_3_jjj_0 -= B_kkk_0_jjj_0 * A_iii_3_kkk_0;

              C_iii_0_jjj_0 -= B_kkk_1_jjj_0 * A_iii_0_kkk_1;
              C_iii_1_jjj_0 -= B_kkk_1_jjj_0 * A_iii_1_kkk_1;
              C_iii_2_jjj_0 -= B_kkk_1_jjj_0 * A_iii_2_kkk_1;
              C_iii_3_jjj_0 -= B_kkk_1_jjj_0 * A_iii_3_kkk_1;

              C_iii_0_jjj_0 -= B_kkk_2_jjj_0 * A_iii_0_kkk_2;
              C_iii_1_jjj_0 -= B_kkk_2_jjj_0 * A_iii_1_kkk_2;
              C_iii_2_jjj_0 -= B_kkk_2_jjj_0 * A_iii_2_kkk_2;
              C_iii_3_jjj_0 -= B_kkk_2_jjj_0 * A_iii_3_kkk_2;

              C_iii_0_jjj_0 -= B_kkk_3_jjj_0 * A_iii_0_kkk_3;
              C_iii_1_jjj_0 -= B_kkk_3_jjj_0 * A_iii_1_kkk_3;
              C_iii_2_jjj_0 -= B_kkk_3_jjj_0 * A_iii_2_kkk_3;
              C_iii_3_jjj_0 -= B_kkk_3_jjj_0 * A_iii_3_kkk_3;

              C_iii_0_jjj_0 -= B_kkk_4_jjj_0 * A_iii_0_kkk_4;
              C_iii_1_jjj_0 -= B_kkk_4_jjj_0 * A_iii_1_kkk_4;
              C_iii_2_jjj_0 -= B_kkk_4_jjj_0 * A_iii_2_kkk_4;
              C_iii_3_jjj_0 -= B_kkk_4_jjj_0 * A_iii_3_kkk_4;

              C_iii_0_jjj_0 -= B_kkk_5_jjj_0 * A_iii_0_kkk_5;
              C_iii_1_jjj_0 -= B_kkk_5_jjj_0 * A_iii_1_kkk_5;
              C_iii_2_jjj_0 -= B_kkk_5_jjj_0 * A_iii_2_kkk_5;
              C_iii_3_jjj_0 -= B_kkk_5_jjj_0 * A_iii_3_kkk_5;

              C_iii_0_jjj_0 -= B_kkk_6_jjj_0 * A_iii_0_kkk_6;
              C_iii_1_jjj_0 -= B_kkk_6_jjj_0 * A_iii_1_kkk_6;
              C_iii_2_jjj_0 -= B_kkk_6_jjj_0 * A_iii_2_kkk_6;
              C_iii_3_jjj_0 -= B_kkk_6_jjj_0 * A_iii_3_kkk_6;

              C_iii_0_jjj_0 -= B_kkk_7_jjj_0 * A_iii_0_kkk_7;
              C_iii_1_jjj_0 -= B_kkk_7_jjj_0 * A_iii_1_kkk_7;
              C_iii_2_jjj_0 -= B_kkk_7_jjj_0 * A_iii_2_kkk_7;
              C_iii_3_jjj_0 -= B_kkk_7_jjj_0 * A_iii_3_kkk_7;

              // -----

              MIX(C, LDC, ii + 0, j + 0) = C_iii_0_jjj_0;
              MIX(C, LDC, ii + 1, j + 0) = C_iii_1_jjj_0;
              MIX(C, LDC, ii + 2, j + 0) = C_iii_2_jjj_0;
              MIX(C, LDC, ii + 3, j + 0) = C_iii_3_jjj_0;
            }

        // k overflow
        for (; k < K; ++k)

          // Cache blocking
          for (ii = i; ii <= (i + NB) - MU; ii += MU)

          // Register blocking
          /* for (iii = ii; iii < ii + MU; ++iii) */
          {
            C_iii_0_jjj_0 = MIX(C, LDC, ii + 0, j + 0);
            C_iii_1_jjj_0 = MIX(C, LDC, ii + 1, j + 0);
            C_iii_2_jjj_0 = MIX(C, LDC, ii + 2, j + 0);
            C_iii_3_jjj_0 = MIX(C, LDC, ii + 3, j + 0);

            B_kkk_0_jjj_0 = MIX(B, LDB, k + 0, j + 0);

            A_iii_0_kkk_0 = MIX(A, LDA, ii + 0, k + 0);
            A_iii_1_kkk_0 = MIX(A, LDA, ii + 1, k + 0);
            A_iii_2_kkk_0 = MIX(A, LDA, ii + 2, k + 0);
            A_iii_3_kkk_0 = MIX(A, LDA, ii + 3, k + 0);

            // ------

            C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
            C_iii_1_jjj_0 -= B_kkk_0_jjj_0 * A_iii_1_kkk_0;
            C_iii_2_jjj_0 -= B_kkk_0_jjj_0 * A_iii_2_kkk_0;
            C_iii_3_jjj_0 -= B_kkk_0_jjj_0 * A_iii_3_kkk_0;

            // -----

            MIX(C, LDC, ii + 0, j + 0) = C_iii_0_jjj_0;
            MIX(C, LDC, ii + 1, j + 0) = C_iii_1_jjj_0;
            MIX(C, LDC, ii + 2, j + 0) = C_iii_2_jjj_0;
            MIX(C, LDC, ii + 3, j + 0) = C_iii_3_jjj_0;
          }
      }
    }

    // i overflow
    for (; i < M; ++i)
    {
      // j blocked
      for (j = 0; j <= N - NB; j += NB)
      {

        // k blocked
        for (k = 0; k <= K - NB; k += NB)

          // Cache blocking
          for (jj = j; jj <= (j + NB) - NU; jj += NU)
            for (kk = k; kk <= (k + NB) - KU; kk += KU)

            // Register blocking
            /* for (kkk = kk; kkk < kk + KU; ++kkk) */
            /*   for (jjj = jj; jjj < jj + NU; ++jjj) */
            {

              C_iii_0_jjj_0 = MIX(C, LDC, i + 0, jj + 0);
              C_iii_0_jjj_1 = MIX(C, LDC, i + 0, jj + 1);

              B_kkk_0_jjj_0 = MIX(B, LDB, kk + 0, jj + 0);
              B_kkk_0_jjj_1 = MIX(B, LDB, kk + 0, jj + 1);
              B_kkk_1_jjj_0 = MIX(B, LDB, kk + 1, jj + 0);
              B_kkk_1_jjj_1 = MIX(B, LDB, kk + 1, jj + 1);
              B_kkk_2_jjj_0 = MIX(B, LDB, kk + 2, jj + 0);
              B_kkk_2_jjj_1 = MIX(B, LDB, kk + 2, jj + 1);
              B_kkk_3_jjj_0 = MIX(B, LDB, kk + 3, jj + 0);
              B_kkk_3_jjj_1 = MIX(B, LDB, kk + 3, jj + 1);
              B_kkk_4_jjj_0 = MIX(B, LDB, kk + 4, jj + 0);
              B_kkk_4_jjj_1 = MIX(B, LDB, kk + 4, jj + 1);
              B_kkk_5_jjj_0 = MIX(B, LDB, kk + 5, jj + 0);
              B_kkk_5_jjj_1 = MIX(B, LDB, kk + 5, jj + 1);
              B_kkk_6_jjj_0 = MIX(B, LDB, kk + 6, jj + 0);
              B_kkk_6_jjj_1 = MIX(B, LDB, kk + 6, jj + 1);
              B_kkk_7_jjj_0 = MIX(B, LDB, kk + 7, jj + 0);
              B_kkk_7_jjj_1 = MIX(B, LDB, kk + 7, jj + 1);

              A_iii_0_kkk_0 = MIX(A, LDA, i + 0, kk + 0);
              A_iii_0_kkk_1 = MIX(A, LDA, i + 0, kk + 1);
              A_iii_0_kkk_2 = MIX(A, LDA, i + 0, kk + 2);
              A_iii_0_kkk_3 = MIX(A, LDA, i + 0, kk + 3);
              A_iii_0_kkk_4 = MIX(A, LDA, i + 0, kk + 4);
              A_iii_0_kkk_5 = MIX(A, LDA, i + 0, kk + 5);
              A_iii_0_kkk_6 = MIX(A, LDA, i + 0, kk + 6);
              A_iii_0_kkk_7 = MIX(A, LDA, i + 0, kk + 7);

              // ------

              C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
              C_iii_0_jjj_1 -= B_kkk_0_jjj_1 * A_iii_0_kkk_0;

              C_iii_0_jjj_0 -= B_kkk_1_jjj_0 * A_iii_0_kkk_1;
              C_iii_0_jjj_1 -= B_kkk_1_jjj_1 * A_iii_0_kkk_1;

              C_iii_0_jjj_0 -= B_kkk_2_jjj_0 * A_iii_0_kkk_2;
              C_iii_0_jjj_1 -= B_kkk_2_jjj_1 * A_iii_0_kkk_2;

              C_iii_0_jjj_0 -= B_kkk_3_jjj_0 * A_iii_0_kkk_3;
              C_iii_0_jjj_1 -= B_kkk_3_jjj_1 * A_iii_0_kkk_3;

              C_iii_0_jjj_0 -= B_kkk_4_jjj_0 * A_iii_0_kkk_4;
              C_iii_0_jjj_1 -= B_kkk_4_jjj_1 * A_iii_0_kkk_4;

              C_iii_0_jjj_0 -= B_kkk_5_jjj_0 * A_iii_0_kkk_5;
              C_iii_0_jjj_1 -= B_kkk_5_jjj_1 * A_iii_0_kkk_5;

              C_iii_0_jjj_0 -= B_kkk_6_jjj_0 * A_iii_0_kkk_6;
              C_iii_0_jjj_1 -= B_kkk_6_jjj_1 * A_iii_0_kkk_6;

              C_iii_0_jjj_0 -= B_kkk_7_jjj_0 * A_iii_0_kkk_7;
              C_iii_0_jjj_1 -= B_kkk_7_jjj_1 * A_iii_0_kkk_7;

              // -----

              MIX(C, LDC, i + 0, jj + 0) = C_iii_0_jjj_0;
              MIX(C, LDC, i + 0, jj + 1) = C_iii_0_jjj_1;
            }

        // k overflow
        for (; k < K; ++k)

          // Cache blocking
          for (jj = j; jj <= (j + NB) - NU; jj += NU)

          // Register blocking
          /* for (jjj = jj; jjj < jj + NU; ++jjj) */
          {

            C_iii_0_jjj_0 = MIX(C, LDC, i + 0, jj + 0);
            C_iii_0_jjj_1 = MIX(C, LDC, i + 0, jj + 1);

            B_kkk_0_jjj_0 = MIX(B, LDB, k + 0, jj + 0);
            B_kkk_0_jjj_1 = MIX(B, LDB, k + 0, jj + 1);

            A_iii_0_kkk_0 = MIX(A, LDA, i + 0, k + 0);

            // ------

            C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
            C_iii_0_jjj_1 -= B_kkk_0_jjj_1 * A_iii_0_kkk_0;

            // -----

            MIX(C, LDC, i + 0, jj + 0) = C_iii_0_jjj_0;
            MIX(C, LDC, i + 0, jj + 1) = C_iii_0_jjj_1;
          }
      }

      // j overflow
      for (; j < N; ++j)
      {

        // k blocked
        for (k = 0; k <= K - NB; k += NB)

          // Cache blocking
          for (kk = k; kk <= (k + NB) - KU; kk += KU)

          // Register blocking
          /* for (kkk = kk; kkk < kk + KU; ++kkk) */
          {
            C_iii_0_jjj_0 = MIX(C, LDC, i + 0, j + 0);

            B_kkk_0_jjj_0 = MIX(B, LDB, kk + 0, j + 0);
            B_kkk_1_jjj_0 = MIX(B, LDB, kk + 1, j + 0);
            B_kkk_2_jjj_0 = MIX(B, LDB, kk + 2, j + 0);
            B_kkk_3_jjj_0 = MIX(B, LDB, kk + 3, j + 0);
            B_kkk_4_jjj_0 = MIX(B, LDB, kk + 4, j + 0);
            B_kkk_5_jjj_0 = MIX(B, LDB, kk + 5, j + 0);
            B_kkk_6_jjj_0 = MIX(B, LDB, kk + 6, j + 0);
            B_kkk_7_jjj_0 = MIX(B, LDB, kk + 7, j + 0);

            A_iii_0_kkk_0 = MIX(A, LDA, i + 0, kk + 0);
            A_iii_0_kkk_1 = MIX(A, LDA, i + 0, kk + 1);
            A_iii_0_kkk_2 = MIX(A, LDA, i + 0, kk + 2);
            A_iii_0_kkk_3 = MIX(A, LDA, i + 0, kk + 3);
            A_iii_0_kkk_4 = MIX(A, LDA, i + 0, kk + 4);
            A_iii_0_kkk_5 = MIX(A, LDA, i + 0, kk + 5);
            A_iii_0_kkk_6 = MIX(A, LDA, i + 0, kk + 6);
            A_iii_0_kkk_7 = MIX(A, LDA, i + 0, kk + 7);

            // ------

            C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
            C_iii_0_jjj_0 -= B_kkk_1_jjj_0 * A_iii_0_kkk_1;
            C_iii_0_jjj_0 -= B_kkk_2_jjj_0 * A_iii_0_kkk_2;
            C_iii_0_jjj_0 -= B_kkk_3_jjj_0 * A_iii_0_kkk_3;
            C_iii_0_jjj_0 -= B_kkk_4_jjj_0 * A_iii_0_kkk_4;
            C_iii_0_jjj_0 -= B_kkk_5_jjj_0 * A_iii_0_kkk_5;
            C_iii_0_jjj_0 -= B_kkk_6_jjj_0 * A_iii_0_kkk_6;
            C_iii_0_jjj_0 -= B_kkk_7_jjj_0 * A_iii_0_kkk_7;

            // -----

            MIX(C, LDC, i + 0, j + 0) = C_iii_0_jjj_0;
          }

        // k overflow
        for (; k < K; ++k)
          MIX(C, LDC, i, j) =
              MIX(C, LDC, i, j) - MIX(B, LDB, k, j) * MIX(A, LDA, i, k);
      }
    }
  }

  // Matrices where M < N
  // FIXME replace all inner unrolled loops with the version from above.
  else
  {

    // jj blocked
    for (j = 0; j <= N - NB; j += NB)
    {

      // i blocked
      for (i = 0; i <= M - NB; i += NB)
      {

        // k blocked
        for (k = 0; k <= K - NB; k += NB)

          // Cache blocked
          for (jj = j; jj <= (j + NB) - NU; jj += NU)
            for (ii = i; ii <= (i + NB) - MU; ii += MU)
              for (kk = k; kk <= (k + NB) - KU; kk += KU)

              // Register blocked
              /* for (kkk = kk; kkk < kk + KU; ++kkk) */
              /*   for (iii = ii; iii < ii + MU; ++iii) */
              /*     for (jjj = jj; jjj < jj + NU; ++jjj) */
              {

                C_iii_0_jjj_0 = MIX(C, LDC, ii + 0, jj + 0);
                C_iii_0_jjj_1 = MIX(C, LDC, ii + 0, jj + 1);
                C_iii_1_jjj_0 = MIX(C, LDC, ii + 1, jj + 0);
                C_iii_1_jjj_1 = MIX(C, LDC, ii + 1, jj + 1);
                C_iii_2_jjj_0 = MIX(C, LDC, ii + 2, jj + 0);
                C_iii_2_jjj_1 = MIX(C, LDC, ii + 2, jj + 1);
                C_iii_3_jjj_0 = MIX(C, LDC, ii + 3, jj + 0);
                C_iii_3_jjj_1 = MIX(C, LDC, ii + 3, jj + 1);

                B_kkk_0_jjj_0 = MIX(B, LDB, kk + 0, jj + 0);
                B_kkk_0_jjj_1 = MIX(B, LDB, kk + 0, jj + 1);
                B_kkk_1_jjj_0 = MIX(B, LDB, kk + 1, jj + 0);
                B_kkk_1_jjj_1 = MIX(B, LDB, kk + 1, jj + 1);
                B_kkk_2_jjj_0 = MIX(B, LDB, kk + 2, jj + 0);
                B_kkk_2_jjj_1 = MIX(B, LDB, kk + 2, jj + 1);
                B_kkk_3_jjj_0 = MIX(B, LDB, kk + 3, jj + 0);
                B_kkk_3_jjj_1 = MIX(B, LDB, kk + 3, jj + 1);
                B_kkk_4_jjj_0 = MIX(B, LDB, kk + 4, jj + 0);
                B_kkk_4_jjj_1 = MIX(B, LDB, kk + 4, jj + 1);
                B_kkk_5_jjj_0 = MIX(B, LDB, kk + 5, jj + 0);
                B_kkk_5_jjj_1 = MIX(B, LDB, kk + 5, jj + 1);
                B_kkk_6_jjj_0 = MIX(B, LDB, kk + 6, jj + 0);
                B_kkk_6_jjj_1 = MIX(B, LDB, kk + 6, jj + 1);
                B_kkk_7_jjj_0 = MIX(B, LDB, kk + 7, jj + 0);
                B_kkk_7_jjj_1 = MIX(B, LDB, kk + 7, jj + 1);

                A_iii_0_kkk_0 = MIX(A, LDA, ii + 0, kk + 0);
                A_iii_0_kkk_1 = MIX(A, LDA, ii + 0, kk + 1);
                A_iii_0_kkk_2 = MIX(A, LDA, ii + 0, kk + 2);
                A_iii_0_kkk_3 = MIX(A, LDA, ii + 0, kk + 3);
                A_iii_0_kkk_4 = MIX(A, LDA, ii + 0, kk + 4);
                A_iii_0_kkk_5 = MIX(A, LDA, ii + 0, kk + 5);
                A_iii_0_kkk_6 = MIX(A, LDA, ii + 0, kk + 6);
                A_iii_0_kkk_7 = MIX(A, LDA, ii + 0, kk + 7);

                A_iii_1_kkk_0 = MIX(A, LDA, ii + 1, kk + 0);
                A_iii_1_kkk_1 = MIX(A, LDA, ii + 1, kk + 1);
                A_iii_1_kkk_2 = MIX(A, LDA, ii + 1, kk + 2);
                A_iii_1_kkk_3 = MIX(A, LDA, ii + 1, kk + 3);
                A_iii_1_kkk_4 = MIX(A, LDA, ii + 1, kk + 4);
                A_iii_1_kkk_5 = MIX(A, LDA, ii + 1, kk + 5);
                A_iii_1_kkk_6 = MIX(A, LDA, ii + 1, kk + 6);
                A_iii_1_kkk_7 = MIX(A, LDA, ii + 1, kk + 7);

                A_iii_2_kkk_0 = MIX(A, LDA, ii + 2, kk + 0);
                A_iii_2_kkk_1 = MIX(A, LDA, ii + 2, kk + 1);
                A_iii_2_kkk_2 = MIX(A, LDA, ii + 2, kk + 2);
                A_iii_2_kkk_3 = MIX(A, LDA, ii + 2, kk + 3);
                A_iii_2_kkk_4 = MIX(A, LDA, ii + 2, kk + 4);
                A_iii_2_kkk_5 = MIX(A, LDA, ii + 2, kk + 5);
                A_iii_2_kkk_6 = MIX(A, LDA, ii + 2, kk + 6);
                A_iii_2_kkk_7 = MIX(A, LDA, ii + 2, kk + 7);

                A_iii_3_kkk_0 = MIX(A, LDA, ii + 3, kk + 0);
                A_iii_3_kkk_1 = MIX(A, LDA, ii + 3, kk + 1);
                A_iii_3_kkk_2 = MIX(A, LDA, ii + 3, kk + 2);
                A_iii_3_kkk_3 = MIX(A, LDA, ii + 3, kk + 3);
                A_iii_3_kkk_4 = MIX(A, LDA, ii + 3, kk + 4);
                A_iii_3_kkk_5 = MIX(A, LDA, ii + 3, kk + 5);
                A_iii_3_kkk_6 = MIX(A, LDA, ii + 3, kk + 6);
                A_iii_3_kkk_7 = MIX(A, LDA, ii + 3, kk + 7);
                A_iii_3_kkk_8 = MIX(A, LDA, ii + 3, kk + 8);

                // ------

                C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
                C_iii_0_jjj_0 -= B_kkk_1_jjj_0 * A_iii_0_kkk_1;
                C_iii_0_jjj_0 -= B_kkk_2_jjj_0 * A_iii_0_kkk_2;
                C_iii_0_jjj_0 -= B_kkk_3_jjj_0 * A_iii_0_kkk_3;
                C_iii_0_jjj_0 -= B_kkk_4_jjj_0 * A_iii_0_kkk_4;
                C_iii_0_jjj_0 -= B_kkk_5_jjj_0 * A_iii_0_kkk_5;
                C_iii_0_jjj_0 -= B_kkk_6_jjj_0 * A_iii_0_kkk_6;
                C_iii_0_jjj_0 -= B_kkk_7_jjj_0 * A_iii_0_kkk_7;

                C_iii_0_jjj_1 -= B_kkk_0_jjj_1 * A_iii_0_kkk_0;
                C_iii_0_jjj_1 -= B_kkk_1_jjj_1 * A_iii_0_kkk_1;
                C_iii_0_jjj_1 -= B_kkk_2_jjj_1 * A_iii_0_kkk_2;
                C_iii_0_jjj_1 -= B_kkk_3_jjj_1 * A_iii_0_kkk_3;
                C_iii_0_jjj_1 -= B_kkk_4_jjj_1 * A_iii_0_kkk_4;
                C_iii_0_jjj_1 -= B_kkk_5_jjj_1 * A_iii_0_kkk_5;
                C_iii_0_jjj_1 -= B_kkk_6_jjj_1 * A_iii_0_kkk_6;
                C_iii_0_jjj_1 -= B_kkk_7_jjj_1 * A_iii_0_kkk_7;

                C_iii_1_jjj_0 -= B_kkk_0_jjj_0 * A_iii_1_kkk_0;
                C_iii_1_jjj_0 -= B_kkk_1_jjj_0 * A_iii_1_kkk_1;
                C_iii_1_jjj_0 -= B_kkk_2_jjj_0 * A_iii_1_kkk_2;
                C_iii_1_jjj_0 -= B_kkk_3_jjj_0 * A_iii_1_kkk_3;
                C_iii_1_jjj_0 -= B_kkk_4_jjj_0 * A_iii_1_kkk_4;
                C_iii_1_jjj_0 -= B_kkk_5_jjj_0 * A_iii_1_kkk_5;
                C_iii_1_jjj_0 -= B_kkk_6_jjj_0 * A_iii_1_kkk_6;
                C_iii_1_jjj_0 -= B_kkk_7_jjj_0 * A_iii_1_kkk_7;

                C_iii_1_jjj_1 -= B_kkk_0_jjj_1 * A_iii_1_kkk_0;
                C_iii_1_jjj_1 -= B_kkk_1_jjj_1 * A_iii_1_kkk_1;
                C_iii_1_jjj_1 -= B_kkk_2_jjj_1 * A_iii_1_kkk_2;
                C_iii_1_jjj_1 -= B_kkk_3_jjj_1 * A_iii_1_kkk_3;
                C_iii_1_jjj_1 -= B_kkk_4_jjj_1 * A_iii_1_kkk_4;
                C_iii_1_jjj_1 -= B_kkk_5_jjj_1 * A_iii_1_kkk_5;
                C_iii_1_jjj_1 -= B_kkk_6_jjj_1 * A_iii_1_kkk_6;
                C_iii_1_jjj_1 -= B_kkk_7_jjj_1 * A_iii_1_kkk_7;

                C_iii_2_jjj_0 -= B_kkk_0_jjj_0 * A_iii_2_kkk_0;
                C_iii_2_jjj_0 -= B_kkk_1_jjj_0 * A_iii_2_kkk_1;
                C_iii_2_jjj_0 -= B_kkk_2_jjj_0 * A_iii_2_kkk_2;
                C_iii_2_jjj_0 -= B_kkk_3_jjj_0 * A_iii_2_kkk_3;
                C_iii_2_jjj_0 -= B_kkk_4_jjj_0 * A_iii_2_kkk_4;
                C_iii_2_jjj_0 -= B_kkk_5_jjj_0 * A_iii_2_kkk_5;
                C_iii_2_jjj_0 -= B_kkk_6_jjj_0 * A_iii_2_kkk_6;
                C_iii_2_jjj_0 -= B_kkk_7_jjj_0 * A_iii_2_kkk_7;

                C_iii_2_jjj_1 -= B_kkk_0_jjj_1 * A_iii_2_kkk_0;
                C_iii_2_jjj_1 -= B_kkk_1_jjj_1 * A_iii_2_kkk_1;
                C_iii_2_jjj_1 -= B_kkk_2_jjj_1 * A_iii_2_kkk_2;
                C_iii_2_jjj_1 -= B_kkk_3_jjj_1 * A_iii_2_kkk_3;
                C_iii_2_jjj_1 -= B_kkk_4_jjj_1 * A_iii_2_kkk_4;
                C_iii_2_jjj_1 -= B_kkk_5_jjj_1 * A_iii_2_kkk_5;
                C_iii_2_jjj_1 -= B_kkk_6_jjj_1 * A_iii_2_kkk_6;
                C_iii_2_jjj_1 -= B_kkk_7_jjj_1 * A_iii_2_kkk_7;

                C_iii_3_jjj_0 -= B_kkk_0_jjj_0 * A_iii_3_kkk_0;
                C_iii_3_jjj_0 -= B_kkk_1_jjj_0 * A_iii_3_kkk_1;
                C_iii_3_jjj_0 -= B_kkk_2_jjj_0 * A_iii_3_kkk_2;
                C_iii_3_jjj_0 -= B_kkk_3_jjj_0 * A_iii_3_kkk_3;
                C_iii_3_jjj_0 -= B_kkk_4_jjj_0 * A_iii_3_kkk_4;
                C_iii_3_jjj_0 -= B_kkk_5_jjj_0 * A_iii_3_kkk_5;
                C_iii_3_jjj_0 -= B_kkk_6_jjj_0 * A_iii_3_kkk_6;
                C_iii_3_jjj_0 -= B_kkk_7_jjj_0 * A_iii_3_kkk_7;

                C_iii_3_jjj_1 -= B_kkk_0_jjj_1 * A_iii_3_kkk_0;
                C_iii_3_jjj_1 -= B_kkk_1_jjj_1 * A_iii_3_kkk_1;
                C_iii_3_jjj_1 -= B_kkk_2_jjj_1 * A_iii_3_kkk_2;
                C_iii_3_jjj_1 -= B_kkk_3_jjj_1 * A_iii_3_kkk_3;
                C_iii_3_jjj_1 -= B_kkk_4_jjj_1 * A_iii_3_kkk_4;
                C_iii_3_jjj_1 -= B_kkk_5_jjj_1 * A_iii_3_kkk_5;
                C_iii_3_jjj_1 -= B_kkk_6_jjj_1 * A_iii_3_kkk_6;
                C_iii_3_jjj_1 -= B_kkk_7_jjj_1 * A_iii_3_kkk_7;

                // -----

                MIX(C, LDC, ii + 0, jj + 0) = C_iii_0_jjj_0;
                MIX(C, LDC, ii + 0, jj + 1) = C_iii_0_jjj_1;
                MIX(C, LDC, ii + 1, jj + 0) = C_iii_1_jjj_0;
                MIX(C, LDC, ii + 1, jj + 1) = C_iii_1_jjj_1;
                MIX(C, LDC, ii + 2, jj + 0) = C_iii_2_jjj_0;
                MIX(C, LDC, ii + 2, jj + 1) = C_iii_2_jjj_1;
                MIX(C, LDC, ii + 3, jj + 0) = C_iii_3_jjj_0;
                MIX(C, LDC, ii + 3, jj + 1) = C_iii_3_jjj_1;
              }

        // k overflow
        for (; k < K; ++k)

          // Cache blocked
          for (jj = j; jj <= (j + NB) - NU; jj += NU)
            for (ii = i; ii <= (i + NB) - MU; ii += MU)

            // Register blocked
            /* for (iii = ii; iii < ii + MU; ++iii) */
            /*   for (jjj = jj; jjj < jj + NU; ++jjj) */
            {

              C_iii_0_jjj_0 = MIX(C, LDC, ii + 0, jj + 0);
              C_iii_0_jjj_1 = MIX(C, LDC, ii + 0, jj + 1);

              C_iii_1_jjj_0 = MIX(C, LDC, ii + 1, jj + 0);
              C_iii_1_jjj_1 = MIX(C, LDC, ii + 1, jj + 1);

              C_iii_2_jjj_0 = MIX(C, LDC, ii + 2, jj + 0);
              C_iii_2_jjj_1 = MIX(C, LDC, ii + 2, jj + 1);

              C_iii_3_jjj_0 = MIX(C, LDC, ii + 3, jj + 0);
              C_iii_3_jjj_1 = MIX(C, LDC, ii + 3, jj + 1);

              B_kkk_0_jjj_0 = MIX(B, LDB, k + 0, jj + 0);
              B_kkk_0_jjj_1 = MIX(B, LDB, k + 0, jj + 1);

              A_iii_0_kkk_0 = MIX(A, LDA, ii + 0, k + 0);
              A_iii_1_kkk_0 = MIX(A, LDA, ii + 1, k + 0);
              A_iii_2_kkk_0 = MIX(A, LDA, ii + 2, k + 0);
              A_iii_3_kkk_0 = MIX(A, LDA, ii + 3, k + 0);

              // ----

              C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
              C_iii_0_jjj_1 -= B_kkk_0_jjj_1 * A_iii_0_kkk_0;
              C_iii_1_jjj_0 -= B_kkk_0_jjj_0 * A_iii_1_kkk_0;
              C_iii_1_jjj_1 -= B_kkk_0_jjj_1 * A_iii_1_kkk_0;
              C_iii_2_jjj_0 -= B_kkk_0_jjj_0 * A_iii_2_kkk_0;
              C_iii_2_jjj_1 -= B_kkk_0_jjj_1 * A_iii_2_kkk_0;
              C_iii_3_jjj_0 -= B_kkk_0_jjj_0 * A_iii_3_kkk_0;
              C_iii_3_jjj_1 -= B_kkk_0_jjj_1 * A_iii_3_kkk_0;

              // ----

              MIX(C, LDC, ii + 0, jj + 0) = C_iii_0_jjj_0;
              MIX(C, LDC, ii + 0, jj + 1) = C_iii_0_jjj_1;
              MIX(C, LDC, ii + 1, jj + 0) = C_iii_1_jjj_0;
              MIX(C, LDC, ii + 1, jj + 1) = C_iii_1_jjj_1;
              MIX(C, LDC, ii + 2, jj + 0) = C_iii_2_jjj_0;
              MIX(C, LDC, ii + 2, jj + 1) = C_iii_2_jjj_1;
              MIX(C, LDC, ii + 3, jj + 0) = C_iii_3_jjj_0;
              MIX(C, LDC, ii + 3, jj + 1) = C_iii_3_jjj_1;
            }
      }

      // i overflow
      for (; i < M; ++i)
      {

        // K blocked
        for (k = 0; k <= K - NB; k += NB)

          // Cache blocked
          for (jj = j; jj <= (j + NB) - NU; jj += NU)
            for (kk = k; kk <= (k + NB) - KU; kk += KU)

            // Register blocked
            /* for (jjj = jj; jjj < jj + NU; ++jjj) */
            /*   for (kkk = kk; kkk < kk + KU; ++kkk) */
            {

              C_iii_0_jjj_0 = MIX(C, LDC, i + 0, jj + 0);
              C_iii_0_jjj_1 = MIX(C, LDC, i + 0, jj + 1);

              B_kkk_0_jjj_0 = MIX(B, LDB, kk + 0, jj + 0);
              B_kkk_0_jjj_1 = MIX(B, LDB, kk + 0, jj + 1);
              B_kkk_1_jjj_0 = MIX(B, LDB, kk + 1, jj + 0);
              B_kkk_1_jjj_1 = MIX(B, LDB, kk + 1, jj + 1);
              B_kkk_2_jjj_0 = MIX(B, LDB, kk + 2, jj + 0);
              B_kkk_2_jjj_1 = MIX(B, LDB, kk + 2, jj + 1);
              B_kkk_3_jjj_0 = MIX(B, LDB, kk + 3, jj + 0);
              B_kkk_3_jjj_1 = MIX(B, LDB, kk + 3, jj + 1);
              B_kkk_4_jjj_0 = MIX(B, LDB, kk + 4, jj + 0);
              B_kkk_4_jjj_1 = MIX(B, LDB, kk + 4, jj + 1);
              B_kkk_5_jjj_0 = MIX(B, LDB, kk + 5, jj + 0);
              B_kkk_5_jjj_1 = MIX(B, LDB, kk + 5, jj + 1);
              B_kkk_6_jjj_0 = MIX(B, LDB, kk + 6, jj + 0);
              B_kkk_6_jjj_1 = MIX(B, LDB, kk + 6, jj + 1);
              B_kkk_7_jjj_0 = MIX(B, LDB, kk + 7, jj + 0);
              B_kkk_7_jjj_1 = MIX(B, LDB, kk + 7, jj + 1);

              A_iii_0_kkk_0 = MIX(A, LDA, i + 0, kk + 0);
              A_iii_0_kkk_1 = MIX(A, LDA, i + 0, kk + 1);
              A_iii_0_kkk_2 = MIX(A, LDA, i + 0, kk + 2);
              A_iii_0_kkk_3 = MIX(A, LDA, i + 0, kk + 3);
              A_iii_0_kkk_4 = MIX(A, LDA, i + 0, kk + 4);
              A_iii_0_kkk_5 = MIX(A, LDA, i + 0, kk + 5);
              A_iii_0_kkk_6 = MIX(A, LDA, i + 0, kk + 6);
              A_iii_0_kkk_7 = MIX(A, LDA, i + 0, kk + 7);

              // ------

              C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
              C_iii_0_jjj_0 -= B_kkk_1_jjj_0 * A_iii_0_kkk_1;
              C_iii_0_jjj_0 -= B_kkk_2_jjj_0 * A_iii_0_kkk_2;
              C_iii_0_jjj_0 -= B_kkk_3_jjj_0 * A_iii_0_kkk_3;
              C_iii_0_jjj_0 -= B_kkk_4_jjj_0 * A_iii_0_kkk_4;
              C_iii_0_jjj_0 -= B_kkk_5_jjj_0 * A_iii_0_kkk_5;
              C_iii_0_jjj_0 -= B_kkk_6_jjj_0 * A_iii_0_kkk_6;
              C_iii_0_jjj_0 -= B_kkk_7_jjj_0 * A_iii_0_kkk_7;

              C_iii_0_jjj_1 -= B_kkk_0_jjj_1 * A_iii_0_kkk_0;
              C_iii_0_jjj_1 -= B_kkk_1_jjj_1 * A_iii_0_kkk_1;
              C_iii_0_jjj_1 -= B_kkk_2_jjj_1 * A_iii_0_kkk_2;
              C_iii_0_jjj_1 -= B_kkk_3_jjj_1 * A_iii_0_kkk_3;
              C_iii_0_jjj_1 -= B_kkk_4_jjj_1 * A_iii_0_kkk_4;
              C_iii_0_jjj_1 -= B_kkk_5_jjj_1 * A_iii_0_kkk_5;
              C_iii_0_jjj_1 -= B_kkk_6_jjj_1 * A_iii_0_kkk_6;
              C_iii_0_jjj_1 -= B_kkk_7_jjj_1 * A_iii_0_kkk_7;

              // -----

              MIX(C, LDC, i + 0, jj + 0) = C_iii_0_jjj_0;
              MIX(C, LDC, i + 0, jj + 1) = C_iii_0_jjj_1;
            }

        // K overflow
        for (; k < K; ++k)

          // Cache blocking
          for (jj = j; jj <= (j + NB) - NU; jj += NU)

          // Register blocking
          /* for (jjj = jj; jjj < jj + NU; ++jjj) */
          {

            C_iii_0_jjj_0 = MIX(C, LDC, i + 0, jj + 0);
            C_iii_0_jjj_1 = MIX(C, LDC, i + 0, jj + 1);

            B_kkk_0_jjj_0 = MIX(B, LDB, k + 0, jj + 0);
            B_kkk_0_jjj_1 = MIX(B, LDB, k + 0, jj + 1);

            A_iii_0_kkk_0 = MIX(A, LDA, i + 0, k + 0);

            // ------

            C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
            C_iii_0_jjj_1 -= B_kkk_0_jjj_1 * A_iii_0_kkk_0;

            // -----

            MIX(C, LDC, i + 0, jj + 0) = C_iii_0_jjj_0;
            MIX(C, LDC, i + 0, jj + 1) = C_iii_0_jjj_1;
          }
      }
    }

    // j overflow
    for (; j < N; ++j)
    {

      // i blocked
      for (i = 0; i <= M - NB; i += NB)
      {

        // k blocked
        for (k = 0; k <= K - NB; k += NB)

          // Cache blocking
          for (ii = i; ii <= (i + NB) - MU; ii += MU)
            for (kk = k; kk <= (k + NB) - KU; kk += KU)

            // Register blocking
            /* for (kkk = kk; kkk < kk + KU; ++kkk) */
            /*   for (iii = ii; iii < ii + MU; ++iii) */
            {

              C_iii_0_jjj_0 = MIX(C, LDC, ii + 0, j + 0);
              C_iii_1_jjj_0 = MIX(C, LDC, ii + 1, j + 0);
              C_iii_2_jjj_0 = MIX(C, LDC, ii + 2, j + 0);
              C_iii_3_jjj_0 = MIX(C, LDC, ii + 3, j + 0);

              B_kkk_0_jjj_0 = MIX(B, LDB, kk + 0, j + 0);
              B_kkk_1_jjj_0 = MIX(B, LDB, kk + 1, j + 0);
              B_kkk_2_jjj_0 = MIX(B, LDB, kk + 2, j + 0);
              B_kkk_3_jjj_0 = MIX(B, LDB, kk + 3, j + 0);
              B_kkk_4_jjj_0 = MIX(B, LDB, kk + 4, j + 0);
              B_kkk_5_jjj_0 = MIX(B, LDB, kk + 5, j + 0);
              B_kkk_6_jjj_0 = MIX(B, LDB, kk + 6, j + 0);
              B_kkk_7_jjj_0 = MIX(B, LDB, kk + 7, j + 0);

              A_iii_0_kkk_0 = MIX(A, LDA, ii + 0, kk + 0);
              A_iii_0_kkk_1 = MIX(A, LDA, ii + 0, kk + 1);
              A_iii_0_kkk_2 = MIX(A, LDA, ii + 0, kk + 2);
              A_iii_0_kkk_3 = MIX(A, LDA, ii + 0, kk + 3);
              A_iii_0_kkk_4 = MIX(A, LDA, ii + 0, kk + 4);
              A_iii_0_kkk_5 = MIX(A, LDA, ii + 0, kk + 5);
              A_iii_0_kkk_6 = MIX(A, LDA, ii + 0, kk + 6);
              A_iii_0_kkk_7 = MIX(A, LDA, ii + 0, kk + 7);

              A_iii_1_kkk_0 = MIX(A, LDA, ii + 1, kk + 0);
              A_iii_1_kkk_1 = MIX(A, LDA, ii + 1, kk + 1);
              A_iii_1_kkk_2 = MIX(A, LDA, ii + 1, kk + 2);
              A_iii_1_kkk_3 = MIX(A, LDA, ii + 1, kk + 3);
              A_iii_1_kkk_4 = MIX(A, LDA, ii + 1, kk + 4);
              A_iii_1_kkk_5 = MIX(A, LDA, ii + 1, kk + 5);
              A_iii_1_kkk_6 = MIX(A, LDA, ii + 1, kk + 6);
              A_iii_1_kkk_7 = MIX(A, LDA, ii + 1, kk + 7);

              A_iii_2_kkk_0 = MIX(A, LDA, ii + 2, kk + 0);
              A_iii_2_kkk_1 = MIX(A, LDA, ii + 2, kk + 1);
              A_iii_2_kkk_2 = MIX(A, LDA, ii + 2, kk + 2);
              A_iii_2_kkk_3 = MIX(A, LDA, ii + 2, kk + 3);
              A_iii_2_kkk_4 = MIX(A, LDA, ii + 2, kk + 4);
              A_iii_2_kkk_5 = MIX(A, LDA, ii + 2, kk + 5);
              A_iii_2_kkk_6 = MIX(A, LDA, ii + 2, kk + 6);
              A_iii_2_kkk_7 = MIX(A, LDA, ii + 2, kk + 7);

              A_iii_3_kkk_0 = MIX(A, LDA, ii + 3, kk + 0);
              A_iii_3_kkk_1 = MIX(A, LDA, ii + 3, kk + 1);
              A_iii_3_kkk_2 = MIX(A, LDA, ii + 3, kk + 2);
              A_iii_3_kkk_3 = MIX(A, LDA, ii + 3, kk + 3);
              A_iii_3_kkk_4 = MIX(A, LDA, ii + 3, kk + 4);
              A_iii_3_kkk_5 = MIX(A, LDA, ii + 3, kk + 5);
              A_iii_3_kkk_6 = MIX(A, LDA, ii + 3, kk + 6);
              A_iii_3_kkk_7 = MIX(A, LDA, ii + 3, kk + 7);
              A_iii_3_kkk_8 = MIX(A, LDA, ii + 3, kk + 8);

              // ------

              C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
              C_iii_0_jjj_0 -= B_kkk_1_jjj_0 * A_iii_0_kkk_1;
              C_iii_0_jjj_0 -= B_kkk_2_jjj_0 * A_iii_0_kkk_2;
              C_iii_0_jjj_0 -= B_kkk_3_jjj_0 * A_iii_0_kkk_3;
              C_iii_0_jjj_0 -= B_kkk_4_jjj_0 * A_iii_0_kkk_4;
              C_iii_0_jjj_0 -= B_kkk_5_jjj_0 * A_iii_0_kkk_5;
              C_iii_0_jjj_0 -= B_kkk_6_jjj_0 * A_iii_0_kkk_6;
              C_iii_0_jjj_0 -= B_kkk_7_jjj_0 * A_iii_0_kkk_7;

              C_iii_1_jjj_0 -= B_kkk_0_jjj_0 * A_iii_1_kkk_0;
              C_iii_1_jjj_0 -= B_kkk_1_jjj_0 * A_iii_1_kkk_1;
              C_iii_1_jjj_0 -= B_kkk_2_jjj_0 * A_iii_1_kkk_2;
              C_iii_1_jjj_0 -= B_kkk_3_jjj_0 * A_iii_1_kkk_3;
              C_iii_1_jjj_0 -= B_kkk_4_jjj_0 * A_iii_1_kkk_4;
              C_iii_1_jjj_0 -= B_kkk_5_jjj_0 * A_iii_1_kkk_5;
              C_iii_1_jjj_0 -= B_kkk_6_jjj_0 * A_iii_1_kkk_6;
              C_iii_1_jjj_0 -= B_kkk_7_jjj_0 * A_iii_1_kkk_7;

              C_iii_2_jjj_0 -= B_kkk_0_jjj_0 * A_iii_2_kkk_0;
              C_iii_2_jjj_0 -= B_kkk_1_jjj_0 * A_iii_2_kkk_1;
              C_iii_2_jjj_0 -= B_kkk_2_jjj_0 * A_iii_2_kkk_2;
              C_iii_2_jjj_0 -= B_kkk_3_jjj_0 * A_iii_2_kkk_3;
              C_iii_2_jjj_0 -= B_kkk_4_jjj_0 * A_iii_2_kkk_4;
              C_iii_2_jjj_0 -= B_kkk_5_jjj_0 * A_iii_2_kkk_5;
              C_iii_2_jjj_0 -= B_kkk_6_jjj_0 * A_iii_2_kkk_6;
              C_iii_2_jjj_0 -= B_kkk_7_jjj_0 * A_iii_2_kkk_7;

              C_iii_3_jjj_0 -= B_kkk_0_jjj_0 * A_iii_3_kkk_0;
              C_iii_3_jjj_0 -= B_kkk_1_jjj_0 * A_iii_3_kkk_1;
              C_iii_3_jjj_0 -= B_kkk_2_jjj_0 * A_iii_3_kkk_2;
              C_iii_3_jjj_0 -= B_kkk_3_jjj_0 * A_iii_3_kkk_3;
              C_iii_3_jjj_0 -= B_kkk_4_jjj_0 * A_iii_3_kkk_4;
              C_iii_3_jjj_0 -= B_kkk_5_jjj_0 * A_iii_3_kkk_5;
              C_iii_3_jjj_0 -= B_kkk_6_jjj_0 * A_iii_3_kkk_6;
              C_iii_3_jjj_0 -= B_kkk_7_jjj_0 * A_iii_3_kkk_7;

              // -----

              MIX(C, LDC, ii + 0, j + 0) = C_iii_0_jjj_0;
              MIX(C, LDC, ii + 1, j + 0) = C_iii_1_jjj_0;
              MIX(C, LDC, ii + 2, j + 0) = C_iii_2_jjj_0;
              MIX(C, LDC, ii + 3, j + 0) = C_iii_3_jjj_0;
            }

        // k overflow
        for (; k < K; ++k)

          // Cache blocking
          for (ii = i; ii <= (i + NB) - MU; ii += MU)

          // Register blocking
          /* for (iii = ii; iii < ii + MU; ++iii) */
          {
            C_iii_0_jjj_0 = MIX(C, LDC, ii + 0, j + 0);
            C_iii_1_jjj_0 = MIX(C, LDC, ii + 1, j + 0);
            C_iii_2_jjj_0 = MIX(C, LDC, ii + 2, j + 0);
            C_iii_3_jjj_0 = MIX(C, LDC, ii + 3, j + 0);

            B_kkk_0_jjj_0 = MIX(B, LDB, k + 0, j + 0);

            A_iii_0_kkk_0 = MIX(A, LDA, ii + 0, k + 0);
            A_iii_1_kkk_0 = MIX(A, LDA, ii + 1, k + 0);
            A_iii_2_kkk_0 = MIX(A, LDA, ii + 2, k + 0);
            A_iii_3_kkk_0 = MIX(A, LDA, ii + 3, k + 0);

            // ------

            C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
            C_iii_1_jjj_0 -= B_kkk_0_jjj_0 * A_iii_1_kkk_0;
            C_iii_2_jjj_0 -= B_kkk_0_jjj_0 * A_iii_2_kkk_0;
            C_iii_3_jjj_0 -= B_kkk_0_jjj_0 * A_iii_3_kkk_0;

            // -----

            MIX(C, LDC, ii + 0, j + 0) = C_iii_0_jjj_0;
            MIX(C, LDC, ii + 1, j + 0) = C_iii_1_jjj_0;
            MIX(C, LDC, ii + 2, j + 0) = C_iii_2_jjj_0;
            MIX(C, LDC, ii + 3, j + 0) = C_iii_3_jjj_0;
          }
      }

      // i overflow
      for (; i < M; ++i)
      {

        // K blocked
        for (k = 0; k <= K - NB; k += NB)

          // Cache blocking
          for (kk = k; kk <= (k + NB) - KU; kk += KU)

          // Register blocking
          /* for (kkk = kk; kkk < kk + KU; ++kkk) */
          {

            C_iii_0_jjj_0 = MIX(C, LDC, i + 0, j + 0);

            B_kkk_0_jjj_0 = MIX(B, LDB, kk + 0, j + 0);
            B_kkk_1_jjj_0 = MIX(B, LDB, kk + 1, j + 0);
            B_kkk_2_jjj_0 = MIX(B, LDB, kk + 2, j + 0);
            B_kkk_3_jjj_0 = MIX(B, LDB, kk + 3, j + 0);
            B_kkk_4_jjj_0 = MIX(B, LDB, kk + 4, j + 0);
            B_kkk_5_jjj_0 = MIX(B, LDB, kk + 5, j + 0);
            B_kkk_6_jjj_0 = MIX(B, LDB, kk + 6, j + 0);
            B_kkk_7_jjj_0 = MIX(B, LDB, kk + 7, j + 0);

            A_iii_0_kkk_0 = MIX(A, LDA, i + 0, kk + 0);
            A_iii_0_kkk_1 = MIX(A, LDA, i + 0, kk + 1);
            A_iii_0_kkk_2 = MIX(A, LDA, i + 0, kk + 2);
            A_iii_0_kkk_3 = MIX(A, LDA, i + 0, kk + 3);
            A_iii_0_kkk_4 = MIX(A, LDA, i + 0, kk + 4);
            A_iii_0_kkk_5 = MIX(A, LDA, i + 0, kk + 5);
            A_iii_0_kkk_6 = MIX(A, LDA, i + 0, kk + 6);
            A_iii_0_kkk_7 = MIX(A, LDA, i + 0, kk + 7);

            // ------

            C_iii_0_jjj_0 -= B_kkk_0_jjj_0 * A_iii_0_kkk_0;
            C_iii_0_jjj_0 -= B_kkk_1_jjj_0 * A_iii_0_kkk_1;
            C_iii_0_jjj_0 -= B_kkk_2_jjj_0 * A_iii_0_kkk_2;
            C_iii_0_jjj_0 -= B_kkk_3_jjj_0 * A_iii_0_kkk_3;
            C_iii_0_jjj_0 -= B_kkk_4_jjj_0 * A_iii_0_kkk_4;
            C_iii_0_jjj_0 -= B_kkk_5_jjj_0 * A_iii_0_kkk_5;
            C_iii_0_jjj_0 -= B_kkk_6_jjj_0 * A_iii_0_kkk_6;
            C_iii_0_jjj_0 -= B_kkk_7_jjj_0 * A_iii_0_kkk_7;

            // -----

            MIX(C, LDC, i + 0, j + 0) = C_iii_0_jjj_0;
          }

        // K overflow
        for (; k < K; ++k)
          MIX(C, LDC, i, j) =
              MIX(C, LDC, i, j) - MIX(B, LDB, k, j) * MIX(A, LDA, i, k);
      }
    }
  }
  return;
}

int sgetf2_2(int M, int N, double *A, int LDA, int *ipiv)
{

  int i, j, k, p_i;

  double   //
      p_v, //
      m_0  //
      ;

  double      //
      A_j_k0, //
      A_j_k1, //
      A_j_k2, //
      A_j_k3, //
      A_j_k4, //
      A_j_k5, //
      A_j_k6, //
      A_j_k7, //

      A_i_k0, //
      A_i_k1, //
      A_i_k2, //
      A_i_k3, //
      A_i_k4, //
      A_i_k5, //
      A_i_k6, //
      A_i_k7  //
      ;
  // Quick return
  if (!M || !N)
    return 0;

  // NOTE I have increased this to iterate *until* N however you could
  // stop at < N - 1 if you cover the bounds case.
  for (i = 0; i < MIN(M, N); ++i)
  {

    p_i = i;
    p_v = MIX(A, LDA, p_i, i);

    // Only try to find a pivot if  the current value is approaching 0
    if (APPROX_EQUAL(p_v, 0.)) // FIXME
    {

      p_i = i + isamax_2(M - i, &MIX(A, LDA, i, i), LDA);
      p_v = MIX(A, LDA, p_i, i);

      if (APPROX_EQUAL(p_v, 0.))
      {
        fprintf(stderr, "ERROR: LU Solve singular matrix\n");
        return -1;
      }
    }

    ipiv[i] = p_i;

    if (i != p_i)
    {
      sswap_2(N, &MIX(A, LDA, i, 0), 1, &MIX(A, LDA, p_i, 0), 1);
    }

    // TODO compute the minimum machine safe integer such that
    // 1 / fmin doesn't overflow. This can be inlind because
    // we are targeting SKYLAKE.

    m_0 = 1 / MIX(A, LDA, i, i);

    // BLAS 1 Scale vector (NOTE column vector)
    for (j = i + 1; j < M; ++j)
    {
      MIX(A, LDA, j, i) = m_0 * MIX(A, LDA, j, i);
    }

    if (i < MIN(M, N))
    {

      // BLAS 2 rank-1 update
      for (j = i + 1; j < M; ++j)
      {

        // Negate to make computations look like FMA :)
        m_0 = -MIX(A, LDA, j, i);
        for (k = i + 1; k < N - 7; k += 8)
        {
          A_j_k0 = MIX(A, LDA, j, k + 0);
          A_j_k1 = MIX(A, LDA, j, k + 1);
          A_j_k2 = MIX(A, LDA, j, k + 2);
          A_j_k3 = MIX(A, LDA, j, k + 3);
          A_j_k4 = MIX(A, LDA, j, k + 4);
          A_j_k5 = MIX(A, LDA, j, k + 5);
          A_j_k6 = MIX(A, LDA, j, k + 6);
          A_j_k7 = MIX(A, LDA, j, k + 7);

          A_i_k0 = MIX(A, LDA, i, k + 0);
          A_i_k1 = MIX(A, LDA, i, k + 1);
          A_i_k2 = MIX(A, LDA, i, k + 2);
          A_i_k3 = MIX(A, LDA, i, k + 3);
          A_i_k4 = MIX(A, LDA, i, k + 4);
          A_i_k5 = MIX(A, LDA, i, k + 5);
          A_i_k6 = MIX(A, LDA, i, k + 6);
          A_i_k7 = MIX(A, LDA, i, k + 7);

          MIX(A, LDA, j, k + 0) = m_0 * A_i_k0 + A_j_k0;
          MIX(A, LDA, j, k + 1) = m_0 * A_i_k1 + A_j_k1;
          MIX(A, LDA, j, k + 2) = m_0 * A_i_k2 + A_j_k2;
          MIX(A, LDA, j, k + 3) = m_0 * A_i_k3 + A_j_k3;
          MIX(A, LDA, j, k + 4) = m_0 * A_i_k4 + A_j_k4;
          MIX(A, LDA, j, k + 5) = m_0 * A_i_k5 + A_j_k5;
          MIX(A, LDA, j, k + 6) = m_0 * A_i_k6 + A_j_k6;
          MIX(A, LDA, j, k + 7) = m_0 * A_i_k7 + A_j_k7;
        }

        for (; k < N; ++k)
        {
          A_j_k0 = MIX(A, LDA, j, k + 0);
          A_i_k0 = MIX(A, LDA, i, k + 0);
          MIX(A, LDA, j, k + 0) = m_0 * A_i_k0 + A_j_k0;
        }
      }
    }
  }

  return 0;
}

// Equivalent to SGETRS (no-transpose).
static int sgetrs_2(int N, double *A, int *ipiv, double *b)
{
  // A now contains L (below diagonal)
  //                U (above diagonal)
  // Swap pivot rows in b
  slaswp_2(1, b, 1, 0, N, ipiv, 1);
  // Forward substitution
  strsm_L_2(N, 1, A, N, b, 1);
  // Backward substitution
  strsm_U_2(N, 1, A, N, b, 1);
  return 0;
}

// Skylake has L1 cache of 32 KB
int lu_solve_2(int N, double *A, int *ipiv, double *b)
{
  int retcode, ib, IB, k;

  const int NB = ideal_block(N, N), //
      M = N,                        //
      LDA = N,                      //
      MIN_MN = N                    //
      ;

  // Use unblocked code
  if (NB <= 1 || NB >= MIN_MN)
  {
    retcode = sgetf2_2(M, N, A, LDA, ipiv);
    if (retcode != 0)
      return retcode;
  }

  // BLocked factor A into [L \ U]
  else
  {
    for (ib = 0; ib < MIN_MN; ib += NB)
    {
      IB = MIN(MIN_MN - ib, NB);

      retcode = sgetf2_2(M - ib, IB, &AIX(ib, ib), LDA, ipiv + ib);

      if (retcode != 0)
        return retcode;

      // Update the pivot indices
      for (k = ib; k < MIN(M, ib + IB) - 7; k += 8)
      {
        ipiv[k + 0] = ipiv[k + 0] + ib;
        ipiv[k + 1] = ipiv[k + 1] + ib;
        ipiv[k + 2] = ipiv[k + 2] + ib;
        ipiv[k + 3] = ipiv[k + 3] + ib;

        ipiv[k + 4] = ipiv[k + 4] + ib;
        ipiv[k + 5] = ipiv[k + 5] + ib;
        ipiv[k + 6] = ipiv[k + 6] + ib;
        ipiv[k + 7] = ipiv[k + 7] + ib;
      }
      for (; k < MIN(M, ib + IB); ++k)
      {
        ipiv[k + 0] = ipiv[k + 0] + ib;
      }

      // Apply interchanges to columns 0 : ib
      slaswp_2(ib, A, LDA, ib, ib + IB, ipiv, 1);

      if (ib + IB < N)
      {
        // Apply interchanges to columns ib + IB : N
        slaswp_2(N - ib - IB, &AIX(0, ib + IB), LDA, ib, ib + IB, ipiv, 1);

        // Compute the block row of U
        strsm_L_2(IB, N - ib - IB, &AIX(ib, ib), LDA, &AIX(ib, ib + IB), LDA);

        // Update trailing submatrix
        sgemm_2(M - ib - IB, N - ib - IB, IB, -ONE, //
                &AIX(ib + IB, ib), LDA,             //
                &AIX(ib, ib + IB), LDA,             //
                ONE,                                //
                &AIX(ib + IB, ib + IB), LDA         //
        );
      }
    }
  }

  // Solve the system with A
  retcode = sgetrs_2(N, A, ipiv, b);

  return retcode;
}

#ifdef TEST_MKL

/**
 * This implementation should remain /equal/ to the previous, however,
 * it uses the Intel OneAPI MKL instead of the handrolled LAPACK routines.
 */
int lu_solve_3(int N, double *A, int *ipiv, double *b)
{
  int retcode, ib, IB, k;

  const int //
      NB = ideal_block(N, N),
      M = N,     //
      LDA = N,   //
      MIN_MN = N //
      ;

  // BLocked factor A into [L \ U]

  if (NB <= 1 || NB >= MIN_MN)
  {
    // Use unblocked code
    retcode = LAPACKE_dgetf2(LAPACK_ROW_MAJOR, M, N, A, LDA, ipiv);
    if (retcode != 0)
      return retcode;
  }
  else
  {
    for (ib = 0; ib < MIN_MN; ib += NB)
    {
      IB = MIN(MIN_MN - ib, NB);

      retcode = LAPACKE_dgetf2(LAPACK_ROW_MAJOR, M - ib, IB, &AIX(ib, ib), LDA,
                               ipiv + ib);

      if (retcode != 0)
        return retcode;

      // Update the pivot indices
      for (k = ib; k < MIN(M, ib + IB) - 7; k += 8)
      {
        ipiv[k + 0] = ipiv[k + 0] + ib;
        ipiv[k + 1] = ipiv[k + 1] + ib;
        ipiv[k + 2] = ipiv[k + 2] + ib;
        ipiv[k + 3] = ipiv[k + 3] + ib;

        ipiv[k + 4] = ipiv[k + 4] + ib;
        ipiv[k + 5] = ipiv[k + 5] + ib;
        ipiv[k + 6] = ipiv[k + 6] + ib;
        ipiv[k + 7] = ipiv[k + 7] + ib;
      }
      for (; k < MIN(M, ib + IB); ++k)
      {
        ipiv[k + 0] = ipiv[k + 0] + ib;
      }

      // Apply interchanges to columns 0 : ib
      slaswp_2(ib, A, LDA, ib, ib + IB, ipiv, 1);

      if (ib + IB < N)
      {
        // Apply interchanges to columns ib + IB : N
        slaswp_2(N - ib - IB, &AIX(0, ib + IB), LDA, ib, ib + IB, ipiv, 1);

        // Compute the block row of U
        strsm_L_2(IB, N - ib - IB, &AIX(ib, ib), LDA, &AIX(ib, ib + IB), LDA);

        // Update trailing submatrix
        sgemm_intel(M - ib - IB, N - ib - IB, IB, -ONE, //
                    &AIX(ib + IB, ib), LDA,             //
                    &AIX(ib, ib + IB), LDA,             //
                    ONE,                                //
                    &AIX(ib + IB, ib + IB), LDA         //
        );
      }
    }
  }

  // Solve the system with A
  // retcode = sgetrs_2(N, A, ipiv, b);
  retcode = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N',
                           N,    // Number of equations
                           1,    // Number of rhs equations
                           A, N, //
                           ipiv, //
                           b, 1  //
  );

  return retcode;
}

/**
 * Solve system only with the Intel OneAPI routine.
 */
int lu_solve_4(int N, double *A, int *ipiv, double *b)
{
  return LAPACKE_dgesv(LAPACK_ROW_MAJOR,
                       N,    // Number of equations
                       1,    // Number of rhs equations
                       A, N, //
                       ipiv, //
                       b, 1  //
  );
}

#endif // TEST_MKL

#ifdef TEST_PERF

void register_functions_LU_SOLVE()
{
  // add_function_LU_SOLVE(&lu_solve_0, "LU Solve Base", 1);
  // add_function_LU_SOLVE(&lu_solve_1, "LU Solve Recursive", 1);
  add_function_LU_SOLVE(&lu_solve_2, "LU Solve Basic C Opts", 1);
#ifdef TEST_MKL
  add_function_LU_SOLVE(&lu_solve_3, "LU Solve Intel DGEMM", 1);
  add_function_LU_SOLVE(&lu_solve_4, "Intel DGESV Row Major", 1);
#endif
}

void register_functions_MMM()
{
  add_function_MMM(&sgemm_1, "MMM Base", 1);
  add_function_MMM(&sgemm_2, "MMM C opts", 1);
  add_function_MMM(&sgemm_intel, "MMM Intel", 1);
}

#endif
