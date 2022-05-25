#include "lu_solve.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "helpers.h"
// #include "mkl.h"

#ifdef TEST_PERF
#include "perf_testers/perf_lu_solve.h"
#endif

/**
 * === State of Affairs ===
 *
 * The bottleneck for this algorithm is dgemm (inappropriately named sgemm
 * here). This is especially the case because the leading dimension of the
 * matrices will most likely be >>512, meaning that a single row occupies at
 * least 1 page in memory (4096 Bytes).
 *
 * Optimization things Gavin plans on trying.
 *
 * - Allocating a static contiguous
 * amount of memory in this file, and copying each matrix into that first,
 * then performing the multiplication, then storing again. This will reduce
 * number pages loaded when accessing different rows of the matrix. Even
 * if we choose an outer block size to optimize for the cache, it doesn't
 * matter if each row is stored on a different page in memory.
 *
 * TODO add more ideas here
 */

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
// Utilities

static void swapd(double *a, int i, int j)
{
  double t = a[i];
  a[i] = a[j];
  a[j] = t;
}

// -----------------
// LAPACK routines
//
// NOTE these functions have no declaration. However, a
// corresponding <name>_XXX  definition is made for each
// iteration XXX. This way performance testing is consistent.
// Forward declaration are here merely to have useful doc
// comments in one place.

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
static int isamax_2(int N, double *A, int stride)
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
      t_7  //
      ;

  int i, ix, iy;

  if (incx == 1 && incy == 1)
  {
    for (i = 0; i < N - 7; i += 8)
    {
      t_0 = X[i + 0];
      t_1 = X[i + 1];
      t_2 = X[i + 2];
      t_3 = X[i + 3];
      t_4 = X[i + 4];
      t_5 = X[i + 5];
      t_6 = X[i + 6];
      t_7 = X[i + 7];

      X[i + 0] = Y[i + 0];
      X[i + 1] = Y[i + 1];
      X[i + 2] = Y[i + 2];
      X[i + 3] = Y[i + 3];
      X[i + 4] = Y[i + 4];
      X[i + 5] = Y[i + 5];
      X[i + 6] = Y[i + 6];
      X[i + 7] = Y[i + 7];

      Y[i + 0] = t_0;
      Y[i + 1] = t_1;
      Y[i + 2] = t_2;
      Y[i + 3] = t_3;
      Y[i + 4] = t_4;
      Y[i + 5] = t_5;
      Y[i + 6] = t_6;
      Y[i + 7] = t_7;
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
          tmp = MIX(A, LDA, i, k + 0);
          MIX(A, LDA, i, k + 0) = MIX(A, LDA, p_i, k + 0);
          MIX(A, LDA, p_i, k + 0) = tmp;
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

/**
 * TODO continue optimizing DGEMM I ran out of time :(
 */
static void sgemm_2(int M, int N, int K, double alpha, double *A, int LDA,
                    double *B, int LDB, double beta, double *C, int LDC)
{

  // NOTE as written below, we specialize to alpha = -1 beta = 1
  assert(APPROX_EQUAL(beta, ONE));
  assert(APPROX_EQUAL(alpha, -ONE));

  // ---------------
  // https://en.wikichip.org/wiki/intel/microarchitectures/skylake_(client)
  //
  // NOTE Skylake has a 32000 B L1 cache
  // With cache lines of 64 B
#define CACHE_BLOCK 16
  const int NB = CACHE_BLOCK;
  const int DBL_LMT = 512;

  double AL[CACHE_BLOCK * CACHE_BLOCK];
  double BL[CACHE_BLOCK * CACHE_BLOCK];
  double CL[CACHE_BLOCK * CACHE_BLOCK];

  // Skylake has 32 FP registers
  // NOTE choose MU + NU + MU * NU <= #_registers
  int MU = 4;
  int NU = 4;
  int KU = 1;

  assert(M % NB == 0);
  assert(N % NB == 0);
  assert(K % NB == 0);

  int i, j, k,      //
      ii, jj, kk,   //
      iii, jjj, kkk //
      ;

  // If the leading dimension is greater than a page,
  // loading and storing into contiguous memory could
  // improve performance.
  if (DBL_LMT < LDA || DBL_LMT < LDB || DBL_LMT < LDC)
  {

    for (j = 0; j < N - NB + 1; j += NB)
      for (i = 0; i < M - NB + 1; i += NB)
        for (k = 0; k < K - NB + 1; k += NB)
        {

          // Load matrix into contiguous memory

          for (ii = 0; ii < NB; ++ii)
            for (jj = 0; jj < NB; ++jj)
            {
              MIX(AL, NB, ii, jj) = MIX(A, LDA, i + ii, k + jj);
              MIX(BL, NB, ii, jj) = MIX(B, LDB, k + ii, j + jj);
              MIX(CL, NB, ii, jj) = MIX(C, LDC, i + ii, j + jj);
            }

          // ----

          // Mini MMM (cache blocking)
          for (jj = 0; jj < NB; jj += NU)
            for (ii = 0; ii < NB; ii += MU)
              for (kk = 0; kk < NB; kk += KU)

                // Micro MMM (register blocking)
                for (kkk = kk; kkk < kk + KU; ++kkk)
                  for (iii = ii; iii < ii + MU; ++iii)
                    for (jjj = jj; jjj < jj + NU; ++jjj)

                      MIX(CL, NB, iii, jjj) =
                          MIX(CL, NB, iii, jjj) -
                          MIX(BL, NB, kkk, jjj) * MIX(AL, NB, iii, kkk);

          // ----

          // Store the local matrices back into memory.

          for (ii = 0; ii < NB; ++ii)
            for (jj = 0; jj < NB; ++jj)
            {
              MIX(A, LDA, i + ii, k + jj) = MIX(AL, NB, ii, jj);
              MIX(B, LDB, k + ii, j + jj) = MIX(BL, NB, ii, jj);
              MIX(C, LDC, i + ii, j + jj) = MIX(CL, NB, ii, jj);
            }
        }

  } // If outside page bounds

  else
  {

    for (j = 0; j < N; j += NB)
      for (i = 0; i < M; i += NB)
        for (k = 0; k < K; k += NB)

          // Mini MMM (cache blocking)
          for (jj = j; jj < j + NB; jj += NU)
            for (ii = i; ii < i + NB; ii += MU)
              for (kk = k; kk < k + NB; kk += KU)

                // Micro MMM (register blocking)
                for (kkk = kk; kkk < kk + KU; ++kkk)
                  for (iii = ii; iii < ii + MU; ++iii)
                    for (jjj = jj; jjj < jj + NU; ++jjj)

                      MIX(C, LDC, iii, jjj) =
                          MIX(C, LDC, iii, jjj) -
                          MIX(B, LDB, kkk, jjj) * MIX(A, LDA, iii, kkk);

  } // Else within page bounds

  // Cleanup cases. This should also be optimized.
  for (; j < N; ++j)
    for (; i < N; ++i)
      for (; k < N; ++k)
        MIX(C, LDC, i, j) =
            MIX(C, LDC, i, j) - MIX(B, LDB, k, j) * MIX(A, LDA, i, k);
}

int sgetf2_2(int M, int N, double *A, int LDA, int *ipiv)
{

  int i, j, k, p_i;

  double   //
      p_v, //
      m_0  //
      ;

  // NOTE I have increased this to iterate *until* N however you could
  // stop at < N - 1 if you cover the bounds case.
  for (i = 0; i < MIN(M, N); ++i)
  {

    p_i = i + isamax_2(M - i, &MIX(A, LDA, i, i), LDA);
    p_v = MIX(A, LDA, p_i, i);

    if (APPROX_EQUAL(p_v, 0.))
    {
      fprintf(stderr, "ERROR: LU Solve singular matrix\n");
      return -1;
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
          MIX(A, LDA, j, k + 0) =
              m_0 * MIX(A, LDA, i, k + 0) + MIX(A, LDA, j, k + 0);
          MIX(A, LDA, j, k + 1) =
              m_0 * MIX(A, LDA, i, k + 1) + MIX(A, LDA, j, k + 1);
          MIX(A, LDA, j, k + 2) =
              m_0 * MIX(A, LDA, i, k + 2) + MIX(A, LDA, j, k + 2);
          MIX(A, LDA, j, k + 3) =
              m_0 * MIX(A, LDA, i, k + 3) + MIX(A, LDA, j, k + 3);

          MIX(A, LDA, j, k + 4) =
              m_0 * MIX(A, LDA, i, k + 4) + MIX(A, LDA, j, k + 4);
          MIX(A, LDA, j, k + 5) =
              m_0 * MIX(A, LDA, i, k + 5) + MIX(A, LDA, j, k + 5);
          MIX(A, LDA, j, k + 6) =
              m_0 * MIX(A, LDA, i, k + 6) + MIX(A, LDA, j, k + 6);
          MIX(A, LDA, j, k + 7) =
              m_0 * MIX(A, LDA, i, k + 7) + MIX(A, LDA, j, k + 7);
        }

        for (; k < N; ++k)
        {
          MIX(A, LDA, j, k) = m_0 * MIX(A, LDA, i, k) + MIX(A, LDA, j, k);
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

  // FIXME XXX we can play with NB
  const int NB = 256, //
      M = N,          //
      LDA = N,        //
      MIN_MN = N      //
      ;

  // BLocked factor A into [L \ U]

  if (NB <= 1 || NB >= MIN_MN)
  {
    // Use unblocked code
    retcode = sgetf2_2(M, N, A, LDA, ipiv);
    if (retcode != 0)
      return retcode;
  }
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

        if (ib + IB < M) // NOTE for a square matrix this will always be true.
        {
          // XXX for most blocks, the height of the multiplication 'M - ib - IB'
          // is going to be much much larger than IB.
          //
          // Update trailing submatrix
          sgemm_2(
              // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              M - ib - IB, N - ib - IB, IB, -ONE, //
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
  retcode = sgetrs_2(N, A, ipiv, b);

  return retcode;
}

#ifdef TEST_PERF

// NOTE I bet we can put these in the templated perf framework and remove the
// weirdness of compiling seperately and linking.
void register_functions_LU_SOLVE()
{
  add_function_LU_SOLVE(&lu_solve_0, "LU Solve Base", 1);
  add_function_LU_SOLVE(&lu_solve_1, "LU Solve Recursive", 1);
  add_function_LU_SOLVE(&lu_solve_2, "LU Solve Unrolled", 1);
}

#endif
