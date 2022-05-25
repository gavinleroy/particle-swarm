#pragma once

#include <stdlib.h>

/** @brief Solve linear systems using LU factorization method.
 *
 * @param A: The segment of memory representing a row-major square
 *           matrix layed out sequentially in memory.
 * @param b: Column vector b in Ax=b
 * @param ipiv: Pivot indices, row i was interchanged with row ipiv[i].
 * @param x: Output buffer for column vector x in Ax=b
 * @param N: The length of one side.
 * @return: 0 on success. <0 for matrix singularity.
 */
int lu_solve(int N, double *A, int *ipiv, double *b);

#ifdef TEST_PERF

void register_functions_LU_SOLVE();

#endif
