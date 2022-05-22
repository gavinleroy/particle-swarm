#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "lu_solve.h"
// #include "perf_lu_solve.h"

int main()
{

  double A[64] = {
      0.166402,    0.222061,  0.118731,  0.921609,  0.0123858, 0.184247,
      0.204479,    0.9318,    0.325253,  0.271533,  0.468016,  0.854823,
      0.592775,    0.57451,   0.960148,  0.762938,  0.711894,  0.371753,
      0.314479,    0.710879,  0.65966,   0.0394215, 0.140754,  0.815794,
      0.000755813, 0.0489883, 0.438471,  0.995461,  0.0889002, 0.186734,
      0.53057,     0.900928,  0.797166,  0.136262,  0.909005,  0.900479,
      0.551321,    0.767061,  0.15897,   0.0586111, 0.584989,  0.836514,
      0.0727873,   0.687834,  0.0800794, 0.249618,  0.0793851, 0.964427,
      0.404865,    0.528879,  0.786431,  0.749641,  0.921302,  0.93087,
      0.74637,     0.222083,  0.0549027, 0.187067,  0.940697,  0.833942,
      0.869757,    0.786435,  0.433628,  0.24035};

  double b[8] = {0.3614445872344271,   0.9188848971981615, 0.06485553257153076,
                 0.4658661689221989,   0.2797787272746787, 0.20231195354569353,
                 0.044427089869083725, 0.5169398670162666};

  double expected_x[8] = {0.612187230540563,   -3.0112661319020124,
                          1.4219085852331368,  -4.382256470614417,
                          -0.5213769648887188, 3.9237104595881513,
                          -0.8441151106088531, 4.565676075132954};

  double *Ab = (double *)aligned_alloc(32, 9 * 8 * sizeof(double));
  double *x = (double *)aligned_alloc(32, 8 * sizeof(double));

  for (int i = 0; i < 8; i++)
    for (int j = 0; j < 8; j++)
      Ab[i * 9 + j] = A[i * 8 + j];

  for (int i = 0; i < 8; i++)
    Ab[i * 9 + 8] = b[i];

  int ret = lu_solve(8, A, b, x); // gaussian_elimination_solve(8, Ab, x);

  assert(ret == 0);

  printf("-ACT- -EXP-\n");
  for (int i = 0; i < 8; i++)
    printf("%f %f\n", x[i], expected_x[i]);

  // perf_test_ge_solve(8, Ab, x);

  return 0;
}
