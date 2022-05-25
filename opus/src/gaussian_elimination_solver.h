#pragma once

int gaussian_elimination_solve(int N, double *Ab, double *x);

#ifdef TEST_PERF

void register_functions_GE_SOLVE();

#endif
