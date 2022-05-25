#pragma once

typedef int (*ge_solve_fun_t)(int N, double *Ab, double *x);

void add_function_GE_SOLVE(ge_solve_fun_t f, char *name, int flop);
void register_functions_GE_SOLVE();
int perf_test_ge_solve(int N, double *Ab, double *x);
