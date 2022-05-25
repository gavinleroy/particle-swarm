#pragma once

typedef int (*fit_surrogate_fun_t)(struct pso_data_constant_inertia *);

void add_function_FIT_SURROGATE(fit_surrogate_fun_t f, char const *name,
                                int flop);
void register_functions_FIT_SURROGATE();
int perf_test_fit_surrogate(struct pso_data_constant_inertia *);
