#pragma once

#include "pso.h"
#include "timing_params.h"

enum profiled_function_t {run_pso_f, plu_factorization_f, plu_solve_f, fit_surrogate_f};

double perf_counter(enum profiled_function_t profiled_function, profiling_hyperparameters * hyperparams);