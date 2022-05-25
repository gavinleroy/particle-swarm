#pragma once

#include "../pso.h"

#ifndef FIT_SURROGATE_VERSION
#define FIT_SURROGATE_VERSION 4
#endif

int fit_surrogate(struct pso_data_constant_inertia *pso);
int prealloc_fit_surrogate(size_t max_n_phi, size_t n_P);

int fit_surrogate_0(struct pso_data_constant_inertia *pso);
int fit_surrogate_1(struct pso_data_constant_inertia *pso);
int fit_surrogate_2(struct pso_data_constant_inertia *pso);
int fit_surrogate_3(struct pso_data_constant_inertia *pso);
int fit_surrogate_4(struct pso_data_constant_inertia *pso);

int prealloc_fit_surrogate_0(size_t max_n_phi, size_t n_P);
int prealloc_fit_surrogate_1(size_t max_n_phi, size_t n_P);
int prealloc_fit_surrogate_2(size_t max_n_phi, size_t n_P);
int prealloc_fit_surrogate_3(size_t max_n_phi, size_t n_P);
int prealloc_fit_surrogate_4(size_t max_n_phi, size_t n_P);
