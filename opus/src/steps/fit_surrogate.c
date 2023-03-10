#include "fit_surrogate.h"

#include <immintrin.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "../helpers.h"
#include "../pso.h"
#include "linear_system_solver.h"

#include "../my_papi.h"

int fit_surrogate_6_GE(struct pso_data_constant_inertia *pso);
int fit_surrogate_6_LU(struct pso_data_constant_inertia *pso);
int fit_surrogate_6_BLOCK_TRI(struct pso_data_constant_inertia *pso);

int prealloc_fit_surrogate_6_GE(size_t max_n_phi, size_t n_P);
int prealloc_fit_surrogate_6_LU(size_t max_n_phi, size_t n_P);
int prealloc_fit_surrogate_6_BLOCK_TRI(size_t max_n_phi, size_t n_P);

int fit_surrogate(struct pso_data_constant_inertia *pso)
{
  PAPI_START("fit_surrogate");
  int ret = FIT_SURROGATE_VERSION(pso);
  PAPI_STOP("fit_surrogate");
  return ret;
}

int prealloc_fit_surrogate(size_t max_n_phi, size_t n_P)
{
  return FIT_SURROGATE_PREALLOC_VERSION(max_n_phi, n_P);
}

#define DEBUG_SURROGATE 0

// TODO: include past_refinement_points in phi !!!

// is either [A | b] for GE and BLOCK_TRI or [A] for LU
static double *fit_surrogate_Ab;
static double *fit_surrogate_P;
// if using LU
static double *fit_surrogate_b;

size_t fit_surrogate_max_N_phi;
double *fit_surrogate_phi_cache;

int prealloc_fit_surrogate_0(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and n x 1 for b
  size_t Ab_size = max_n_A * max_n_A + max_n_A;

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  return 0;
}

int fit_surrogate_0(struct pso_data_constant_inertia *pso)
{

  // the size of phi is the total number of _distinct_ points where
  // f has been evaluated
  // currently : n = x_distinct_s
  size_t n_phi = pso->x_distinct_s;

  // the size of P is n x d+1
  size_t n_P = pso->dimensions + 1;

  // the size of the matrix in the linear system is n+d+1
  size_t n_A = n_phi + n_P;

  double *Ab = fit_surrogate_Ab;

  /********
   * Prepare left hand side A
   ********/

  // phi_p,q = || u_p - u_q ||^3

  for (size_t k1 = 0; k1 < pso->x_distinct_s; k1++)
  {
    double *u_p = pso->x_distinct + k1 * pso->dimensions;

    for (size_t k2 = 0; k2 < pso->x_distinct_s; k2++)
    {
      double *u_q = pso->x_distinct + k2 * pso->dimensions;
      double d2 = dist2(pso->dimensions, u_p, u_q);
      double d = sqrt(d2);
      double d3 = d2 * d;
      Ab[k1 * (n_A + 1) + k2] = d3;
    }
  }

  // P and tP are blocks in A
  // P_{i,j} := A_{i,n_phi + j} = A[i * n_A + n_phi + j]
  // tP_{i,j} := A_{n_phi + i, j} = A[(n_phi + i) * n_A + j]

  for (size_t k = 0; k < pso->x_distinct_s; k++)
  {
    double *u = pso->x_distinct + k * pso->dimensions;

    // P(p,0) = 1;
    Ab[k * (n_A + 1) + n_phi + 0] = 1;
    // tP(0,p) = 1;
    Ab[(n_phi + 0) * (n_A + 1) + k] = 1;

    for (int j = 0; j < pso->dimensions; j++)
    {
      // P(p,1+j) = u[j];
      Ab[k * (n_A + 1) + n_phi + j + 1] = u[j];
      // tP(1+j,p) = u[j];
      Ab[(n_phi + 1 + j) * (n_A + 1) + k] = u[j];
    }
  }

  // lower right block is zeros
  for (size_t i = n_phi; i < n_A; i++)
  {
    for (size_t j = n_phi; j < n_A; j++)
    {
      Ab[i * (n_A + 1) + j] = 0;
    }
  }

  /********
   * Prepare right hand side b
   ********/
  for (size_t k = 0; k < n_phi; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = pso->x_distinct_eval[k];
  }

  for (size_t k = n_phi; k < n_A; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = 0;
  }

#if DEBUG_SURROGATE
  print_rect_matrixd(Ab, n_A, n_A + 1, "Ab");
#endif

  PAPI_START("system_solver");
  int ret = gaussian_elimination_solve(n_A, Ab, pso->lambda_p);
  PAPI_STOP("system_solver");

  if (ret < 0)
  {
    return -1;
  }

#if DEBUG_SURROGATE
  print_vectord(pso->lambda_p, n_A, "x");
#endif

  return 0;
}

int prealloc_fit_surrogate_1(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and n x 1 for b
  size_t Ab_size = max_n_A * max_n_A + max_n_A;

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  return 0;
}

/*
 * replace all occurences of pso-> by scalar replacement
 */

int fit_surrogate_1(struct pso_data_constant_inertia *pso)
{
  // TODO: include past_refinement_points in phi !!!

  size_t dimensions = pso->dimensions;
  size_t popsize = pso->population_size;
  size_t t = pso->time;

  // the size of phi is the total number of _distinct_ points where
  // f has been evaluated
  // currently : n = x_distinct_s
  size_t n_phi = pso->x_distinct_s;
  double *x_distincts = pso->x_distinct;
  double *fxd = pso->x_distinct_eval;

  // the size of P is n x d+1
  size_t n_P = dimensions + 1;

  // the size of the matrix in the linear system is n+d+1
  size_t n_A = n_phi + n_P;

  double *Ab = fit_surrogate_Ab;

  /********
   * Prepare left hand side A
   ********/

  // phi_p,q = || u_p - u_q ||^3

  for (size_t k1 = 0; k1 < n_phi; k1++)
  {
    double *u_p = x_distincts + k1 * dimensions;

    for (size_t k2 = 0; k2 < n_phi; k2++)
    {
      double *u_q = x_distincts + k2 * dimensions;
      double d2 = dist2(dimensions, u_p, u_q);
      double d = sqrt(d2);
      double d3 = d2 * d;
      Ab[k1 * (n_A + 1) + k2] = d3;
    }
  }

  // P and tP are blocks in A
  // P_{i,j} := A_{i,n_phi + j} = A[i * n_A + n_phi + j]
  // tP_{i,j} := A_{n_phi + i, j} = A[(n_phi + i) * n_A + j]

  for (size_t k = 0; k < n_phi; k++)
  {
    double *u = x_distincts + k * dimensions;

    // P(p,0) = 1;
    Ab[k * (n_A + 1) + n_phi + 0] = 1;
    // tP(0,p) = 1;
    Ab[(n_phi + 0) * (n_A + 1) + k] = 1;

    for (int j = 0; j < dimensions; j++)
    {
      // P(p,1+j) = u[j];
      Ab[k * (n_A + 1) + n_phi + j + 1] = u[j];
      // tP(1+j,p) = u[j];
      Ab[(n_phi + 1 + j) * (n_A + 1) + k] = u[j];
    }
  }

  // lower right block is zeros
  for (size_t i = n_phi; i < n_A; i++)
  {
    for (size_t j = n_phi; j < n_A; j++)
    {
      Ab[i * (n_A + 1) + j] = 0;
    }
  }

  /********
   * Prepare right hand side b
   ********/
  for (size_t k = 0; k < n_phi; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = fxd[k];
  }

  for (size_t k = n_phi; k < n_A; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = 0;
  }

#if DEBUG_SURROGATE
  print_rect_matrixd(Ab, n_A, n_A + 1, "Ab");
#endif

  if (gaussian_elimination_solve(n_A, Ab, pso->lambda_p) < 0)
  {
    return -1;
  }

#if DEBUG_SURROGATE
  print_vectord(pso->lambda_p, n_A, "x");
#endif

  return 0;
}

/*
 * Cache distances computation: cache matrix phi
 */

int prealloc_fit_surrogate_2(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and n x 1 for b
  size_t Ab_size = max_n_A * max_n_A + max_n_A;

  fit_surrogate_max_N_phi = max_n_phi;
  fit_surrogate_phi_cache = malloc(max_n_phi * max_n_phi * sizeof(double));

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  return 0;
}

int fit_surrogate_2(struct pso_data_constant_inertia *pso)
{
  // TODO: include past_refinement_points in phi !!!

  size_t dimensions = pso->dimensions;
  size_t popsize = pso->population_size;
  size_t t = pso->time;

  // the size of phi is the total number of _distinct_ points where
  // f has been evaluated
  // currently : n = x_distinct_s
  size_t n_phi = pso->x_distinct_s;
  double *x_distincts = pso->x_distinct;
  double *fxd = pso->x_distinct_eval;

  // the size of P is n x d+1
  size_t n_P = dimensions + 1;

  // the size of the matrix in the linear system is n+d+1
  size_t n_A = n_phi + n_P;

  double *Ab = fit_surrogate_Ab;

  size_t max_N_phi = fit_surrogate_max_N_phi;
  double *phi_cache = fit_surrogate_phi_cache;
  /********
   * Prepare left hand side A
   ********/

  // phi_p,q = || u_p - u_q ||^3

  // j = 0..time are precomputed,
  // need to compute for time+1

  size_t id_new_points = pso->x_distinct_idx_of_last_batch;
#if DEBUG_SURROGATE
  fprintf(stderr, "Compute distances of %zu new points\n",
          n_phi - id_new_points);
#endif

  for (int k1 = id_new_points; k1 < n_phi; k1++)
  {
    double *u_p = x_distincts + k1 * dimensions;
    for (size_t k2 = 0; k2 < id_new_points; k2++)
    {
      double *u_q = x_distincts + k2 * dimensions;
      double d2 = dist2(dimensions, u_p, u_q);
      double d = sqrt(d2);
      double d3 = d2 * d;
      phi_cache[k1 * max_N_phi + k2] = d3;
      phi_cache[k2 * max_N_phi + k1] = d3;
    }

    // d(p_k1, p_k1) = 0
    phi_cache[k1 * max_N_phi + k1] = 0;

    // Avoid double computation for id_new_points < k < n_phi
    // ie avoid pairs k2 < k1
    for (size_t k2 = k1 + 1; k2 < n_phi; k2++)
    {
      double *u_q = x_distincts + k2 * dimensions;
      double d2 = dist2(dimensions, u_p, u_q);
      double d = sqrt(d2);
      double d3 = d2 * d;
      phi_cache[k1 * max_N_phi + k2] = d3;
      phi_cache[k2 * max_N_phi + k1] = d3;
    }
  }
  // start new batch
  pso->x_distinct_idx_of_last_batch = n_phi;

  // Copy the distances from phi_cache to the phi block in A
  // TODO: use symmetry when building matrix? or keep seq accesses?
  for (size_t k1 = 0; k1 < n_phi; k1++)
  {
    double *u = x_distincts + k1 * dimensions;

    for (size_t k2 = 0; k2 < n_phi; k2++)
    {
      Ab[k1 * (n_A + 1) + k2] = phi_cache[k1 * max_N_phi + k2];
    }
  }

  // P and tP are blocks in A
  // P_{i,j} := A_{i,n_phi + j} = A[i * n_A + n_phi + j]
  // tP_{i,j} := A_{n_phi + i, j} = A[(n_phi + i) * n_A + j]

  for (size_t k = 0; k < n_phi; k++)
  {
    double *u = x_distincts + k * dimensions;

    // P(p,0) = 1;
    Ab[k * (n_A + 1) + n_phi + 0] = 1;
    // tP(0,p) = 1;
    Ab[(n_phi + 0) * (n_A + 1) + k] = 1;

    for (int j = 0; j < dimensions; j++)
    {
      // P(p,1+j) = u[j];
      Ab[k * (n_A + 1) + n_phi + j + 1] = u[j];
      // tP(1+j,p) = u[j];
      Ab[(n_phi + 1 + j) * (n_A + 1) + k] = u[j];
    }
  }

  // lower right block is zeros
  for (size_t i = n_phi; i < n_A; i++)
  {
    for (size_t j = n_phi; j < n_A; j++)
    {
      Ab[i * (n_A + 1) + j] = 0;
    }
  }

  /********
   * Prepare right hand side b
   ********/
  for (size_t k = 0; k < n_phi; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = fxd[k];
  }

  for (size_t k = n_phi; k < n_A; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = 0;
  }

#if DEBUG_SURROGATE
  print_rect_matrixd(Ab, n_A, n_A + 1, "Ab");
#endif

  if (gaussian_elimination_solve(n_A, Ab, pso->lambda_p) < 0)
  {
    return -1;
  }

#if DEBUG_SURROGATE
  print_vectord(pso->lambda_p, n_A, "x");
#endif

  return 0;
}

/*
 * Cache distances computation: same as before, but use memcpy
 */

int prealloc_fit_surrogate_3(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and n x 1 for b
  size_t Ab_size = max_n_A * max_n_A + max_n_A;

  fit_surrogate_max_N_phi = max_n_phi;
  fit_surrogate_phi_cache = malloc(max_n_phi * max_n_phi * sizeof(double));

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  return 0;
}

int fit_surrogate_3(struct pso_data_constant_inertia *pso)
{
  // TODO: include past_refinement_points in phi !!!

  size_t dimensions = pso->dimensions;
  size_t popsize = pso->population_size;
  size_t t = pso->time;

  // the size of phi is the total number of _distinct_ points where
  // f has been evaluated
  // currently : n = x_distinct_s
  size_t n_phi = pso->x_distinct_s;
  double *x_distincts = pso->x_distinct;
  double *fxd = pso->x_distinct_eval;

  // the size of P is n x d+1
  size_t n_P = dimensions + 1;

  // the size of the matrix in the linear system is n+d+1
  size_t n_A = n_phi + n_P;

  double *Ab = fit_surrogate_Ab;

  size_t max_N_phi = fit_surrogate_max_N_phi;
  double *phi_cache = fit_surrogate_phi_cache;
  /********
   * Prepare left hand side A
   ********/

  // phi_p,q = || u_p - u_q ||^3

  // j = 0..time are precomputed,
  // need to compute for time+1

  size_t id_new_points = pso->x_distinct_idx_of_last_batch;
#if DEBUG_SURROGATE
  fprintf(stderr, "Compute distances of %zu new points\n",
          n_phi - id_new_points);
#endif

  for (int k1 = id_new_points; k1 < n_phi; k1++)
  {
    double *u_p = x_distincts + k1 * dimensions;
    for (size_t k2 = 0; k2 < id_new_points; k2++)
    {
      double *u_q = x_distincts + k2 * dimensions;
      double d2 = dist2(dimensions, u_p, u_q);
      double d = sqrt(d2);
      double d3 = d2 * d;
      phi_cache[k1 * max_N_phi + k2] = d3;
      phi_cache[k2 * max_N_phi + k1] = d3;
    }

    // d(p_k1, p_k1) = 0
    phi_cache[k1 * max_N_phi + k1] = 0;

    // Avoid double computation for id_new_points < k < n_phi
    // ie avoid pairs k2 < k1
    for (size_t k2 = k1 + 1; k2 < n_phi; k2++)
    {
      double *u_q = x_distincts + k2 * dimensions;
      double d2 = dist2(dimensions, u_p, u_q);
      double d = sqrt(d2);
      double d3 = d2 * d;
      phi_cache[k1 * max_N_phi + k2] = d3;
      phi_cache[k2 * max_N_phi + k1] = d3;
    }
  }
  // start new batch
  pso->x_distinct_idx_of_last_batch = n_phi;

  // Copy the distances from phi_cache to the phi block in A
  for (size_t k1 = 0; k1 < n_phi; k1++)
  {
    double *phi_row = phi_cache + k1 * max_N_phi;
    double *Ab_row = Ab + k1 * (n_A + 1);
    memcpy(Ab_row, phi_row, n_phi * sizeof(double));
  }

  // P and tP are blocks in A
  // P_{i,j} := A_{i,n_phi + j} = A[i * n_A + n_phi + j]
  // tP_{i,j} := A_{n_phi + i, j} = A[(n_phi + i) * n_A + j]

  for (size_t k = 0; k < n_phi; k++)
  {
    double *u = x_distincts + k * dimensions;

    // P(p,0) = 1;
    Ab[k * (n_A + 1) + n_phi + 0] = 1;
    // tP(0,p) = 1;
    Ab[(n_phi + 0) * (n_A + 1) + k] = 1;

    for (int j = 0; j < dimensions; j++)
    {
      // P(p,1+j) = u[j];
      Ab[k * (n_A + 1) + n_phi + j + 1] = u[j];
      // tP(1+j,p) = u[j];
      Ab[(n_phi + 1 + j) * (n_A + 1) + k] = u[j];
    }
  }

  // lower right block is zeros
  for (size_t i = n_phi; i < n_A; i++)
  {
    for (size_t j = n_phi; j < n_A; j++)
    {
      Ab[i * (n_A + 1) + j] = 0;
    }
  }

  /********
   * Prepare right hand side b
   ********/
  for (size_t k = 0; k < n_phi; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = fxd[k];
  }

  for (size_t k = n_phi; k < n_A; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = 0;
  }

#if DEBUG_SURROGATE
  print_rect_matrixd(Ab, n_A, n_A + 1, "Ab");
#endif

  if (gaussian_elimination_solve(n_A, Ab, pso->lambda_p) < 0)
  {
    return -1;
  }

#if DEBUG_SURROGATE
  print_vectord(pso->lambda_p, n_A, "x");
#endif

  return 0;
}

/*
 * Cache distances + early exit if no new elements
 */

int prealloc_fit_surrogate_4(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and n x 1 for b
  size_t Ab_size = max_n_A * max_n_A + max_n_A;

  fit_surrogate_max_N_phi = max_n_phi;
  fit_surrogate_phi_cache = malloc(max_n_phi * max_n_phi * sizeof(double));

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  return 0;
}

int fit_surrogate_4(struct pso_data_constant_inertia *pso)
{
  // TODO: include past_refinement_points in phi !!!

  size_t dimensions = pso->dimensions;
  size_t popsize = pso->population_size;
  size_t t = pso->time;

  // the size of phi is the total number of _distinct_ points where
  // f has been evaluated
  // currently : n = x_distinct_s
  size_t n_phi = pso->x_distinct_s;
  double *x_distincts = pso->x_distinct;
  double *fxd = pso->x_distinct_eval;

  // the size of P is n x d+1
  size_t n_P = dimensions + 1;

  // the size of the matrix in the linear system is n+d+1
  size_t n_A = n_phi + n_P;

  double *Ab = fit_surrogate_Ab;

  size_t max_N_phi = fit_surrogate_max_N_phi;
  double *phi_cache = fit_surrogate_phi_cache;
  /********
   * Prepare left hand side A
   ********/

  // phi_p,q = || u_p - u_q ||^3

  // j = 0..time are precomputed,
  // need to compute for time+1

  size_t id_new_points = pso->x_distinct_idx_of_last_batch;

  if (id_new_points == n_phi)
  {
    // There are no new points ! The surrogate is already fit !
#if DEBUG_SURROGATE
    printf("Skip fit_surrogate: no new evaluation position!\n");
#endif
    return 0;
  }

#if DEBUG_SURROGATE
  printf("Compute distances of %zu new points\n", n_phi - id_new_points);
#endif

  for (int k1 = id_new_points; k1 < n_phi; k1++)
  {
    double *u_p = x_distincts + k1 * dimensions;
    for (size_t k2 = 0; k2 < id_new_points; k2++)
    {
      double *u_q = x_distincts + k2 * dimensions;
      double d2 = dist2(dimensions, u_p, u_q);
      double d = sqrt(d2);
      double d3 = d2 * d;
      phi_cache[k1 * max_N_phi + k2] = d3;
      phi_cache[k2 * max_N_phi + k1] = d3;
    }

    // d(p_k1, p_k1) = 0
    phi_cache[k1 * max_N_phi + k1] = 0;

    // Avoid double computation for id_new_points < k < n_phi
    // ie avoid pairs k2 < k1
    for (size_t k2 = k1 + 1; k2 < n_phi; k2++)
    {
      double *u_q = x_distincts + k2 * dimensions;
      double d2 = dist2(dimensions, u_p, u_q);
      double d = sqrt(d2);
      double d3 = d2 * d;
      phi_cache[k1 * max_N_phi + k2] = d3;
      phi_cache[k2 * max_N_phi + k1] = d3;
    }
  }

  // start new batch
  pso->x_distinct_idx_of_last_batch = n_phi;

  // Copy the distances from phi_cache to the phi block in A
  for (size_t k1 = 0; k1 < n_phi; k1++)
  {
    double *phi_row = phi_cache + k1 * max_N_phi;
    double *Ab_row = Ab + k1 * (n_A + 1);
    memcpy(Ab_row, phi_row, n_phi * sizeof(double));
  }

  // P and tP are blocks in A
  // P_{i,j} := A_{i,n_phi + j} = A[i * n_A + n_phi + j]
  // tP_{i,j} := A_{n_phi + i, j} = A[(n_phi + i) * n_A + j]

  for (size_t k = 0; k < n_phi; k++)
  {
    double *u = x_distincts + k * dimensions;

    // P(p,0) = 1;
    Ab[k * (n_A + 1) + n_phi + 0] = 1;
    // tP(0,p) = 1;
    Ab[(n_phi + 0) * (n_A + 1) + k] = 1;

    for (int j = 0; j < dimensions; j++)
    {
      // P(p,1+j) = u[j];
      Ab[k * (n_A + 1) + n_phi + j + 1] = u[j];
      // tP(1+j,p) = u[j];
      Ab[(n_phi + 1 + j) * (n_A + 1) + k] = u[j];
    }
  }

  // lower right block is zeros
  for (size_t i = n_phi; i < n_A; i++)
  {
    for (size_t j = n_phi; j < n_A; j++)
    {
      Ab[i * (n_A + 1) + j] = 0;
    }
  }

  /********
   * Prepare right hand side b
   ********/
  for (size_t k = 0; k < n_phi; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = fxd[k];
  }

  for (size_t k = n_phi; k < n_A; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = 0;
  }

#if DEBUG_SURROGATE
  print_rect_matrixd(Ab, n_A, n_A + 1, "Ab");
#endif

  if (gaussian_elimination_solve(n_A, Ab, pso->lambda_p) < 0)
  {
    return -1;
  }

#if DEBUG_SURROGATE
  print_vectord(pso->lambda_p, n_A, "x");
#endif

  return 0;
}

/*
 * Make cache more dense:
 * phi(0,1) phi(0,2) phi(1,2) phi(0,3) phi(1,3) phi(2,3) phi(0,4) ...
 *
 * phi(i,j) with i<j stored at (j * (j-1) + i)
 */
int prealloc_fit_surrogate_5(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and n x 1 for b
  size_t Ab_size = max_n_A * (max_n_A + 1);
  size_t phi_cache_size = max_n_phi * (max_n_phi - 1) / 2;

  fit_surrogate_max_N_phi = max_n_phi;

  fit_surrogate_phi_cache = malloc(phi_cache_size * sizeof(double));

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  return 0;
}

int fit_surrogate_5(struct pso_data_constant_inertia *pso)
{
  // TODO: include past_refinement_points in phi !!!

  size_t dimensions = pso->dimensions;
  size_t popsize = pso->population_size;
  size_t t = pso->time;

  // the size of phi is the total number of _distinct_ points where
  // f has been evaluated
  // currently : n = x_distinct_s
  size_t n_phi = pso->x_distinct_s;
  double *x_distincts = pso->x_distinct;
  double *fxd = pso->x_distinct_eval;

  // the size of P is n x d+1
  size_t n_P = dimensions + 1;

  // the size of the matrix in the linear system is n+d+1
  size_t n_A = n_phi + n_P;
  size_t n_Ab = n_A + 1;

  double *Ab = fit_surrogate_Ab;

  size_t max_N_phi = fit_surrogate_max_N_phi;
  double *phi_cache = fit_surrogate_phi_cache;
  /********
   * Prepare left hand side A
   ********/

  // phi_p,q = || u_p - u_q ||^3

  // j = 0..time are precomputed,
  // need to compute for time+1

  size_t id_new_points = pso->x_distinct_idx_of_last_batch;

  if (id_new_points == n_phi)
  {
    // There are no new points ! The surrogate is already fit !
#if DEBUG_SURROGATE
    printf("Skip fit_surrogate: no new evaluation position!\n");
#endif
    return 0;
  }

#if DEBUG_SURROGATE
  printf("Compute distances of %zu new points\n", n_phi - id_new_points);
#endif

  for (size_t j = id_new_points; j < n_phi; j++)
  {
    double *u_j = x_distincts + j * dimensions;
    double *d3_to_u_j_cached = phi_cache + j * (j - 1) / 2;
    for (size_t i = 0; i < j; i++)
    {
      double *u_i = x_distincts + i * dimensions;
      double d2 = dist2(dimensions, u_j, u_i);
      double d = sqrt(d2);
      double d3 = d2 * d;
      d3_to_u_j_cached[i] = d3;
    }
  }

  // start new batch
  pso->x_distinct_idx_of_last_batch = n_phi;

  // Copy the distances from phi_cache to the phi block in A
  for (size_t j = 0; j < n_phi; j++)
  {
    double *d3_to_u_j_cached = phi_cache + j * (j - 1) / 2;

    for (size_t i = 0; i < j; i++)
    {
      double phi_i_j = d3_to_u_j_cached[i];
      Ab[i * n_Ab + j] = phi_i_j;
      Ab[j * n_Ab + i] = phi_i_j;
    }
    Ab[j * n_Ab + j] = 0.;
  }

  // P and tP are blocks in A
  // P_{i,j} := A_{i,n_phi + j} = A[i * n_A + n_phi + j]
  // tP_{i,j} := A_{n_phi + i, j} = A[(n_phi + i) * n_A + j]

  for (size_t k = 0; k < n_phi; k++)
  {
    double *u = x_distincts + k * dimensions;

    // P(p,0) = 1;
    Ab[k * (n_A + 1) + n_phi + 0] = 1;
    // tP(0,p) = 1;
    Ab[(n_phi + 0) * (n_A + 1) + k] = 1;

    for (int j = 0; j < dimensions; j++)
    {
      // P(p,1+j) = u[j];
      Ab[k * (n_A + 1) + n_phi + j + 1] = u[j];
      // tP(1+j,p) = u[j];
      Ab[(n_phi + 1 + j) * (n_A + 1) + k] = u[j];
    }
  }

  // lower right block is zeros
  for (size_t i = n_phi; i < n_A; i++)
  {
    for (size_t j = n_phi; j < n_A; j++)
    {
      Ab[i * (n_A + 1) + j] = 0;
    }
  }

  /********
   * Prepare right hand side b
   ********/
  for (size_t k = 0; k < n_phi; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = fxd[k];
  }

  for (size_t k = n_phi; k < n_A; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = 0;
  }

#if DEBUG_SURROGATE
  print_rect_matrixd(Ab, n_A, n_A + 1, "Ab");
#endif

  if (gaussian_elimination_solve(n_A, Ab, pso->lambda_p) < 0)
  {
    return -1;
  }

#if DEBUG_SURROGATE
  print_vectord(pso->lambda_p, n_A, "x");
#endif

  return 0;
}

int prealloc_fit_surrogate_6(size_t max_n_phi, size_t n_P)
{
#if LINEAR_SYSTEM_SOLVER_USED == GE_SOLVER
  return prealloc_fit_surrogate_6_GE(max_n_phi, n_P);
#elif LINEAR_SYSTEM_SOLVER_USED == LU_SOLVER
  return prealloc_fit_surrogate_6_LU(max_n_phi, n_P);
#elif LINEAR_SYSTEM_SOLVER_USED == BLOCK_TRI_SOLVER
  return prealloc_fit_surrogate_6_BLOCK_TRI(max_n_phi, n_P);
#endif
}

int fit_surrogate_6(struct pso_data_constant_inertia *pso)
{
#if LINEAR_SYSTEM_SOLVER_USED == GE_SOLVER
  return fit_surrogate_6_GE(pso);
#elif LINEAR_SYSTEM_SOLVER_USED == LU_SOLVER
  return fit_surrogate_6_LU_blocked(pso);
#elif LINEAR_SYSTEM_SOLVER_USED == BLOCK_TRI_SOLVER
  return fit_surrogate_6_BLOCK_TRI(pso);
#endif
}

/*
 * Cooperation between prealloc_fit_surrogate_6 and check_distinct
 */

int prealloc_fit_surrogate_6_GE(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and n x 1 for b
  size_t Ab_size = max_n_A * (max_n_A + 1);
  size_t phi_cache_size = max_n_phi * (max_n_phi - 1) / 2;

  fit_surrogate_max_N_phi = max_n_phi;

  fit_surrogate_phi_cache = malloc(phi_cache_size * sizeof(double));

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  return 0;
}

int fit_surrogate_6_GE(struct pso_data_constant_inertia *pso)
{

  size_t dimensions = pso->dimensions;
  size_t popsize = pso->population_size;
  size_t t = pso->time;

  // the size of phi is the total number of _distinct_ points where
  // f has been evaluated
  // currently : n = x_distinct_s
  size_t n_phi = pso->x_distinct_s;
  double *x_distincts = pso->x_distinct;
  double *fxd = pso->x_distinct_eval;

  // the size of P is n x d+1
  size_t n_P = dimensions + 1;

  // the size of the matrix in the linear system is n+d+1
  size_t n_A = n_phi + n_P;
  size_t n_Ab = n_A + 1;

  double *Ab = fit_surrogate_Ab;

  size_t max_N_phi = fit_surrogate_max_N_phi;
  double *phi_cache = fit_surrogate_phi_cache;
  /********
   * Prepare left hand side A
   ********/

  // phi_p,q = || u_p - u_q ||^3
  // all of those are already computed in check_if_distinct!
  size_t prev_n_phi = pso->x_distinct_idx_of_last_batch;
  if (prev_n_phi == n_phi)
  {
    // There are no new points ! The surrogate is already fit !
#if DEBUG_SURROGATE
    printf("Skip fit_surrogate: no new evaluation position!\n");
#endif
    return 0;
  }
  else
  {
    pso->x_distinct_idx_of_last_batch = n_phi;
  }

  // Copy the distances from phi_cache to the phi block in A
  for (size_t j = 0; j < n_phi; j++)
  {
    double *d3_to_u_j_cached = phi_cache + j * (j - 1) / 2;

    for (size_t i = 0; i < j; i++)
    {
      double phi_i_j = d3_to_u_j_cached[i];
      Ab[i * n_Ab + j] = phi_i_j;
      Ab[j * n_Ab + i] = phi_i_j;
    }
    Ab[j * n_Ab + j] = 0.;
  }

  // P and tP are blocks in A
  // P_{i,j} := A_{i,n_phi + j} = A[i * n_A + n_phi + j]
  // tP_{i,j} := A_{n_phi + i, j} = A[(n_phi + i) * n_A + j]

  for (size_t k = 0; k < n_phi; k++)
  {
    double *u = x_distincts + k * dimensions;

    // P(p,0) = 1;
    Ab[k * (n_A + 1) + n_phi + 0] = 1;
    // tP(0,p) = 1;
    Ab[(n_phi + 0) * (n_A + 1) + k] = 1;

    for (int j = 0; j < dimensions; j++)
    {
      // P(p,1+j) = u[j];
      Ab[k * (n_A + 1) + n_phi + j + 1] = u[j];
      // tP(1+j,p) = u[j];
      Ab[(n_phi + 1 + j) * (n_A + 1) + k] = u[j];
    }
  }

  // lower right block is zeros
  for (size_t i = n_phi; i < n_A; i++)
  {
    for (size_t j = n_phi; j < n_A; j++)
    {
      Ab[i * (n_A + 1) + j] = 0;
    }
  }

  /********
   * Prepare right hand side b
   ********/
  for (size_t k = 0; k < n_phi; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = fxd[k];
  }

  for (size_t k = n_phi; k < n_A; k++)
  {
    // set b_k
    Ab[k * (n_A + 1) + n_A] = 0;
  }

#if DEBUG_SURROGATE
  print_rect_matrixd(Ab, n_A, n_A + 1, "Ab");
#endif

  if (gaussian_elimination_solve(n_A, Ab, pso->lambda_p) < 0)
  {
    return -1;
  }

#if DEBUG_SURROGATE
  print_vectord(pso->lambda_p, n_A, "x");
#endif

  return 0;
}

/*
 * Cooperation between prealloc_fit_surrogate_6_LU and check_distinct
 */

int prealloc_fit_surrogate_6_LU(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and b is stored separately
  size_t Ab_size = max_n_A * max_n_A;
  size_t b_size = max_n_A;
  size_t phi_cache_size = max_n_phi * (max_n_phi - 1) / 2;

  fit_surrogate_max_N_phi = max_n_phi;

  fit_surrogate_phi_cache = malloc(phi_cache_size * sizeof(double));

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  fit_surrogate_b = malloc(b_size * sizeof(double));

  lu_initialize_memory(n_P + max_n_phi);
  return 0;
}

int fit_surrogate_6_LU(struct pso_data_constant_inertia *pso)
{

  size_t dimensions = pso->dimensions;
  size_t popsize = pso->population_size;
  size_t t = pso->time;

  // the size of phi is the total number of _distinct_ points where
  // f has been evaluated
  // currently : n = x_distinct_s
  size_t n_phi = pso->x_distinct_s;
  double *x_distincts = pso->x_distinct;
  double *fxd = pso->x_distinct_eval;

  // the size of P is n x d+1
  size_t n_P = dimensions + 1;

  // the size of the matrix in the linear system is n+d+1
  size_t n_A = n_phi + n_P;
  size_t n_Ab = n_A + 1;

  double *A = fit_surrogate_Ab;

  double *b = fit_surrogate_b;

  size_t max_N_phi = fit_surrogate_max_N_phi;
  double *phi_cache = fit_surrogate_phi_cache;
  /********
   * Prepare left hand side A
   ********/

  // phi_p,q = || u_p - u_q ||^3
  // all of those are already computed in check_if_distinct!
  size_t prev_n_phi = pso->x_distinct_idx_of_last_batch;
  if (prev_n_phi == n_phi)
  {
    // There are no new points ! The surrogate is already fit !
#if DEBUG_SURROGATE
    printf("Skip fit_surrogate: no new evaluation position!\n");
#endif
    return 0;
  }
  else
  {
    pso->x_distinct_idx_of_last_batch = n_phi;
  }

  // Copy the distances from phi_cache to the phi block in A
  for (size_t j = 0; j < n_phi; j++)
  {
    double *d3_to_u_j_cached = phi_cache + j * (j - 1) / 2;

    for (size_t i = 0; i < j; i++)
    {
      double phi_i_j = d3_to_u_j_cached[i];
      A[i * n_A + j] = phi_i_j;
      A[j * n_A + i] = phi_i_j;
    }
    A[j * n_A + j] = 0.;
  }

  // P and tP are blocks in A
  // P_{i,j} := A_{i,n_phi + j} = A[i * n_A + n_phi + j]
  // tP_{i,j} := A_{n_phi + i, j} = A[(n_phi + i) * n_A + j]

  for (size_t k = 0; k < n_phi; k++)
  {
    double *u = x_distincts + k * dimensions;

    // P(p,0) = 1;
    A[k * n_A + n_phi + 0] = 1;
    // tP(0,p) = 1;
    A[(n_phi + 0) * n_A + k] = 1;

    for (int j = 0; j < dimensions; j++)
    {
      // P(p,1+j) = u[j];
      A[k * n_A + n_phi + j + 1] = u[j];
      // tP(1+j,p) = u[j];
      A[(n_phi + 1 + j) * n_A + k] = u[j];
    }
  }

  // lower right block is zeros
  for (size_t i = n_phi; i < n_A; i++)
  {
    for (size_t j = n_phi; j < n_A; j++)
    {
      A[i * n_A + j] = 0;
    }
  }

  /********
   * Prepare right hand side b
   ********/
  for (size_t k = 0; k < n_phi; k++)
  {
    // set b_k
    b[k] = fxd[k];
  }

  for (size_t k = n_phi; k < n_A; k++)
  {
    // set b_k
    b[k] = 0;
  }

#if DEBUG_SURROGATE
  print_matrixd(A, n_A, "A");
  print_vectord(b, n_A, "b");
#endif

  if (lu_solve(n_A, A, b) < 0)
  {
    return -1;
  }

  // b is overwitten with the result of Ax = b in lu_solve
  pso->lambda_p = b;

#if DEBUG_SURROGATE
  print_vectord(pso->lambda_p, n_A, "x");
#endif

  return 0;
}

/*
 * Cooperation between prealloc_fit_surrogate_6 and check_distinct
 */

int prealloc_fit_surrogate_6_BLOCK_TRI(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and n x 1 for b
  size_t Ab_size = max_n_A * (max_n_A + 1);
  size_t phi_cache_size = max_n_phi * (max_n_phi - 1) / 2;

  fit_surrogate_max_N_phi = max_n_phi;

  fit_surrogate_phi_cache = malloc(phi_cache_size * sizeof(double));

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  return 0;
}

// assumes Ab, n_Ab, n_phi, n_P are in scope.
// P_{i,j} := A_{i, j} = A[i * n_Ab + j]
// (P^t)_{i,j} := A_{n_phi + i, n_P + j} = A[(n_phi + i) * n_Ab + n_P + j]
#define BLOCK_TRI_P(i, j) (Ab[(i)*n_Ab + (j)])
#define BLOCK_TRI_Pt(i, j) (Ab[((i) + n_phi) * n_Ab + (j) + n_P])

#define BLOCK_TRI_Phi(i, j) (Ab[(i)*n_Ab + (j) + n_P])
#define BLOCK_TRI_Zeros(i, j) (Ab[((i) + n_phi) * n_Ab + (j)])

#define BLOCK_TRI_b(i) (Ab[(i)*n_Ab + n_A])

int fit_surrogate_6_BLOCK_TRI(struct pso_data_constant_inertia *pso)
{

  size_t dimensions = pso->dimensions;
  size_t popsize = pso->population_size;
  size_t t = pso->time;

  // the size of phi is the total number of _distinct_ points where
  // f has been evaluated
  // currently : n = x_distinct_s
  size_t n_phi = pso->x_distinct_s;
  double *x_distincts = pso->x_distinct;
  double *fxd = pso->x_distinct_eval;

  // the size of P is n x d+1
  size_t n_P = dimensions + 1;

  // the size of the matrix in the linear system is n+d+1
  size_t n_A = n_phi + n_P;
  size_t n_Ab = n_A + 1;

  double *Ab = fit_surrogate_Ab;

  size_t max_N_phi = fit_surrogate_max_N_phi;
  double *phi_cache = fit_surrogate_phi_cache;

  /********
   * Prepare Phi
   ********/

  // phi_p,q = || u_p - u_q ||^3
  // all of those are already computed in check_if_distinct!
  size_t prev_n_phi = pso->x_distinct_idx_of_last_batch;
  if (prev_n_phi == n_phi)
  {
    // There are no new points ! The surrogate is already fit !
#if DEBUG_SURROGATE
    printf("Skip fit_surrogate: no new evaluation position!\n");
#endif
    return 0;
  }
  else
  {
    pso->x_distinct_idx_of_last_batch = n_phi;
  }

  // Copy the distances from phi_cache to the phi block in A
  // Phi_{i, j} := A_{i, j + n_P} = A[i * n_Ab + n_P + j]
  for (size_t j = 0; j < n_phi; j++)
  {
    double *d3_to_u_j_cached = phi_cache + j * (j - 1) / 2;

    for (size_t i = 0; i < j; i++)
    {
      double phi_i_j = d3_to_u_j_cached[i];
      BLOCK_TRI_Phi(i, j) = phi_i_j;
      BLOCK_TRI_Phi(j, i) = phi_i_j;
    }
    BLOCK_TRI_Phi(j, j) = 0.;
  }

  /********
   * Prepare P and P^t
   ********/

  // P and tP are blocks in A

  for (size_t k = 0; k < n_phi; k++)
  {
    double *u = x_distincts + k * dimensions;

    // P(p,0) = 1;
    BLOCK_TRI_P(k, 0) = 1;
    // tP(0,p) = 1;
    BLOCK_TRI_Pt(0, k) = 1;

    for (int j = 0; j < dimensions; j++)
    {
      // P(p,1+j) = u[j];
      BLOCK_TRI_P(k, 1 + j) = u[j];
      // tP(1+j,p) = u[j];
      BLOCK_TRI_Pt(1 + j, k) = u[j];
    }
  }

  /********
   * Prepare 0 block of n_P x n_P
   ********/

  for (size_t i = 0; i < n_P; i++)
  {
    for (size_t j = 0; j < n_P; j++)
    {
      BLOCK_TRI_Zeros(i, j) = 0;
    }
  }

  /********
   * Prepare right hand side b
   ********/
  for (size_t k = 0; k < n_phi; k++)
  {
    // set b_k
    BLOCK_TRI_b(k) = 0;
  }

  for (size_t k = n_phi; k < n_A; k++)
  {
    // set b_k
    BLOCK_TRI_b(k) = fxd[k - n_phi];
  }

#if DEBUG_SURROGATE
  printf("matrix size: %d\n", n_Ab);
  print_rect_matrixd(Ab, n_A, n_Ab, "Ab");
#endif

  if (triangular_system_solve(n_A, n_P, Ab, pso->lambda_p) < 0)
  {
    return -1;
  }

#if DEBUG_SURROGATE
  print_vectord(pso->lambda_p, n_A, "x");
#endif

  return 0;
}

int fit_surrogate_6_LU_blocked(struct pso_data_constant_inertia *pso)
{
  size_t dimensions = pso->dimensions;
  size_t popsize = pso->population_size;
  size_t t = pso->time;

  // the size of phi is the total number of _distinct_ points where
  // f has been evaluated
  // currently : n = x_distinct_s
  size_t n_phi = pso->x_distinct_s;
  double *x_distincts = pso->x_distinct;
  double *fxd = pso->x_distinct_eval;

  // the size of P is n x d+1
  size_t n_P = dimensions + 1;

  // the size of the matrix in the linear system is n+d+1
  size_t n_A = n_phi + n_P;
  size_t n_Ab = n_A + 1;

  double *A = fit_surrogate_Ab;

  double *b = fit_surrogate_b;

  size_t max_N_phi = fit_surrogate_max_N_phi;
  double *phi_cache = fit_surrogate_phi_cache;
  /********
   * Prepare left hand side A
   ********/

  // phi_p,q = || u_p - u_q ||^3
  // all of those are already computed in check_if_distinct!
  size_t prev_n_phi = pso->x_distinct_idx_of_last_batch;
  if (prev_n_phi == n_phi)
  {
    // There are no new points ! The surrogate is already fit !
#if DEBUG_SURROGATE
    printf("Skip fit_surrogate: no new evaluation position!\n");
#endif
    return 0;
  }
  else
  {
    pso->x_distinct_idx_of_last_batch = n_phi;
  }

  // Copy the distances from phi_cache to the phi block in A

  /** auto-generated code **/

  /* Fill the matrix rowwise/colwise by blocks at the same time, starting at the
   * diagonal entry i */
  size_t i = 0;
  while (i + 4 - 1 < n_phi)
  {
    size_t i0 = i + 0;
    size_t i1 = i + 1;
    size_t i2 = i + 2;
    size_t i3 = i + 3;

    size_t cache_idx_row_i0 = i0 * (i0 - 1) / 2;
    size_t cache_idx_row_i1 = i1 * (i1 - 1) / 2;
    size_t cache_idx_row_i2 = i2 * (i2 - 1) / 2;
    size_t cache_idx_row_i3 = i3 * (i3 - 1) / 2;

    /* Fill in the corresponding blocks in the upper and lower triangle
     - these blocks are transposed of each other.
     - i is a multiple of 4
     - block on diagonal is special case
    */
    for (size_t j = 0; j < i; j += 4)
    {
      size_t j0 = j + 0;
      size_t j1 = j + 1;
      size_t j2 = j + 2;
      size_t j3 = j + 3;

      // load cached values
      double phi_0_0 = phi_cache[cache_idx_row_i0 + j0];
      double phi_0_1 = phi_cache[cache_idx_row_i0 + j1];
      double phi_0_2 = phi_cache[cache_idx_row_i0 + j2];
      double phi_0_3 = phi_cache[cache_idx_row_i0 + j3];
      double phi_1_0 = phi_cache[cache_idx_row_i1 + j0];
      double phi_1_1 = phi_cache[cache_idx_row_i1 + j1];
      double phi_1_2 = phi_cache[cache_idx_row_i1 + j2];
      double phi_1_3 = phi_cache[cache_idx_row_i1 + j3];
      double phi_2_0 = phi_cache[cache_idx_row_i2 + j0];
      double phi_2_1 = phi_cache[cache_idx_row_i2 + j1];
      double phi_2_2 = phi_cache[cache_idx_row_i2 + j2];
      double phi_2_3 = phi_cache[cache_idx_row_i2 + j3];
      double phi_3_0 = phi_cache[cache_idx_row_i3 + j0];
      double phi_3_1 = phi_cache[cache_idx_row_i3 + j1];
      double phi_3_2 = phi_cache[cache_idx_row_i3 + j2];
      double phi_3_3 = phi_cache[cache_idx_row_i3 + j3];

      // fill block in lower triangle
      A[i0 * n_A + j0] = phi_0_0;
      A[i0 * n_A + j1] = phi_0_1;
      A[i0 * n_A + j2] = phi_0_2;
      A[i0 * n_A + j3] = phi_0_3;
      A[i1 * n_A + j0] = phi_1_0;
      A[i1 * n_A + j1] = phi_1_1;
      A[i1 * n_A + j2] = phi_1_2;
      A[i1 * n_A + j3] = phi_1_3;
      A[i2 * n_A + j0] = phi_2_0;
      A[i2 * n_A + j1] = phi_2_1;
      A[i2 * n_A + j2] = phi_2_2;
      A[i2 * n_A + j3] = phi_2_3;
      A[i3 * n_A + j0] = phi_3_0;
      A[i3 * n_A + j1] = phi_3_1;
      A[i3 * n_A + j2] = phi_3_2;
      A[i3 * n_A + j3] = phi_3_3;

      // fill block in upper triangle
      A[j0 * n_A + i0] = phi_0_0;
      A[j0 * n_A + i1] = phi_1_0;
      A[j0 * n_A + i2] = phi_2_0;
      A[j0 * n_A + i3] = phi_3_0;
      A[j1 * n_A + i0] = phi_0_1;
      A[j1 * n_A + i1] = phi_1_1;
      A[j1 * n_A + i2] = phi_2_1;
      A[j1 * n_A + i3] = phi_3_1;
      A[j2 * n_A + i0] = phi_0_2;
      A[j2 * n_A + i1] = phi_1_2;
      A[j2 * n_A + i2] = phi_2_2;
      A[j2 * n_A + i3] = phi_3_2;
      A[j3 * n_A + i0] = phi_0_3;
      A[j3 * n_A + i1] = phi_1_3;
      A[j3 * n_A + i2] = phi_2_3;
      A[j3 * n_A + i3] = phi_3_3;
    }
    /* The last block (j = i) is on the diagonal, it is symmetric with a zero
     * diagonal. */
    {
      size_t j0 = i + 0;
      size_t j1 = i + 1;
      size_t j2 = i + 2;
      size_t j3 = i + 3;
      double phi_1_0 = phi_cache[cache_idx_row_i1 + j0];
      double phi_2_0 = phi_cache[cache_idx_row_i2 + j0];
      double phi_2_1 = phi_cache[cache_idx_row_i2 + j1];
      double phi_3_0 = phi_cache[cache_idx_row_i3 + j0];
      double phi_3_1 = phi_cache[cache_idx_row_i3 + j1];
      double phi_3_2 = phi_cache[cache_idx_row_i3 + j2];

      A[i0 * n_A + j0] = 0.;
      A[i0 * n_A + j1] = phi_1_0;
      A[i0 * n_A + j2] = phi_2_0;
      A[i0 * n_A + j3] = phi_3_0;

      A[i1 * n_A + j0] = phi_1_0;
      A[i1 * n_A + j1] = 0.;
      A[i1 * n_A + j2] = phi_2_1;
      A[i1 * n_A + j3] = phi_3_1;

      A[i2 * n_A + j0] = phi_2_0;
      A[i2 * n_A + j1] = phi_2_1;
      A[i2 * n_A + j2] = 0.;
      A[i2 * n_A + j3] = phi_3_2;

      A[i3 * n_A + j0] = phi_3_0;
      A[i3 * n_A + j1] = phi_3_1;
      A[i3 * n_A + j2] = phi_3_2;
      A[i3 * n_A + j3] = 0.;
    }
    i += 4;
  }

  /* finish the rows and the columns */
  while (i < n_phi)
  {
    size_t cache_idx_row_i = i * (i - 1) / 2;

    size_t j = 0;
    while (j + 4 - 1 < i)
    {
      size_t j0 = j + 0;
      size_t j1 = j + 1;
      size_t j2 = j + 2;
      size_t j3 = j + 3;
      double phi_i_0 = phi_cache[cache_idx_row_i + j0];
      double phi_i_1 = phi_cache[cache_idx_row_i + j1];
      double phi_i_2 = phi_cache[cache_idx_row_i + j2];
      double phi_i_3 = phi_cache[cache_idx_row_i + j3];

      // fill line in lower triangle
      A[i * n_A + j0] = phi_i_0;
      A[i * n_A + j1] = phi_i_1;
      A[i * n_A + j2] = phi_i_2;
      A[i * n_A + j3] = phi_i_3;

      // fill line in upper triangle
      A[j0 * n_A + i] = phi_i_0;
      A[j1 * n_A + i] = phi_i_1;
      A[j2 * n_A + i] = phi_i_2;
      A[j3 * n_A + i] = phi_i_3;
      j += 4;
    }
    while (j < i)
    {
      double phi_i_j = phi_cache[cache_idx_row_i + j];
      A[i * n_A + j] = phi_i_j;
      A[j * n_A + i] = phi_i_j;
      j++;
    }
    A[i * n_A + i] = 0.;
    i++;
  }
  /** end of auto-generated code **/

  // P and tP are blocks in A
  // P_{i,j} := A_{i,n_phi + j} = A[i * n_A + n_phi + j]
  // tP_{i,j} := A_{n_phi + i, j} = A[(n_phi + i) * n_A + j]

  for (size_t k = 0; k < n_phi; k++)
  {
    double *u = x_distincts + k * dimensions;

    // P(p,0) = 1;
    A[k * n_A + n_phi + 0] = 1;
    // tP(0,p) = 1;
    A[(n_phi + 0) * n_A + k] = 1;

    for (int j = 0; j < dimensions; j++)
    {
      // P(p,1+j) = u[j];
      A[k * n_A + n_phi + j + 1] = u[j];
      // tP(1+j,p) = u[j];
      A[(n_phi + 1 + j) * n_A + k] = u[j];
    }
  }

  // lower right block is zeros
  for (size_t i = n_phi; i < n_A; i++)
  {
    for (size_t j = n_phi; j < n_A; j++)
    {
      A[i * n_A + j] = 0;
    }
  }

  /********
   * Prepare right hand side b
   ********/
  for (size_t k = 0; k < n_phi; k++)
  {
    // set b_k
    b[k] = fxd[k];
  }

  for (size_t k = n_phi; k < n_A; k++)
  {
    // set b_k
    b[k] = 0;
  }

#if DEBUG_SURROGATE
  print_matrixd(A, n_A, "A");
  print_vectord(b, n_A, "b");
#endif

  PAPI_START("system_solver");
  int ret = lu_solve(n_A, A, b);
  PAPI_STOP("system_solver");

  if (ret < 0)
  {
    return -1;
  }

  // b is overwitten with the result of Ax = b in lu_solve
  pso->lambda_p = b;

#if DEBUG_SURROGATE
  print_vectord(pso->lambda_p, n_A, "x");
#endif

  return 0;
}