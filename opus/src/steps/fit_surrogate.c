#include "fit_surrogate.h"

#include "math.h"
#include "stdlib.h"
#include "string.h"

#include "../gaussian_elimination_solver.h"
#include "../helpers.h"
#include "../pso.h"

#define QUOTE(x) #x
#define STR(x) QUOTE(x)

int fit_surrogate(struct pso_data_constant_inertia *pso)
{
  fit_surrogate_4(pso);
}

int prealloc_fit_surrogate(size_t max_n_phi, size_t n_P)
{
  prealloc_fit_surrogate_4(max_n_phi, n_P);
}

#define DEBUG_SURROGATE 0

// TODO: include past_refinement_points in phi !!!

static double *fit_surrogate_Ab;
static double *fit_surrogate_x;
static double *fit_surrogate_P;

static size_t fit_surrogate_max_N_phi;
static double *fit_surrogate_phi_cache;

int prealloc_fit_surrogate_0(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and n x 1 for b
  size_t Ab_size = max_n_A * max_n_A + max_n_A;

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_x = malloc(max_n_A * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  return 0;
}

int fit_surrogate_0(struct pso_data_constant_inertia *pso)
{
  // TODO: note that the matrix and vector barely change between the
  //  iterations. Maybe there could be a way to re-use them?

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

  // phi_p,q = || u_p - u_q ||
  // currently the {u_p} = {x_i(t=j)} i=0..pop_size, j=0..time+1

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

  double *x = fit_surrogate_x;

  if (gaussian_elimination_solve(n_A, Ab, pso->lambda_p) < 0)
  {
    return -1;
  }

#if DEBUG_SURROGATE
  print_vectord(pso->lambda_p, n_A, "x");
#endif
}

int prealloc_fit_surrogate_1(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and n x 1 for b
  size_t Ab_size = max_n_A * max_n_A + max_n_A;

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_x = malloc(max_n_A * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  return 0;
}

int fit_surrogate_1(struct pso_data_constant_inertia *pso)
{
  // TODO: include past_refinement_points in phi !!!

  // TODO: note that the matrix and vector barely change between the
  //  iterations. Maybe there could be a way to re-use them?

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
  double *x = fit_surrogate_x;

  /********
   * Prepare left hand side A
   ********/

  // phi_p,q = || u_p - u_q ||
  // currently the {u_p} = {x_i(t=j)} i=0..pop_size, j=0..time+1

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

int prealloc_fit_surrogate_2(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and n x 1 for b
  size_t Ab_size = max_n_A * max_n_A + max_n_A;

  fit_surrogate_max_N_phi = max_n_phi;
  fit_surrogate_phi_cache = malloc(max_n_phi * max_n_phi * sizeof(double));

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_x = malloc(max_n_A * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  return 0;
}

int fit_surrogate_2(struct pso_data_constant_inertia *pso)
{
  // TODO: include past_refinement_points in phi !!!

  // TODO: note that the matrix and vector barely change between the
  //  iterations. Maybe there could be a way to re-use them?

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
  double *x = fit_surrogate_x;

  size_t max_N_phi = fit_surrogate_max_N_phi;
  double *phi_cache = fit_surrogate_phi_cache;
  /********
   * Prepare left hand side A
   ********/

  // phi_p,q = || u_p - u_q ||
  // currently the {u_p} = {x_i(t=j)} i=0..pop_size, j=0..time+1

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

int prealloc_fit_surrogate_3(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and n x 1 for b
  size_t Ab_size = max_n_A * max_n_A + max_n_A;

  fit_surrogate_max_N_phi = max_n_phi;
  fit_surrogate_phi_cache = malloc(max_n_phi * max_n_phi * sizeof(double));

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_x = malloc(max_n_A * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  return 0;
}

int fit_surrogate_3(struct pso_data_constant_inertia *pso)
{
  // TODO: include past_refinement_points in phi !!!

  // TODO: note that the matrix and vector barely change between the
  //  iterations. Maybe there could be a way to re-use them?

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
  double *x = fit_surrogate_x;

  size_t max_N_phi = fit_surrogate_max_N_phi;
  double *phi_cache = fit_surrogate_phi_cache;
  /********
   * Prepare left hand side A
   ********/

  // phi_p,q = || u_p - u_q ||
  // currently the {u_p} = {x_i(t=j)} i=0..pop_size, j=0..time+1

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

int prealloc_fit_surrogate_4(size_t max_n_phi, size_t n_P)
{
  size_t max_n_A = max_n_phi + n_P;
  // Ab size: n x n for A and n x 1 for b
  size_t Ab_size = max_n_A * max_n_A + max_n_A;

  fit_surrogate_max_N_phi = max_n_phi;
  fit_surrogate_phi_cache = malloc(max_n_phi * max_n_phi * sizeof(double));

  fit_surrogate_Ab = malloc(Ab_size * sizeof(double));
  fit_surrogate_x = malloc(max_n_A * sizeof(double));
  fit_surrogate_P = malloc(max_n_phi * n_P * sizeof(double));
  return 0;
}

int fit_surrogate_4(struct pso_data_constant_inertia *pso)
{
  // TODO: include past_refinement_points in phi !!!

  // TODO: note that the matrix and vector barely change between the
  //  iterations. Maybe there could be a way to re-use them?

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
  double *x = fit_surrogate_x;

  size_t max_N_phi = fit_surrogate_max_N_phi;
  double *phi_cache = fit_surrogate_phi_cache;
  /********
   * Prepare left hand side A
   ********/

  // phi_p,q = || u_p - u_q ||
  // currently the {u_p} = {x_i(t=j)} i=0..pop_size, j=0..time+1

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
