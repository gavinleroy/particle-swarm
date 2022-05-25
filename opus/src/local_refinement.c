
#include "local_refinement.h"

#include <float.h>
#include <stdlib.h>

void local_optimization(local_optimization_function f, // R^d -> R
                        size_t dimensions,             // R
                        double const *center,          // R^d
                        double xi,                     // R
                        double const *a,               // R^d
                        double const *b,               // R^d
                        void const *const additionnal_f_args, double *x_min)
{
  // hyperparameter
  const size_t divisions = 10;

  /* Performing the intersection:
      [lo; hi] = [x-eta; x+eta] \and [a; b]
  */
  double *space_lo = (double *)malloc(sizeof(double) * dimensions);
  double *space_hi = (double *)malloc(sizeof(double) * dimensions);

  for (size_t it = 0; it < dimensions; ++it)
  {
    if (a[it] <= center[it] - xi)
      space_lo[it] = center[it] - xi;
    else
      space_lo[it] = a[it];
  }

  for (size_t it = 0; it < dimensions; ++it)
  {
    if (b[it] >= center[it] + xi)
      space_lo[it] = center[it] + xi;
    else
      space_lo[it] = b[it];
  }

  /*  Dividing the space [lo; hi] into a grid.
      Each grid cube has side length (hi - lo) / divisions
      We then take the center as the representative of the cell
  */

  // c_k = grid_centers[k * divisions .. k * divisions + dimensions]
  double *grid_centers =
      (double *)malloc(divisions * dimensions * sizeof(double));
  for (size_t div_it = 0; div_it < divisions; ++div_it)
  {
    for (size_t dim_it = 0; dim_it < dimensions; ++dim_it)
    {
      double interval_length =
          (space_hi[dim_it] - space_lo[dim_it]) / divisions;

      grid_centers[div_it * dimensions + dim_it] =
          space_lo[dim_it] + interval_length * (div_it + 0.5);
    }
  }

  /* We search for the grid_center that has the highest value under f
      Such grid_center will be our x_opt (return value)
  */
  size_t best_center_index = 0;
  double best_center_value = DBL_MAX;
  for (size_t div_it = 0; div_it < divisions; ++div_it)
  {
    double v = f(&grid_centers[div_it * dimensions], additionnal_f_args);
    if (v < best_center_value)
    {
      best_center_index = div_it;
      best_center_value = v;
    }
  }

  // copy best center to output
  for (int i = 0; i < dimensions; i++)
  {
    x_min[i] = grid_centers[best_center_index * dimensions + i];
  }

  // De-allocate

  free(space_lo);
  free(space_hi);
  free(grid_centers);
}
