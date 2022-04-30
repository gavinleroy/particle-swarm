#include "pso.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>


#include "local_refinement.h"
#include "plu_factorization.h"


static double rand_between(double a, double b)
{
    // see http://c-faq.com/lib/randrange.html
    return a + (double)rand() / ((double)RAND_MAX / (b - a));
}


static double clamp(double v, double lo, double hi)
{
    if (v < lo)
    {
        return lo;
    }
    else if (v > hi)
    {
        return hi;
    }
    else
    {
        return v;
    }
}

static double dist2(size_t dim, double const * x, double const * y)
{
    double s = 0;
    for (size_t i = 0 ; i < dim ; i++)
    {
        double v = x[i] - y[i];
        s += v * v;
    }
    return s;
}

static double dist(size_t dim, double const * x, double const * y)
{
    return sqrt(dist2(dim, x, y));
}


// for other inertia choices, see https://ieeexplore.ieee.org/document/6089659
struct pso_data_constant_inertia
{
    blackbox_fun f;
    // positions x_i, saved for all times
    // i.e. x[t][i] is the position vector x_i(t)
    double *** x;
    // f(x_i), for all times
    double ** x_eval;

    // velocity v_i
    double ** v;
    // best position recorded for particle i y_i
    double ** y;
    double * y_eval;
    
    double * y_hat;

    // Keep track of local refinement points.
    // max lenght of list = tmax
    // current length n_past_refinement_points
    double ** past_refinement_points;
    double * past_refinement_points_eval; 


    double * bound_low;
    double * bound_high;
    double * vmin;
    double * vmax;

    // parameters of the surrogate
    double * lambda;
    double * p;

    double inertia;
    double social;
    double cognition;

    double local_refinement_box_size;
    double min_minimizer_distance;

    double y_hat_eval;
    
    int dimensions;
    int population_size;
    int n_trials;

    int n_past_refinement_points;

    int time_max;
    int time;
};




double surrogate_eval(struct pso_data_constant_inertia const * pso, double const * x)
{
    //TODO: add the past_local_refinements
    
    double res = 0;
    
    int p = 0;
    for (int t = 0 ; t < pso->time + 1 ; t++)
    {
        for (int i = 0 ; i < pso->population_size ; i++)
        {
            double * u = pso->x[t][i];
            double d = dist(pso->dimensions, u, x);
            res += pso->lambda[p] * d;


            p++;
        }
    }

    for (int j = 0 ; j < pso->dimensions ; j++)
    {
        res += pso->p[j+1] * x[j];
    }
    res += pso->p[0];

    return res;
}

// surrogate inferface taking a void pointer
// to pass to local_optimizer
double surrogate_eval_void(double const * x, void const * args)
{
    struct pso_data_constant_inertia const * pso = args;
    return surrogate_eval(pso, x);
}


void fit_surrogate(struct pso_data_constant_inertia * pso)
{

    //TODO: include past_refinement_points in phi !!!

    //TODO: note that the matrix and vector barely change between the
    // iterations. Maybe there could be a way to re-use them?

    // the size of phi is the total number of points where
    // f has been evaluated
    // currently : n = population_size * (time + 1)
    size_t n_phi = pso->population_size * (pso->time + 1); // + pso->n_past_refinement_points

    // the size of P is n x d+1
    size_t n_P = pso->dimensions + 1;

    // the size of the matrix in the linear system is n+d+1
    size_t n_A = n_phi + n_P;

    double * A = malloc(n_A * n_A * sizeof(double));

    // phi_p,q = || u_p - u_q ||
    // currently the {u_p} = {x_i(t=j)} i=0..pop_size, j=0..time+1

    size_t p = 0;
    for(int t1=0 ; t1 < pso->time+1 ; t1++)
    {
        for (int i1 = 0 ; i1 < pso->population_size ; i1++)
        {
            double * u_p = pso->x[t1][i1];

            size_t q = 0;
            for (int t2 = 0 ; t2 < pso->time+1 ; t2++)
            {
                for (int i2 = 0 ; i2 < pso->population_size ; i2++)
                {
                    double * u_q = pso->x[t2][i2];

                    // Phi is the top left block: phi_ij := A_ij
                    A[p * n_A + q] = dist(pso->dimensions, u_p, u_q);

                    q++;
                }
            }
            p++;
        }
    }


    // P and tP are blocks in A
    // P_{i,j} := A_{i,n_phi + j} = A[i * n_A + n_phi + j]
    // tP_{i,j} := A_{n_phi + i, j} = A[(n_phi + i) * n_A + j]

    p = 0;
    for (int t = 0 ; t < pso->time + 1 ; t++)
    {
        for (int i = 0 ; i < pso->population_size ; i++)
        {
            double * u = pso->x[t][i];
        
            // P(p,0) = 1;
            A[p * n_A + n_phi + 0] = 1;
            // tP(0,p) = 1;
            A[(n_phi + 0) * n_A + p] = 1;

            for (int j = 0 ; j < pso->dimensions ; j++)
            {
                //P(p,1+j) = u[j];
                A[p * n_A + n_phi + j + 1] = u[j];
                //tP(1+j,p) = u[j];
                A[(n_phi + 1 + j) * n_A + p] = u[j];
            }
            p++;
        }
    }

    // lower right block is zeros
    for(int i = n_phi ; i < n_A ; i++)
    {
        for (int j = n_phi ; j < n_A ; j++)
        {
            A[i * n_A + j] = 0;
        }
    }

    // right hand side
    double * b = malloc(n_A * sizeof(double));

    p = 0;
    for (int t = 0 ; t < pso->time_max + 1 ; t++)
    {
        for (int i = 0 ; i < pso->population_size ; i ++)
        {
            b[p] = pso->x_eval[t][i];
        }
    }
    for (int i = n_phi ; i < n_A ; i++)
    {
        b[i] = 0;
    }

    // solve A x = b using partial pivotting LU
    plu_factorization plu;
    alloc_plu_factorization(pso->dimensions, &plu);
    plu_factorize(pso->dimensions, A, &plu);

    double * x = malloc(n_A * sizeof(double));
    plu_solve(pso->dimensions, &plu, b, x);

    realloc(pso->lambda, n_phi);
    for (int i = 0 ; i < n_phi ; i++)
    {
        pso->lambda[i] = x[i];
    }


    for(int i = 0 ; i < n_P ; i++)
    {
        pso->p[i] = x[n_phi + i];
    }

    free(x);
    free_plu_factorization(&plu);
    free(b);
    free(A);
}






void pso_constant_inertia_init(
    struct pso_data_constant_inertia * pso,
    blackbox_fun f,
    double inertia,
    double social, double cognition,
    double local_refinement_box_size,
    double min_minimizer_distance,
    int dimensions,
    int population_size, int time_max, int n_trials,
    double * bounds_low, double * bounds_high,
    double * vmin, double * vmax,
    double ** initial_positions
)
{
    pso->f = f;
    pso->inertia = inertia;
    pso->dimensions = dimensions;
    pso->social = social, pso->cognition = cognition;
    pso->local_refinement_box_size = local_refinement_box_size;
    pso->min_minimizer_distance = min_minimizer_distance;
    pso->vmin = vmin, pso->vmax = vmax;
    pso->population_size = population_size;
    pso->time_max = time_max, pso->n_trials = n_trials;

    pso->n_past_refinement_points = 0;
    pso->time = 0;


    pso->x = (double***)malloc(pso->time_max * sizeof(double **));
    for (int t = 0 ; t < pso->time_max ; t++)
    {
        pso->x[t] = (double**)malloc(pso->population_size * sizeof(double*));
        // vector gets allocated for trial anyway
        //for (int i = 0 ; i < pso->population_size ; i++)
        //{
        //    pso->x[t][i] = (double*)malloc(pso->dimensions * sizeof(double));
        //}
    }

    pso->x_eval = (double**)malloc(pso->time_max * sizeof(double *));
    for (int t = 0 ; t < pso->time_max ; t++)
    {
        pso->x_eval[t] = (double*)malloc(pso->population_size * sizeof(double));
    }

    pso->v = (double**)malloc(pso->population_size * sizeof(double *));
    for (int i = 0 ; i < pso->population_size ; i++)
    {
        pso->v[i] = (double*)malloc(pso->dimensions * sizeof(double));
    }

    pso->y = (double**)malloc(pso->population_size * sizeof(double *));
    for (int i = 0 ; i < pso->population_size ; i++)
    {
        pso->y[i] = (double*)malloc(pso->dimensions * sizeof(double));
    }

    pso->y_eval= (double*)malloc(pso->population_size * sizeof(double *));

    pso->past_refinement_points = malloc(pso->time_max * sizeof(double*));
    // alloc-ed in the loop
    //for (int i = 0 ; i < pso->time_max ; i++)
    //{
    //    pso->past_refinement_points[i] = malloc(pso->dimensions * sizeof(double));
    //}

    pso->past_refinement_points_eval = malloc(pso->population_size * sizeof(double));


    pso->bound_low = (double*)malloc(pso->dimensions * sizeof(double));
    pso->bound_high = (double*)malloc(pso->dimensions * sizeof(double));

    pso->vmin = (double*)malloc(pso->dimensions * sizeof(double));
    pso->vmax = (double*)malloc(pso->dimensions * sizeof(double));


    // will realloc in fit_surrogate
    pso->lambda = NULL;
    pso->p = malloc((pso->dimensions + 1) * sizeof(double));


    // setup x
    for (int i = 0 ; i < population_size ; i ++)
    {
        for (int j = 0 ; j < pso->dimensions ; j++)
        {
            pso->x[0][i][j] = initial_positions[i][j];
        }
    }

    // setup bounds in space
    for(int k = 0 ; k < pso->dimensions ; k++)
        pso->bound_low[k] = bounds_low[k];
    for(int k = 0 ; k < pso->dimensions ; k++)
        pso->bound_high[k] = bounds_high[k];
    // setup bounds on velocity
    for(int k = 0 ; k < pso->dimensions ; k++)
        pso->vmin[k] = vmin[k];
    for(int k = 0 ; k < pso->dimensions ; k++)
        pso->vmax[k] = vmax[k];
}

int is_far_from_previous_evaluations(struct pso_data_constant_inertia const * pso, double * x)
{
    double delta2 = pso->min_minimizer_distance * pso->min_minimizer_distance;

    // check previous particle positions
    for(int t = 0 ; t < pso->time + 1 ; t++)
    {
        for(int i = 0 ; i < pso->population_size ; i++)
        {
            double d2 = dist2(pso->dimensions, x, pso->x[t][i]);
            if (d2 < delta2)
            {
                return 0;
            }
        }
    }

    // check previous local minimizers positions
    for(int k = 0 ; k < pso->n_past_refinement_points ; k++)
    {
        double d2 = dist2(pso->dimensions, x, pso->past_refinement_points[k]);
        if (d2 < delta2)
        {
            return 0;
        } 
    }

    return 1;
}


void pso_constant_inertia_first_steps(struct pso_data_constant_inertia * pso)
{
    //TODO: step 1 and 2 "space-filling design"

    // Step 3. Initialize particle velocities
    for(int i = 0 ; i < pso->population_size ; i++)
    {
        for (int k=0 ; k < pso->dimensions ; k++)
        {
            double uk = rand_between(pso->bound_low[k], pso->bound_high[k]);
            pso->v[i][k] = 1. / 2 * (uk - pso->x[0][i][k]);
        }
    }

    // Step 4. Initialise y, y_eval, and x_eval for each particle 
    for (int i = 0 ; i < pso->population_size ; i++)
    {
        for (int k = 0 ; k < pso->dimensions ;  k++)
        {
            pso->y[i][k] = pso->x[0][i][k];
        }

        double x_eval = pso->f(pso->x[0][i]);
        pso->y_eval[i] = x_eval;
        pso->x_eval[0][i] = x_eval;

    }

    // find ŷ
    
    double * y_hat = pso->y[0];
    double y_hat_eval = pso->y_eval[0];
    for (int i = 1 ; i < pso->population_size ; i++)
    {
        if (pso->y_eval[i] < y_hat_eval)
        {
            y_hat = pso->y[i];
            y_hat_eval = pso->y_eval[i];
        }
    }
    pso->y_hat = y_hat;
    pso->y_hat_eval = y_hat_eval;
}


bool pso_constant_inertia_loop(struct pso_data_constant_inertia * pso)
{
    int const t = pso->time;


    //TODO: update steps below from normal to OPUS

    // Step 5.
    // Fit surrogate
    // f already evaluated on x[0..t][0..i-1]
    fit_surrogate(pso);


    // Step 6
    // Determine new particle positions
    double * v_trial = malloc(pso->dimensions * sizeof(double));
    double * x_trial = malloc(pso->dimensions * sizeof(double));
    
    double * v_trial_best = malloc(pso->dimensions * sizeof(double));
    double * x_trial_best = malloc(pso->dimensions * sizeof(double));
    double x_trial_best_seval = DBL_MAX;

    for(int i = 0 ; i < pso->population_size ; i++)
    {
        for (int l = 0 ; l < pso->n_trials ; l++)
        {
            
            for (int j = 0 ; j < pso->dimensions ; j++)
            {
                // compute v_i(t+1) from v_i(t)

                double w1 = (double)rand() / RAND_MAX;
                double w2 = (double)rand() / RAND_MAX;
                double v =
                    pso->inertia * pso->v[i][j] +
                    pso->cognition * w1 * (pso->y[i][j] - pso->x[t][i][j]) +
                    pso->social * w2 * (pso->y_hat[j] - pso->x[t][i][j]);

                v_trial[j] = clamp(v, pso->vmin[j], pso->vmax[j]);


                x_trial[j] = clamp(pso->x[pso->time][i][j] + v_trial[j], pso->bound_low[j], pso->bound_high[j]);
            }

            double x_trial_seval = surrogate_eval(pso, x_trial);
            if(x_trial_seval < x_trial_best_seval)
            {
                // keep x_trial as x_trial_best by swapping the two buffers: the new x_trial
                // will get overwritten in the next iteration
                x_trial_best_seval = x_trial_seval;
                
                double * t;
                
                t = x_trial;
                x_trial = x_trial_best;
                x_trial_best = t;

                t = v_trial;
                v_trial = v_trial_best;
                v_trial_best = t;
            }
        }

        // set next position and update velocity
        // note: x[t+1][i] was unset, but v[i] is freed
        // before being replaced by v_trial
        pso->x[t+1][i] = x_trial_best;
        free(pso->v[i]);
        pso->v[i] = v_trial_best;

        free(x_trial);
        free(v_trial);
    }






    // Step 7
    // Evaluate swarm positions
    for (int i = 0 ; i < pso->population_size ; i++)
    {
        pso->x_eval[t+1][i] = pso->f(pso->x[t+1][i]);
    }

    // Step 8. Update the best positions per particle and overall

    for (int i = 0 ; i < pso->population_size ; i++)
    {
        double x_eval = pso->x_eval[t+1][i];

        if (x_eval < pso->y_eval[i])
        {
            // y_i <- x_i
            pso->y[i] = pso->x[t+1][i];
            pso->y_eval[i] = x_eval;

            // is x_i(t+1) better than ŷ(t) ?
            if (x_eval < pso->y_hat_eval)
            {
                pso->y_hat = pso->y[i];
            }
        }
    }

    pso->time++;

    // Step 9
    // Refit surrogate with time = t+1
    fit_surrogate(pso);


    // Step 10
    // Local refinement
    double * x_local = malloc(pso->dimensions * sizeof(double));

    local_optimization(
        &surrogate_eval_void,
        pso->dimensions,
        pso->y_hat,
        pso->local_refinement_box_size,
        pso->bound_low,
        pso->bound_high,
        pso,
        x_local
    );
    
    // Step 11
    // Determine if minimizer of surrogate is far from previous points
    if(is_far_from_previous_evaluations(pso, x_local))
    {
        double x_local_eval = pso->f(x_local);

        pso->past_refinement_points[pso->n_past_refinement_points] = x_local;
        pso->past_refinement_points_eval[pso->n_past_refinement_points] = x_local_eval;
        pso->n_past_refinement_points++;

        // update overall best if applicable
        if(x_local_eval < pso->y_hat_eval)
        {
            pso->y_hat = x_local;
            pso->y_hat_eval = x_local_eval;
        }
    }
    else
    {
        free(x_local);
    }


    return (pso->time < pso->time_max);
}


void run_pso(
    blackbox_fun f,
    double inertia,
    double social, double cognition,
    double local_refinement_box_size,
    double min_minimizer_distance,
    int dimensions,
    int population_size, int time_max, int n_trials,
    double * bounds_low, double * bounds_high,
    double * vmin, double * vmax,
    double ** initial_positions
)
{
    struct pso_data_constant_inertia pso;
    pso_constant_inertia_init(
        &pso,
        f,
        inertia, social, cognition,
        local_refinement_box_size,
        min_minimizer_distance,
        dimensions,
        population_size, time_max, n_trials,
        bounds_low, bounds_high,
        vmin, vmax,
        initial_positions
    );

    pso_constant_inertia_first_steps(&pso);

    do {
        printf("t=%d  ŷ=[", pso.time);
        for (int j = 0 ; j < dimensions ; j++)
        {
            printf("%f", pso.y_hat[j]);
            if (j < dimensions - 1) printf(", ");
        }
        printf("]  f(ŷ)=%f\n", pso.y_hat_eval);
    } while(pso_constant_inertia_loop(&pso));
}