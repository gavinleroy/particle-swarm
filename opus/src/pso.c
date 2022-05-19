#include "pso.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#include "helpers.h"
#include "local_refinement.h"
#include "plu_factorization.h"


#define DEBUG_TRIALS 0
#define DEBUG_SURROGATE 0

// #define DEBUG_STDOUT

// #if USE_ROUNDING_BLOOM_FILTER
#include "rounding_bloom.h"
// #endif

#define LOG_SURROGATE 1

#if LOG_SURROGATE
#include "logging.h"
#endif

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

int is_far_from_previous_evaluations(struct pso_data_constant_inertia const * pso, double * x, double min_dist)
{
    // - find way to check "if minimizer of surrogate is far from previous points"
    //    + see https://en.wikipedia.org/wiki/Nearest_neighbor_search
    //    + in high dim naive search can be best  

    double delta2 = min_dist * min_dist;

    // check previous particle positions
    for(int t = 0 ; t < pso->time + 1 ; t++)
    {
        for(int i = 0 ; i < pso->population_size ; i++)
        {
            double d2 = dist2(pso->dimensions, x, PSO_X(pso, t, i));
            if (d2 < delta2)
            {
                return 0;
            }
        }
    }

    // check previous local minimizers positions
    for(int k = 0 ; k < pso->n_past_refinement_points ; k++)
    {
        double d2 = dist2(pso->dimensions, x, PSO_PAST_REFINEMENT(pso, k));
        if (d2 < delta2)
        {
            return 0;
        } 
    }

    return 1;
}


void is_x_distinct(int dim, size_t point_cloud_s, double const * point_cloud, double const * x)
{
}

double surrogate_eval(struct pso_data_constant_inertia const * pso, double const * x)
{
    //TODO: add the past_local_refinements
    //Note: will also require adding them to the bloom filter / distinctiveness check
    
    //TODO: I think there is something wrong with the surrogate:
    // I don't have the same values when evaluating in the Jupyter notebook
    // To investigate!
    // VO

    double res = 0;
    
    
    for (int k = 0 ; k < pso->x_distinct_s ; k++)
    {
        size_t p = pso->x_distinct[k];
        double * u = pso->x + p * pso->dimensions;
        double d = dist(pso->dimensions, u, x);
        res += pso->lambda[k] * d * d * d;
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


int fit_surrogate(struct pso_data_constant_inertia * pso)
{
    int ret = 0;
    //TODO: include past_refinement_points in phi !!!

    //TODO: note that the matrix and vector barely change between the
    // iterations. Maybe there could be a way to re-use them?

    // the size of phi is the total number of _distinct_ points where
    // f has been evaluated
    // currently : n = x_distinct_s
    size_t n_phi = pso->x_distinct_s; // + pso->n_past_refinement_points

    // the size of P is n x d+1
    size_t n_P = pso->dimensions + 1;

    // the size of the matrix in the linear system is n+d+1
    size_t n_A = n_phi + n_P;

    double * A = malloc(n_A * n_A * sizeof(double));

    // phi_p,q = || u_p - u_q ||
    // currently the {u_p} = {x_i(t=j)} i=0..pop_size, j=0..time+1

    for (size_t k1 = 0 ; k1 < pso->x_distinct_s ; k1++)
    {
        size_t p = pso->x_distinct[k1];
        double * u_p = pso->x + p * pso->dimensions;

        for(size_t k2 = 0 ; k2 < pso->x_distinct_s ; k2++)
        {
            size_t q = pso->x_distinct[k2];
            double * u_q = pso->x + q * pso->dimensions;

            A[k1 * n_A + k2] = dist(pso->dimensions, u_p, u_q);
        }
    }



    // P and tP are blocks in A
    // P_{i,j} := A_{i,n_phi + j} = A[i * n_A + n_phi + j]
    // tP_{i,j} := A_{n_phi + i, j} = A[(n_phi + i) * n_A + j]


    for (size_t k = 0 ; k < pso->x_distinct_s ; k++)
    {
        size_t p = pso->x_distinct[k];
        double * u = pso->x + p * pso->dimensions;
        
        // P(p,0) = 1;
        A[k * n_A + n_phi + 0] = 1;
        // tP(0,p) = 1;
        A[(n_phi + 0) * n_A + k] = 1;

        for (int j = 0 ; j < pso->dimensions ; j++)
        {
            //P(p,1+j) = u[j];
            A[k * n_A + n_phi + j + 1] = u[j];
            //tP(1+j,p) = u[j];
            A[(n_phi + 1 + j) * n_A + k] = u[j];
        }
    }

    // lower right block is zeros
    for(size_t i = n_phi ; i < n_A ; i++)
    {
        for (size_t j = n_phi ; j < n_A ; j++)
        {
            A[i * n_A + j] = 0;
        }
    }


    // right hand side
    double * b = malloc(n_A * sizeof(double));

    for (size_t k = 0 ; k < pso->x_distinct_s ; k++)
    {
        size_t p = pso->x_distinct[k];
        b[k] = pso->x_eval[p];
    }

    for (size_t i = n_phi ; i < n_A ; i++)
    {
        b[i] = 0;
    }

    #if DEBUG_SURROGATE
    print_matrixd(A, n_A, "A");
    print_vectord(b, n_A, "b");
    #endif

    // solve A x = b using partial pivotting LU
    plu_factorization plu;
    alloc_plu_factorization(n_A, &plu);
    if (plu_factorize(n_A, A, &plu) < 0)
    {
        ret = -1;
        goto fit_surrogate_plu_factorization_failure;
    }

    double * x = malloc(n_A * sizeof(double));
    plu_solve(n_A, &plu, b, x);


    #if DEBUG_SURROGATE
    print_vectord(x, n_A, "x");
    #endif

    pso->lambda = realloc(pso->lambda, n_phi * sizeof(double));
    for (size_t i = 0 ; i < n_phi ; i++)
    {
        pso->lambda[i] = x[i];
    }


    for(size_t i = 0 ; i < n_P ; i++)
    {
        pso->p[i] = x[n_phi + i];
    }

    free(x);
fit_surrogate_plu_factorization_failure:
    free_plu_factorization(&plu);
    free(b);
    free(A);

    return ret;
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
    pso->population_size = population_size;
    pso->time_max = time_max, pso->n_trials = n_trials;

    pso->n_past_refinement_points = 0;
    pso->time = 0;


    pso->x = malloc(pso->time_max * pso->population_size * pso->dimensions * sizeof(double));

    pso->x_eval = malloc(pso->time_max * pso->population_size *  sizeof(double));

    pso->v = malloc(pso->population_size * pso->dimensions * sizeof(double));

    pso->y = malloc(pso->population_size * pso->dimensions * sizeof(double));

    pso->y_eval= malloc(pso->population_size * sizeof(double *));

    pso->past_refinement_points = malloc(pso->time_max * pso->dimensions * sizeof(double));
    pso->past_refinement_points_eval = malloc(pso->population_size * sizeof(double));


    pso->bound_low = (double*)malloc(pso->dimensions * sizeof(double));
    pso->bound_high = (double*)malloc(pso->dimensions * sizeof(double));

    pso->vmin = (double*)malloc(pso->dimensions * sizeof(double));
    pso->vmax = (double*)malloc(pso->dimensions * sizeof(double));

    for (int j = 0 ; j < dimensions ; j++)
    {
        pso->vmin[j] = vmin[j];
        pso->vmax[j] = vmax[j];
    }

    pso->x_distinct = malloc(pso->time_max * pso->population_size * sizeof(size_t));
    pso->x_distinct_s = 0;

    #if USE_ROUNDING_BLOOM_FILTER
    pso->bloom = malloc(sizeof(struct rounding_bloom));
    int bloom_entries = pso->time_max * pso->population_size;
    if (bloom_entries < 1000) bloom_entries = 1000;
    double bloom_false_pos_rate = 0.001;
    double bloom_rounding_eps = 0.1;
    rounding_bloom_init(pso->bloom, bloom_entries, bloom_false_pos_rate, bloom_rounding_eps, dimensions, bounds_low);
    #else
    // No data structure needed for naive search
    #endif

    // will realloc in fit_surrogate
    pso->lambda = NULL;
    pso->p = malloc((pso->dimensions + 1) * sizeof(double));


    // setup x
    for (int i = 0 ; i < population_size ; i ++)
    {
        // add it to the point cloud
        for (int j = 0 ; j < pso->dimensions ; j++)
        {
            PSO_X(pso, 0, i)[j] = initial_positions[i][j];
        }

        // Check if x is distinct
        #if USE_ROUNDING_BLOOM_FILTER
        // add and check proximity to previous points
        if (! rounding_bloom_check_add(pso->bloom, dimensions, PSO_X(pso, 0, i), 1))
        {
            pso->x_distinct[pso->x_distinct_s] = i;
            pso->x_distinct_s++;
        }
        #else
        // naive implementation with distance computation
        assert(false);
        #endif


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



void pso_constant_inertia_first_steps(struct pso_data_constant_inertia * pso)
{
    //TODO: step 1 and 2 "space-filling design"

    // 

    // Step 3. Initialize particle velocities
    for(int i = 0 ; i < pso->population_size ; i++)
    {
        for (int k=0 ; k < pso->dimensions ; k++)
        {
            double uk = rand_between(pso->bound_low[k], pso->bound_high[k]);
            PSO_V(pso, i)[k] = 1. / 2 * (uk - PSO_X(pso, 0, i)[k]);
        }
    }

    // Step 4. Initialise y, y_eval, and x_eval for each particle 
    for (int i = 0 ; i < pso->population_size ; i++)
    {
        for (int k = 0 ; k < pso->dimensions ;  k++)
        {
            PSO_Y(pso, i)[k] = PSO_X(pso, 0, i)[k];
        }

        double x_eval = pso->f(PSO_X(pso, 0, i), pso->dimensions);
        pso->y_eval[i] = x_eval;
        PSO_FX(pso, 0, i) = x_eval;

    }
    #ifdef DEBUG_STDOUT
        printf("before: find y\n");
    #endif
    // find ŷ
    
    double * y_hat = PSO_Y(pso, 0);
    double y_hat_eval = pso->y_eval[0];
    for (int i = 1 ; i < pso->population_size ; i++)
    {
        if (pso->y_eval[i] < y_hat_eval)
        {
            y_hat = PSO_Y(pso, i);
            y_hat_eval = pso->y_eval[i];
        }
    }
    pso->y_hat = y_hat;
    pso->y_hat_eval = y_hat_eval;
}


bool pso_constant_inertia_loop(struct pso_data_constant_inertia * pso)
{
    int const t = pso->time;

    // Step 5.
    // Fit surrogate
    // f already evaluated on x[0..t][0..i-1]
     if (fit_surrogate(pso) < 0)
     {
         fprintf(stderr, "ERROR: Failed to fit surrogate\n");
         exit(1);
     }


    #if LOG_SURROGATE
    {
        char fname[256] = {0};
        snprintf(fname, sizeof(fname), "surrogate_step5_t_%05d.struct", t);
        log_surrogate(fname, pso->lambda, pso->p, pso->x, t, pso->dimensions, pso->population_size);
    }
    #endif

    // Step 6
    // Determine new particle positions
    double * v_trial = malloc(pso->dimensions * sizeof(double));
    double * x_trial = malloc(pso->dimensions * sizeof(double));
    
    double * v_trial_best = malloc(pso->dimensions * sizeof(double));
    double * x_trial_best = malloc(pso->dimensions * sizeof(double));

    for(int i = 0 ; i < pso->population_size ; i++)
    {

        #if DEBUG_TRIALS
        printf("step6:start trials for particle %d\n", i);
        #endif

        double x_trial_best_seval = DBL_MAX;
        for (int l = 0 ; l < pso->n_trials ; l++)
        {
            
            for (int j = 0 ; j < pso->dimensions ; j++)
            {
                // compute v_i(t+1) from v_i(t)

                double w1 = (double)rand() / RAND_MAX;
                double w2 = (double)rand() / RAND_MAX;
                double v =
                    pso->inertia * PSO_V(pso, i)[j] +
                    pso->cognition * w1 * (PSO_Y(pso, i)[j] - PSO_X(pso, t, i)[j]) +
                    pso->social * w2 * (pso->y_hat[j] - PSO_X(pso, t, i)[j]);

                v_trial[j] = clamp(v, pso->vmin[j], pso->vmax[j]);


                x_trial[j] = clamp(PSO_X(pso, pso->time, i)[j] + v_trial[j], pso->bound_low[j], pso->bound_high[j]);
            }

            double x_trial_seval = surrogate_eval(pso, x_trial);

            #if DEBUG_TRIALS
            char trial_name[16] = {0};
            snprintf(trial_name, sizeof(trial_name), "x%d_trial%d", i, l);
            print_vectord(x_trial, pso->dimensions, trial_name);
            #ifdef DEBUG_STDOUT
                printf("score of trial %d = %f", l, x_trial_seval);
            #endif
            #endif

            if(x_trial_seval < x_trial_best_seval)
            {
                #if DEBUG_TRIALS
                printf(" (new best!)");
                #endif
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
            #if DEBUG_TRIALS
            printf("\n");
            #endif
        }

        // set next position and update velocity

        for (int j = 0 ; j < pso->dimensions ; j++)
        {
            PSO_X(pso, t+1, i)[j] = x_trial_best[j];
            PSO_V(pso, i)[j] = v_trial_best[j];
        }

        #if DEBUG_TRIALS
        printf("set x[t+1][i]:\n");
        char new_x_name[16] = {0};
        snprintf(new_x_name, sizeof(new_x_name), "x%d_(t=%d)", i, t+1);
        print_vectord(PSO_X(pso, t+1, i), pso->dimensions, new_x_name);
        #endif

    }

    free(x_trial);
    free(v_trial);
    free(x_trial_best);
    free(v_trial_best);



    // Step 7
    // Evaluate swarm positions
    for (int i = 0 ; i < pso->population_size ; i++)
    {
        PSO_FX(pso, t+1, i) = pso->f(PSO_X(pso, t+1, i), pso->dimensions);
    }

    // Step 8. Update the best positions per particle and overall

    for (int i = 0 ; i < pso->population_size ; i++)
    {
        double x_eval = PSO_FX(pso, t+1, i);

        if (x_eval < pso->y_eval[i])
        {
            // y_i <- x_i
            for (int k = 0 ; k < pso->dimensions ; k++)
            {
                PSO_Y(pso, i)[k] = PSO_X(pso, t+1, i)[k];
            }
            pso->y_eval[i] = x_eval;

            // is x_i(t+1) better than ŷ(t) ?
            if (x_eval < pso->y_hat_eval)
            {
                pso->y_hat = PSO_Y(pso, i);
            }
        }
    }

    pso->time++;

    // Step 9
    // Refit surrogate with time = t+1

    // first update the set of distinct points
    #if USE_ROUNDING_BLOOM_FILTER
    for(int i = 0 ; i < pso->population_size ; i++)
    {
        // add and check proximity to previous points
        if (! rounding_bloom_check_add(pso->bloom, pso->dimensions, PSO_X(pso, t+1, i), 1))
        {
            pso->x_distinct[pso->x_distinct_s] = (t+1) * pso->population_size + i;
            pso->x_distinct_s++;
        }
    }
    #else
    // naive implementation with distance computation
    assert(false);
    #endif


    if (fit_surrogate(pso) < 0)
    {
        fprintf(stderr, "ERROR: Failed to fit surrogate\n");
        exit(1);
    }

    #if LOG_SURROGATE
    {
        char fname[256] = {0};
        snprintf(fname, sizeof(fname), "surrogate_step9_t_%05d.struct", t);
        log_surrogate(fname, pso->lambda, pso->p, pso->x, t, pso->dimensions, pso->population_size);
    }
    #endif

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
    if(is_far_from_previous_evaluations(pso, x_local, pso->min_minimizer_distance))
    {
        double x_local_eval = pso->f(x_local, pso->dimensions);

        for (int k = 0 ; k < pso->dimensions ; k++)
        {
            PSO_PAST_REFINEMENT(pso, pso->n_past_refinement_points)[k] = x_local[k];
        }
        PSO_PAST_REFINEMENT_EVAL(pso, pso->n_past_refinement_points) = x_local_eval;

        // update overall best if applicable
        if(x_local_eval < pso->y_hat_eval)
        {
            pso->y_hat = PSO_PAST_REFINEMENT(pso, pso->n_past_refinement_points);
            pso->y_hat_eval = x_local_eval;
        }

        pso->n_past_refinement_points++;
    }

    free(x_local);




    return (pso->time < pso->time_max - 1);
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

    #ifdef DEBUG_STDOUT
    printf("t=%d  ŷ=[", pso.time);
    for (int j = 0 ; j < dimensions ; j++)
    {
        printf("%f", pso.y_hat[j]);
        if (j < dimensions - 1) printf(", ");
    }
    printf("]  f(ŷ)=%f\n", pso.y_hat_eval);
    #endif


    while(pso.time < pso.time_max - 1)
    {
        pso_constant_inertia_loop(&pso);

        #ifdef DEBUG_STDOUT
        printf("t=%d  ŷ=[", pso.time);
        for (int j = 0 ; j < dimensions ; j++)
        {
            printf("%f", pso.y_hat[j]);
            if (j < dimensions - 1) printf(", ");
        }
        printf("]  f(ŷ)=%f\n", pso.y_hat_eval);
        #endif
    }
}