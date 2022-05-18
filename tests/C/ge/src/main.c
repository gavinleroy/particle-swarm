#include "helpers.h"
#include "gaussian_elimination_solver.h"


#define DIM 3
int main(int argc, char ** argv)
{
    double A[DIM * DIM] = {1., 2., 3., 4., 5., 6, 7, 8, 10};
    double b[DIM] = {5., 6., 7.};

    print_matrixd(A, 3, "A");
    print_vectord(b, 3, "b");

    double Ab[DIM * (DIM+1)] = {0};
    for (int i = 0; i < DIM ; i++)
    {
        for(int j = 0 ; j < DIM; j++)
        {
            Ab[i * (DIM+1) + j] = A[i * DIM + j];
        }
        Ab[i * (DIM+1) + DIM] = b[i];
    }
    
    print_rect_matrixd(Ab, DIM, DIM+1, "Ab_before");


    double x[DIM] = {0}; 
    gaussian_elimination_solve(DIM, Ab, x);

    print_rect_matrixd(Ab, DIM, DIM+1, "Ab_after");
    print_vectord(x, DIM, "x");

}