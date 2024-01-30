#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

int main(int argc, char *argv[])
{

    // <Handle the inputs here>
    if(argc != 2) {
        printf("Please enter one file name.");
        return 1;
    }
    const char *filename = argv[1];
    CSRMatrix A;
    ReadMMtoCSR(filename, &A);

    // Initializing all the vector b (in Ax=b)
    double *b = (double *)malloc(A.num_rows * sizeof(double));
    // Set all elements of b to 1
    for (int i = 0; i < A.num_rows; ++i)
    {
        b[i] = 1.0;
    }

    // <The rest of your code goes here>
    printf("The matrix name: %s\n", filename);
    printf("The dimension of the matrix: %d by %d\n", A.num_rows, A.num_cols);
    printf("Number of non-zeros: %d\n", A.num_non_zeros);
    
    //allocated necessary memory
    double *x = (double *)malloc(A.num_cols * sizeof(double));
    double *residual = (double *)malloc(A.num_cols * sizeof(double));

    //initial all x-values to 0.0 (This will serve as the initial guess for the Jacobi method algorithm)
    for (int i = 0; i < A.num_rows; ++i)
    {
        x[i] = 1.0;
        residual[i] = 0;
    }
    solve(&A, b, x);
    compute_residual(&A, b, x, residual);
    printf("Residual Norm: %.2e\n", compute_norm(residual, A.num_rows));
   
    //free all the allocated memory
    free(x);
    free(b);
    free(residual);
    clearMatrix(&A);
}