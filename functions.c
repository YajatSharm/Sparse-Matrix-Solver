#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include "functions.h"


void printMatrix(CSRMatrix *matrix)
{
    for (int i = 0; i < matrix->num_rows + 1; i++)
    {
        printf("%d ", matrix->row_ptr[i]);
    }
    printf("\n");
    for (int i = 0; i < matrix->num_non_zeros; i++)
    {
        printf("%d ", matrix->col_ind[i]);
    }
    printf("\n");

    for (int i = 0; i < matrix->num_non_zeros; i++)
    {
        printf("%f ", matrix->csr_data[i]);
    }
    printf("\n");
}


void clearMatrix(CSRMatrix *matrix)
{
    free(matrix->row_ptr);
    free(matrix->col_ind);
    free(matrix->csr_data);
}

void bubbleSort(int *row_ptr, int *col_ind, double *csr_data, int size)
{
    int i, j;
    int temp_row_ptr;
    int temp_col_ind;
    double temp_csr_data;
    for (i = 0; i < size - 1; i++)
    {
        for (j = 0; j < size - i - 1; j++)
        {
            if (row_ptr[j] > row_ptr[j + 1])
            {
                temp_row_ptr = row_ptr[j];
                temp_col_ind = col_ind[j];
                temp_csr_data = csr_data[j];
                row_ptr[j] = row_ptr[j + 1];
                row_ptr[j + 1] = temp_row_ptr;
                col_ind[j] = col_ind[j + 1];
                col_ind[j + 1] = temp_col_ind;
                csr_data[j] = csr_data[j + 1];
                csr_data[j + 1] = temp_csr_data;
            }
        }
    }
}

void ReadMMtoCSR(const char *filename, CSRMatrix *matrix)
{
    //matrix market file to read data from
    FILE *MM = fopen(filename, "r");
    //file to store coordinates of non-zero entries to be used for graphing
    FILE *storage  = fopen("Matrix_data.txt", "w");

    char buffer[100];
    if (MM == NULL)
    {
        printf("The Matrix Market file could not be opened :(\n");
    }
    else
    {
        fgets(buffer, 100, MM);
        bool isSymmetric = false;
        
        //check if the first line of the file contains the word "symmetric"
        if (strstr(buffer, "symmetric"))
        {
            isSymmetric = true;
        }

        //loop through all the comments
        do
        {
            fgets(buffer, sizeof(buffer), MM);

        } while (buffer[0] == '%');
        int temp_index = 0;
        /*
        find the length of the last line that was read
        This is to prevent a bug from affecting the output because the previous 
        loop ends up reading the first non-comment line as well
        */
        while (buffer[temp_index] != '\0')
        {
            temp_index++;
        }

        //move the pointer back the length of the first non-comment line which contains
        //The number of rows, columns, and non-zeros (not the matrix, this will be changed because the matrix is symmetric) in the file
        fseek(MM, -temp_index, SEEK_CUR);
        fscanf(MM, "%d %d %d", &matrix->num_rows, &matrix->num_cols, &matrix->num_non_zeros);

        int row, col;
        if (isSymmetric)
        {
            //Allocate the necessary amount of memory
            matrix->num_non_zeros = (2 * matrix->num_non_zeros) - matrix->num_cols; //Assume that the diagonal entries are populated
            matrix->row_ptr = (int *)malloc(matrix->num_non_zeros * sizeof(int));
            matrix->col_ind = (int *)malloc(matrix->num_non_zeros * sizeof(int));
            matrix->csr_data = (double *)malloc(matrix->num_non_zeros * sizeof(double));
            for (int i = 0; i < matrix->num_non_zeros; i++)
            {
                //get the data from file and ensure that the row and column are 1 less than what's in the file 
                //because indexing starts from 0 in C
                fscanf(MM, "%d %d %lf", &row, &col, &matrix->csr_data[i]);
                matrix->row_ptr[i] = row - 1;
                matrix->col_ind[i] = col - 1;

                //only duplicate entries that are not on the diagonal
                if (matrix->row_ptr[i] != matrix->col_ind[i])
                {
                    //swap rows and columns
                    matrix->row_ptr[i + 1] = matrix->col_ind[i];
                    matrix->col_ind[i + 1] = matrix->row_ptr[i];
                    matrix->csr_data[i + 1] = matrix->csr_data[i];
                    i++; //skip the next index in the array because it is already populated now
                }
            }
        }
        else
        {
            //less memory is needed if the matrix is not symmetric
            matrix->row_ptr = (int *)malloc(matrix->num_non_zeros * sizeof(int));
            matrix->col_ind = (int *)malloc(matrix->num_non_zeros * sizeof(int));
            matrix->csr_data = (double *)malloc(matrix->num_non_zeros * sizeof(double));
            for (int i = 0; i < matrix->num_non_zeros; i++)
            {
                //Similar to the symmetric case but no need to duplicate entries
                fscanf(MM, "%d %d %lf", &row, &col, &matrix->csr_data[i]);
                matrix->row_ptr[i] = row - 1;
                matrix->col_ind[i] = col - 1;
            }
        }

        //sort according to row index in an ascending order
        bubbleSort(matrix->row_ptr, matrix->col_ind, matrix->csr_data, matrix->num_non_zeros);
    
        fprintf(storage, "%d %d ", matrix->num_rows, matrix->num_cols);
        fprintf(storage, "\n");
        int row_ptr_index = 1;
        for (int i = 0; i <= matrix->num_non_zeros; i++)
        {
            //print the row and column indexes of non-zero entries into the storage file 
            //This will be used for graphing later on
            fprintf(storage, "%d %d ", matrix->row_ptr[i], matrix->col_ind[i]);
            fprintf(storage, "\n");

            //change row_ptr array into a compressed row format
            if (matrix->row_ptr[i] != matrix->row_ptr[i + 1])
            {
                matrix->row_ptr[row_ptr_index] = i + 1;
                row_ptr_index++;
            } else {

            }
        }

        //get rid of garbage values that and only keep row pointers
        matrix->row_ptr = (int *)realloc(matrix->row_ptr, (row_ptr_index + 1) * sizeof(int));
    }
    //close all files
    fclose(MM);
    fclose(storage);
}

void spmv_csr(const CSRMatrix *A, const double *x, double *y)
{
    double sum = 0.0;
    
    //Sparse matrix multiplication algorithm that only deals with non-zero entries in the matrix
    for (int i = 0; i < A->num_rows; i++)
    {
        sum = 0.0;
        for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++)
        {
            sum += (A->csr_data[j] * x[A->col_ind[j]]);
        }
        y[i] = sum;
    }
}

void solve(const CSRMatrix *A, const double *b, double *x)
{
    //declare variables for finding CPU time 
    clock_t start_time;
    clock_t end_time;
    double cpu_time_used;
    start_time = clock();

    //implementation of Jacobi algorithm
    double tol = 1e-16;
    double numerator_sum = 0.0;
    double denominator = 0.0;
    double max_diff = 0.0;
    double *x_old = (double *)malloc(A->num_rows * sizeof(double));
    for(int i = 0; i < A->num_rows; i++) {
        x_old[i] = 0.0;
    }
    for (int i = 0; i < 1000000; i++)
    {
        max_diff = 0.0;
        for (int row_index = 0; row_index < A->num_rows; ++row_index)
        {
            numerator_sum = 0.0;
        
            x_old[row_index] = x[row_index]; //save the x value from the previous iteration to be used when finding the difference later
        
            for (int col_index = A->row_ptr[row_index]; col_index < A->row_ptr[row_index + 1]; ++col_index)
            {
                if (row_index != A->col_ind[col_index])
                {
                    numerator_sum += A->csr_data[col_index] * x[A->col_ind[col_index]];
                }
                else
                {
                    denominator = (A->csr_data[col_index] + 1e-16); //prevent division by 0 
                }
            }
            x[row_index] = (b[row_index] - numerator_sum) / denominator;
            //Check whether the new difference is greater than the previous one
            max_diff = fmax(max_diff, fabs(x[row_index] - x_old[row_index]));
        }
        //Compare the maximum differene value to the tolerance
        //This will ensure that the algorithm produces the most accurate result possible
        if (max_diff < tol)
        {
            break;
        } 
    }
    free(x_old);

    // Print the time taken to solve Ax=b
    end_time = clock();
    cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("CPU time taken to solve Ax=b: %f seconds\n", cpu_time_used);
}

void compute_residual(const CSRMatrix *A, const double *b, const double *x, double *r)
{
    double *mult_result = (double *)malloc(A->num_rows * sizeof(double));
    for(int i = 0; i < A->num_rows; i++) {
        mult_result[i] = 0.0;
    }
    
    //find the result of matrix multiplication with the x-values that have been found from the solve method and save it to mult_result
    spmv_csr(A, x, mult_result);
    for (int i = 0; i < A->num_rows; i++)
    {   
        r[i] = mult_result[i] - b[i];
    }
    free(mult_result);
}

double compute_norm(double *r, int size)
{
    double sum = 0.0;
    //find the norm of the residual vector.
    for (int i = 0; i < size; i++)
    {
        sum += r[i] * r[i];
    }
    return sqrt(sum);
}


