#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LU.h"
#include "FFT.h"
#include "SOR.h"
#include "MonteCarlo.h"
#include "LU.h"
#include "Random.h"
#include "Stopwatch.h"
#include "SparseCompRow.h"
#include "array.h"
#include "constants.h"
#include "bitflip.h"
#include <mpfr.h>

    double SparseCompRow_num_flops(int N, int nz, int num_iterations)
    {
        /* Note that if nz does not divide N evenly, then the
           actual number of nonzeros used is adjusted slightly.
        */
        int actual_nz = (nz/N) * N;
        return ((double)actual_nz) * 2.0 * ((double) num_iterations);
    }


    void SparseCompRow_matmult( int M, double *y, double *val, int *row,
        int *col, double *x, int NUM_ITERATIONS)
    {
        int reps;
        int r;
        int i;

        for (reps=0; reps<NUM_ITERATIONS; reps++)
        {

            for (r=0; r<M; r++)
            {
                double sum = 0.0;
                int rowR = row[r];
                int rowRp1 = row[r+1];
                for (i=rowR; i<rowRp1; i++)
                    sum += x[ col[i] ] * val[i];
                y[r] = sum;
            }
        }


    }

    void SparseCompRow_matmult_approx_bitflip( int M, double *y, double *val, int *row,
        int *col, double *x, int NUM_ITERATIONS)
    {
        int reps;
        int r;
        int i;

        for (reps=0; reps<NUM_ITERATIONS; reps++)
        {

            for (r=0; r<M; r++)
            {
                double sum = 0.0;
                int rowR = row[r];
                int rowRp1 = row[r+1];
                for (i=rowR; i<rowRp1; i++){
                    double temp_x = bitflip_float(x[ col[i] ]);
                    double temp_val = bitflip_float(val[i]);
                  //  printf(" %lf %lf  \n", temp_x , bitflip_float(temp_x));
#ifdef AGGRESSIVE
                    sum+= bitflip_float(temp_x * temp_val);
#else
                    sum+= (temp_x * temp_val);
#endif
                    //sum += x[ col[i] ] * val[i];
                  }
#ifdef AGGRESSIVE
                y[r] =bitflip_float(sum);
#else
                y[r] =(sum);
#endif
            }
        }


    }

    void SparseCompRow_matmult_approx_mpfr( int M, double *y, double *val, int *row,
        int *col, double *x, int NUM_ITERATIONS)
    {


        int reps;
        int r;
        int i;
        mpfr_t x_temp;
        mpfr_t val_temp;
        mpfr_t temp_variable;
        mpfr_init2(x_temp, FP_APPROX_FRACTION_BIT);
        mpfr_init2(val_temp, FP_APPROX_FRACTION_BIT);
        mpfr_init2(temp_variable, FP_APPROX_FRACTION_BIT);
        for (reps=0; reps<NUM_ITERATIONS; reps++)
        {

            for (r=0; r<M; r++)
            {
                double sum = 0.0;
                int rowR = row[r];
                int rowRp1 = row[r+1];
                for (i=rowR; i<rowRp1; i++){
                    double temp_x = x[ col[i] ];
                    double temp_val = val[i];
                    mpfr_set_d(x_temp,  x[ col[i] ], MPFR_RNDZ);
                    mpfr_set_d(val_temp, val[i], MPFR_RNDZ);
                    mpfr_mul(temp_variable, val_temp, x_temp, MPFR_RNDZ);
                    sum+= mpfr_get_d(temp_variable, MPFR_RNDZ);
                    //sum += x[ col[i] ] * val[i];
                  }
                y[r] = sum;
            }
        }


    }

    double kernel_measureSparseMatMult(int N, int nz,
            double min_time, Random R)
    {
        /* initialize vector multipliers and storage for result */
        /* y = A*y;  */

        double *x = RandomVector(N, R);
        double *y = (double*) malloc(sizeof(double)*N);
        double *y_error_mpfr = (double*) malloc(sizeof(double)*N);
        double *y_error_bitflip = (double*) malloc(sizeof(double)*N);
        double result = 0.0;


        // initialize square sparse matrix
        //
        // for this test, we create a sparse matrix with M/nz nonzeros
        // per row, with spaced-out evenly between the begining of the
        // row to the main diagonal.  Thus, the resulting pattern looks
        // like
        //             +-----------------+
        //             +*                +
        //             +***              +
        //             +* * *            +
        //             +** *  *          +
        //             +**  *   *        +
        //             +* *   *   *      +
        //             +*  *   *    *    +
        //             +*   *    *    *  +
        //             +-----------------+
        //
        // (as best reproducible with integer artihmetic)
        // Note that the first nr rows will have elements past
        // the diagonal.


        int nr = nz/N;      /* average number of nonzeros per row  */
        int anz = nr *N;    /* _actual_ number of nonzeros         */


        double *val = RandomVector(anz, R);
        int *col = (int*) malloc(sizeof(int)*nz);
        int *row = (int*) malloc(sizeof(int)*(N+1));
        int r=0;
        int cycles=1;

        Stopwatch Q = new_Stopwatch();

        row[0] = 0;
        for (r=0; r<N; r++)
        {
            /* initialize elements for row r */

            int rowr = row[r];
            int step = r/ nr;
            int i=0;

            row[r+1] = rowr + nr;
            if (step < 1) step = 1;   /* take at least unit steps */


            for (i=0; i<nr; i++)
                col[rowr+i] = i*step;

        }



        SparseCompRow_matmult_approx_mpfr(N, y_error_mpfr, val, row, col, x, cycles);

        SparseCompRow_matmult_approx_bitflip(N, y_error_bitflip, val, row, col, x, cycles);

        SparseCompRow_matmult(N, y, val, row, col, x, cycles);

        double sum_err =0 ,err;
        int i;
        for(i=0;i<N;i++)
      	{
      		err = y_error_mpfr[i] - y[i];
      		sum_err += (err*err);
      	}
      	printf("MPFR %lf\n",sqrt(sum_err/N));

        sum_err =0.0;
        for(i=0;i<N;i++)
        {
          err = y_error_bitflip[i] - y[i];
          sum_err += (err*err);
        }

        printf("bitflip %lf\n",sqrt(sum_err/N));

        Stopwatch_delete(Q);
        free(row);
        free(col);
        free(val);
        free(y);
        free(x);

        return result;
    }

    int main(int argc, char *argv[])
    {
            /* default to the (small) cache-contained version */

            double min_time = RESOLUTION_DEFAULT;

            int FFT_size = FFT_SIZE;
            int SOR_size =  SOR_SIZE;
            int Sparse_size_M = SPARSE_SIZE_M;
            int Sparse_size_nz = SPARSE_SIZE_nz;
            int LU_size = LU_SIZE;


            /* run the benchmark */

            double res[6] = {0.0};
            Random R = new_Random_seed(RANDOM_SEED);


            if (argc > 1)
            {
    			int current_arg = 1;

    			if (strcmp(argv[1], "-help")==0  ||
    					strcmp(argv[1], "-h") == 0)
    			{
    				fprintf(stderr, "Usage: [-large] [minimum_time]\n");
    				exit(0);
    			}

    			if (strcmp(argv[1], "-large")==0)
    			{
    				FFT_size = LG_FFT_SIZE;
    				SOR_size = LG_SOR_SIZE;
    				Sparse_size_M = LG_SPARSE_SIZE_M;
    				Sparse_size_nz = LG_SPARSE_SIZE_nz;
    				LU_size = LG_LU_SIZE;

    				current_arg++;
    			}

    			if (current_arg < argc)
    			{
    				min_time = atof(argv[current_arg]);
    			}

            }


  //  		printf("Using %10.2f seconds min time per kenel.\n", min_time);
/*
            res[1] = kernel_measureFFT( FFT_size, min_time, R);
            res[2] = kernel_measureSOR( SOR_size, min_time, R);
            res[3] = kernel_measureMonteCarlo(min_time, R);
  */
            res[4] = kernel_measureSparseMatMult( Sparse_size_M,
                    Sparse_size_nz, min_time, R);
          /*
            res[5] = kernel_measureLU( LU_size, min_time, R);

            res[0] = (res[1] + res[2] + res[3] + res[4] + res[5]) / 5;


            printf("Composite Score:        %8.2f\n" ,res[0]);
            printf("FFT             Mflops: %8.2f    (N=%d)\n", res[1], FFT_size);
            printf("SOR             Mflops: %8.2f    (%d x %d)\n",
    				res[2], SOR_size, SOR_size);
            printf("MonteCarlo:     Mflops: %8.2f\n", res[3]);
            printf("Sparse matmult  Mflops: %8.2f    (N=%d, nz=%d)\n", res[4],
    					Sparse_size_M, Sparse_size_nz);
            printf("LU              Mflops: %8.2f    (M=%d, N=%d)\n", res[5],
    				LU_size, LU_size);
*/

            Random_delete(R);

            return 0;

    }
