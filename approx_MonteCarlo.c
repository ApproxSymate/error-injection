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

static const int SEED = 113;


    double MonteCarlo_num_flops(int Num_samples)
    {
        /* 3 flops in x^2+y^2 and 1 flop in random routine */

        return ((double) Num_samples)* 4.0;

    }



    double MonteCarlo(int Num_samples)
    {


        Random R = new_Random_seed(SEED);


        int under_curve = 0;
        int count;

        for (count=0; count<Num_samples; count++)
        {
            double x= Random_nextDouble(R);
            double y= Random_nextDouble(R);

            if ( x*x + y*y <= 1.0)
                 under_curve ++;

        }

        Random_delete(R);
        double result = ((double) under_curve / Num_samples) * 4.0;

        printf("under_curve %d, result %f \n", under_curve, result);
        return ((double) under_curve / Num_samples) * 4.0;
    }


    double MonteCarlo_mpfr(int Num_samples)
    {


        Random R = new_Random_seed(SEED);

		mpfr_t x_temp;
        mpfr_t y_temp;
        mpfr_t temp_variable1;
        mpfr_t temp_variable2;
        mpfr_init2(x_temp, FP_APPROX_FRACTION_BIT);
        mpfr_init2(y_temp, FP_APPROX_FRACTION_BIT);
		mpfr_init2(temp_variable1, FP_FULL_FRACTION_BIT);
		mpfr_init2(temp_variable2, FP_FULL_FRACTION_BIT);

        int under_curve = 0;
        int count;

        for (count=0; count<Num_samples; count++)
        {
            double x= Random_nextDouble(R);
            double y= Random_nextDouble(R);
            mpfr_set_d(x_temp,  x, MPFR_RNDN);
            mpfr_set_d(y_temp, y, MPFR_RNDN);
			mpfr_mul(temp_variable1, x_temp, x_temp, MPFR_RNDN);
			mpfr_mul(temp_variable2, y_temp, y_temp, MPFR_RNDN);
			mpfr_add(temp_variable1,temp_variable1,temp_variable2,MPFR_RNDN);
			if(mpfr_get_d(temp_variable1, MPFR_RNDN) <= 1.0)
					under_curve ++ ;
           // if ( x*x + y*y <= 1.0)
            //     under_curve ++;

        }

        Random_delete(R);

        mpfr_t temp_result;
        mpfr_init2(temp_result, FP_APPROX_FRACTION_BIT);
        double result = ((double) under_curve / Num_samples) * 4.0;
        mpfr_set_d (temp_result, result, MPFR_RNDN);
        printf("under_curve %d, result %f, result mpfr %f \n", under_curve, result, mpfr_get_d(temp_result,MPFR_RNDN ));
        //under_curve 7011702, result 3.343440  ;
        //under_curve 6587369, result 3.141102  ;
        //under_curve 7011702, result 3.343440  ;
        return mpfr_get_d(temp_result,MPFR_RNDN );
    }


    double MonteCarlo_bitflip(int Num_samples)
    {


        Random R = new_Random_seed(SEED);


        int under_curve = 0;
        int count;

        for (count=0; count<Num_samples; count++)
        {
            double x= bitflip_float(Random_nextDouble(R));
#ifdef AGGRESSIVE
            double y= bitflip_float(Random_nextDouble(R));
#else
            double y= (Random_nextDouble(R));
#endif
            if ( x*x + y*y <= 1.0)
                 under_curve ++;

        }

        Random_delete(R);

        return ((double) under_curve / Num_samples) * 4.0;
    }

    double kernel_measureMonteCarlo(  double min_time, Random R)
    {
		int N_samples = 10000;

    double result = 0.0;
    Stopwatch Q = new_Stopwatch();
double res_mpfr;
double res_bitflip ;
double res_ref;
    int cycles=1;
    int loops = 0;
    double mpfr_error = 0.0;
    cycles = 1492;
    while(1)
    {
        Stopwatch_start(Q);
      //  MonteCarlo_integrate(cycles);
         res_mpfr = MonteCarlo_mpfr(cycles);

         res_bitflip = MonteCarlo_bitflip(cycles);

         res_ref= MonteCarlo(cycles);
         //printf("cycles %d", cycles);

        Stopwatch_stop(Q);
        /*
        if (cycles > 256 )
          break;
        if (Stopwatch_read(Q) >= min_time) break;
*/
        break;
        cycles *= 2;
        //mpfr_error += fabs(res_mpfr - res_ref)/res_ref ;
        loops ++;
    }




      	printf("MPFR %lf\n", fabs(res_mpfr - res_ref)/res_ref);

        printf("bitflip %lf\n",fabs(res_bitflip - res_ref)/res_ref);

		printf("ref %.10f \n", res_ref);

        return 0.0;
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

            //~ res[1] = kernel_measureFFT( FFT_size, min_time, R);
            //~ res[2] = kernel_measureSOR( SOR_size, min_time, R);
            res[3] = kernel_measureMonteCarlo(min_time, R);

            //~ res[4] = kernel_measureSparseMatMult( Sparse_size_M,
                    //~ Sparse_size_nz, min_time, R);


            //~ res[5] = kernel_measureLU( LU_size, min_time, R);

            //~ res[0] = (res[1] + res[2] + res[3] + res[4] + res[5]) / 5;


            //~ printf("Composite Score:        %8.2f\n" ,res[0]);
            //~ printf("FFT             Mflops: %8.2f    (N=%d)\n", res[1], FFT_size);
            //~ printf("SOR             Mflops: %8.2f    (%d x %d)\n",
    				//~ res[2], SOR_size, SOR_size);
            //~ printf("MonteCarlo:     Mflops: %8.2f\n", res[3]);
            //~ printf("Sparse matmult  Mflops: %8.2f    (N=%d, nz=%d)\n", res[4],
    					//~ Sparse_size_M, Sparse_size_nz);
            //~ printf("LU              Mflops: %8.2f    (M=%d, N=%d)\n", res[5],
    				//~ LU_size, LU_size);


            Random_delete(R);

            return 0;

    }
