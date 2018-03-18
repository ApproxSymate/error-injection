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

#define FP_APPROX_FRACTION_BIT 7
#define FP_FULL_FRACTION_BIT 53

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
            mpfr_set_d(x_temp,  x, MPFR_RNDZ);
            mpfr_set_d(y_temp, y, MPFR_RNDZ);			
			mpfr_mul(temp_variable1, x_temp, x_temp, MPFR_RNDZ);
			mpfr_mul(temp_variable2, y_temp, y_temp, MPFR_RNDZ);
			mpfr_add(temp_variable1,temp_variable1,temp_variable2,MPFR_RNDZ);
			if(mpfr_get_d(temp_variable1, MPFR_RNDZ) <= 1.0)
					under_curve ++ ;
           // if ( x*x + y*y <= 1.0)
            //     under_curve ++;
            
        }

        Random_delete(R);

        return ((double) under_curve / Num_samples) * 4.0;
    }


    double MonteCarlo_bitflip(int Num_samples)
    {


        Random R = new_Random_seed(SEED);


        int under_curve = 0;
        int count;

        for (count=0; count<Num_samples; count++)
        {
            double x= bitflip_float(Random_nextDouble(R));
            double y= bitflip_float(Random_nextDouble(R));
		
            if ( x*x + y*y <= 1.0)
                 under_curve ++;
            
        }

        Random_delete(R);

        return ((double) under_curve / Num_samples) * 4.0;
    }

    double kernel_measureMonteCarlo(  double min_time, Random R)
    {
		int N_samples = 10000;
   
        double res_mpfr = MonteCarlo_mpfr(N_samples);

        double res_bitflip = MonteCarlo_bitflip(N_samples);

        double res_ref= MonteCarlo(N_samples);

      	printf("MPFR %lf\n",fabs(res_mpfr - res_ref)/res_ref);

        printf("bitflip %lf\n",fabs(res_bitflip - res_ref)/res_ref);

		printf("ref %lf \n", res_ref);
		
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
