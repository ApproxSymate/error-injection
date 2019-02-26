#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MonteCarlo.h"
#include "Random.h"
#include "bitflip.h"
#include <mpfr.h>

static const int SEED = 113;



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


#define MANTISSA_BITS 4
    double MonteCarlo_zero_masking(int Num_samples)
    {


        Random R = new_Random_seed(SEED);


        int under_curve = 0;
        int count;

        for (count=0; count<Num_samples; count++)
        {
            double x= Random_nextDouble(R);
            double y= Random_nextDouble(R);
            double x_approx = bit_mask(x, MANTISSA_BITS);
            double y_approx = bit_mask(y, MANTISSA_BITS);
            if ( x_approx*x_approx + y_approx*y_approx <= 1.0)
                 under_curve ++;

        }

        Random_delete(R);
        double result = ((double) under_curve / Num_samples) * 4.0;

        printf("zero under_curve %d, result %f \n", under_curve, result);
        result = bit_mask(result, MANTISSA_BITS);
        printf(" approx result %f \n", result);
        return result;
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

        mpfr_t temp_result;
        mpfr_init2(temp_result, FP_APPROX_FRACTION_BIT);
        double result = ((double) under_curve / Num_samples) * 4.0;
        mpfr_set_d (temp_result, result, MPFR_RNDZ);
        //printf("under_curve %d, result %f, result mpfr %f \n", under_curve, result, mpfr_get_d(temp_result,MPFR_RNDZ ));
        //under_curve 7011702, result 3.343440  ;
        //under_curve 6587369, result 3.141102  ;
        //under_curve 7011702, result 3.343440  ;
        return mpfr_get_d(temp_result,MPFR_RNDZ );
    }




    double kernel_measureMonteCarlo(  double min_time, Random R)
    {
		int N_samples = 10000;

    double result = 0.0;
    double res_mpfr;
    double res_bitmasking ;
    double res_ref;
    int cycles=1;
    int loops = 0;
    double mpfr_error = 0.0;
    cycles = 1492;


    res_mpfr = MonteCarlo_mpfr(cycles);
    res_bitmasking = MonteCarlo_zero_masking(cycles);
    res_ref= MonteCarlo(cycles);



      	 printf("MPFR %lf\n", fabs(res_mpfr - res_ref)/res_ref);

         printf("zero out mantissa bits %lf\n",fabs(res_bitmasking - res_ref)/res_ref);

		      printf("ref %.10f \n", res_ref);

        return 0.0;
    }

    int main(int argc, char *argv[])
    {

            Random R = new_Random_seed(SEED);


            kernel_measureMonteCarlo(0, R);


            return 0;

    }
