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



#define LEN 7
 int config_vals[LEN];
  mpfr_t omega_over_four_SOR_execute;
  mpfr_t omega_SOR_execute;
  mpfr_t one_minus_omega_SOR_execute;
mpfr_t Gi_SOR_execute_mpfr;
mpfr_t Gim1_SOR_execute_mpfr;
mpfr_t Gip1_SOR_execute_mpfr;
    mpfr_t temp_var_1;
    mpfr_t temp_var_2;
            mpfr_t temp_var_3;
            mpfr_t temp_var_4;
            mpfr_t temp_var_5;
            mpfr_t temp_var_6;
            mpfr_t temp_var_7;
int init_mpfr() {
  mpfr_init2(omega_over_four_SOR_execute, config_vals[2]);
mpfr_init2(omega_SOR_execute, config_vals[1]);
  mpfr_init2(one_minus_omega_SOR_execute, config_vals[3]);
  mpfr_init2(Gi_SOR_execute_mpfr, config_vals[4]);
  mpfr_init2(Gim1_SOR_execute_mpfr, config_vals[5]);
  mpfr_init2(Gip1_SOR_execute_mpfr, config_vals[6]);
    mpfr_init2 (temp_var_1, config_vals[0]);
    mpfr_init2 (temp_var_2, config_vals[0]);
            mpfr_init2 (temp_var_3, config_vals[0]);
            mpfr_init2 (temp_var_4, config_vals[0]);
            mpfr_init2 (temp_var_5, config_vals[0]);
            mpfr_init2 (temp_var_6, config_vals[0]);
            mpfr_init2 (temp_var_7, config_vals[0]);
}

int init_readconfig() {

  // For reading precision contents of config_file into array
   FILE *myFile;
     myFile = fopen("SOR_config_file.txt", "r");

        if (myFile == NULL) {
					printf("Error Reading File\n");
                exit (0);
                }

        int s;
        for (s = 0; s < LEN; s++) {
            fscanf(myFile, "%d,", &config_vals[s]);
                          }

        fclose(myFile);
        init_mpfr();
        return 0;
}




    double SOR_num_flops(int M, int N, int num_iterations)
    {
        double Md = (double) M;
        double Nd = (double) N;
        double num_iterD = (double) num_iterations;

        return (Md-1)*(Nd-1)*num_iterD*6.0;
    }

    void SOR_execute(int M, int N, double omega, double **G, int
            num_iterations)
    {

        double omega_over_four = omega * 0.25;
        double one_minus_omega = 1.0 - omega;

        /* update interior points */

        int Mm1 = M-1;
        int Nm1 = N-1;
        int p;
        int i;
        int j;
        double *Gi;
        double *Gim1;
        double *Gip1;

        for (p=0; p<num_iterations; p++)
        {
            for (i=1; i<Mm1; i++)
            {
                Gi = G[i];
                Gim1 = G[i-1];
                Gip1 = G[i+1];
                for (j=1; j<Nm1; j++)
                    Gi[j] = omega_over_four * (Gim1[j] + Gip1[j] + Gi[j-1]
                                + Gi[j+1]) + one_minus_omega * Gi[j];
            }
        }
    }

   void SOR_execute_mpfr(int M, int N, double omega, double **G, int
            num_iterations)
    {
      mpfr_set_d(omega_SOR_execute,omega,MPFR_RNDZ);
    mpfr_mul_d(omega_over_four_SOR_execute, omega_SOR_execute, 0.25, MPFR_RNDZ);
      mpfr_d_sub(one_minus_omega_SOR_execute, 1.0, omega_SOR_execute, MPFR_RNDZ);
  int Mm1 = M - 1;
  int Nm1 = N - 1;
  int p;
  int i;
  int j;
  double *Gi;
  double *Gim1;
  double *Gip1;
  for (p = 0; p < num_iterations; p++){
  {
    for (i = 1; i < Mm1; i++){
    {

		                Gi = G[i];
                Gim1 = G[i-1];
                Gip1 = G[i+1];

      for (j = 1; j < Nm1; j++){


 mpfr_set_d(Gip1_SOR_execute_mpfr,Gip1[j],MPFR_RNDZ);
 mpfr_set_d(Gim1_SOR_execute_mpfr,Gim1[j],MPFR_RNDZ);

            mpfr_add(temp_var_3, Gim1_SOR_execute_mpfr, Gip1_SOR_execute_mpfr, MPFR_RNDZ);
 mpfr_set_d(Gi_SOR_execute_mpfr,Gi[j - 1],MPFR_RNDZ);

            mpfr_add(temp_var_4, (temp_var_3), Gi_SOR_execute_mpfr, MPFR_RNDZ);
 mpfr_set_d(Gi_SOR_execute_mpfr,Gi[j + 1],MPFR_RNDZ);

            mpfr_add(temp_var_5, (temp_var_4), Gi_SOR_execute_mpfr, MPFR_RNDZ);

            mpfr_mul(temp_var_6, omega_over_four_SOR_execute, (temp_var_5), MPFR_RNDZ);
 mpfr_set_d(Gi_SOR_execute_mpfr,Gi[j],MPFR_RNDZ);

            mpfr_mul(temp_var_7, one_minus_omega_SOR_execute, Gi_SOR_execute_mpfr, MPFR_RNDZ);
            mpfr_add(Gi_SOR_execute_mpfr, (temp_var_6), (temp_var_7), MPFR_RNDZ);
      Gi[j] = mpfr_get_d(Gi_SOR_execute_mpfr, MPFR_RNDZ);

            }

    }

        }

  }

    }

}


    void SOR_execute_bitflip(int M, int N, double omega, double **G, int
            num_iterations)
    {

		double omega_temp = bitflip_float(omega);
#ifdef AGGRESSIVE
        double omega_over_four = bitflip_float(omega_temp * 0.25);
        double one_minus_omega = bitflip_float(1.0 - omega_temp);
#else
        double omega_over_four = (omega_temp * 0.25);
        double one_minus_omega = (1.0 - omega_temp);
#endif

        /* update interior points */

        int Mm1 = M-1;
        int Nm1 = N-1;
        int p;
        int i;
        int j;
        double *Gi;
        double *Gim1;
        double *Gip1;

        for (p=0; p<num_iterations; p++)
        {
            for (i=1; i<Mm1; i++)
            {
                Gi = G[i];
                Gim1 = G[i-1];
                Gip1 = G[i+1];
                for (j=1; j<Nm1; j++)
#ifdef AGGRESSIVE
                    Gi[j] =bitflip_float(omega_over_four * (Gim1[j] + Gip1[j] + Gi[j-1]
                                + Gi[j+1]) + one_minus_omega * Gi[j]);
#else
                  Gi[j] =(omega_over_four * (Gim1[j] + Gip1[j] + Gi[j-1]
                      +  Gi[j+1]) + one_minus_omega * Gi[j]);
#endif
            }
        }
    }


  int main ()
    {
	int N = SOR_SIZE,i,j;

 init_readconfig();

        Random R = new_Random_seed(RANDOM_SEED);
        Random R_mpfr= new_Random_seed(RANDOM_SEED);
        Random R_bitflip = new_Random_seed(RANDOM_SEED);

	double min_time = RESOLUTION_DEFAULT;

        double **G = RandomMatrix(N, N, R);
        double **G_mpfr= RandomMatrix(N, N, R_mpfr);
        double **G_bitflip= RandomMatrix(N, N, R_bitflip);

        double result = 0.0;

        Stopwatch Q = new_Stopwatch();
        int cycles=1;
       while(1)
        {
            Stopwatch_start(Q);
            SOR_execute(N, N, 1.25, G, cycles);
            SOR_execute_mpfr(N, N, 1.25, G_mpfr, cycles);
//            SOR_execute_bitflip(N, N, 1.25, G_bitflip, cycles);
            Stopwatch_stop(Q);

            if (Stopwatch_read(Q) >= min_time) break;

            cycles*=2;
        }

             //~ SOR_execute(N, N, 1.25, G, cycles);
            //~ SOR_execute_mpfr(N, N, 1.25, G_mpfr, cycles);
            //~ SOR_execute_bitflip(N, N, 1.25, G_bitflip, cycles);

         double sum = 0;
		double err =0;
        for(i=0;i<N;i++)
                for(j=0;j<N;j++){
                        err = (G[i][j]-G_mpfr[i][j])*(G[i][j]-G_mpfr[i][j]);
                        if(err==err && err >= 0){
							sum = sum + err;

                        }
                }

        printf("%lf",sqrt(sum/(N*N)));

    }
