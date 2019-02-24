#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LU.h"
#include "approx_fft.h"
#include "bitflip.h"
#include "MonteCarlo.h"
#include "LU.h"
#include "Random.h"
#include "Stopwatch.h"
#include "SparseCompRow.h"
#include "array.h"
#include "constants.h"
#include <mpfr.h>


double PI;//  3.1415926535897932

/*-----------------------------------------------------------------------*/

static int int_log2(int n);

double FFT_num_flops(int N)
{

     double Nd = (double) N;
     double logN = (double) int_log2(N);

     return (5.0*Nd-2)*logN + 2*(Nd+1);
}

static int int_log2 (int n)
{
    int k = 1;
    int log = 0;
    for(/*k=1*/; k < n; k *= 2, log++);
    if (n != (1 << log))
    {
      printf("FFT: Data length is not a power of 2!: %d ",n);
      exit(1);
    }
    return log;
}

static void FFT_transform_internal (int N, double *data, int direction) {
    int n = N/2;
    int bit = 0;
    int logn;
    int dual = 1;

    if (n == 1) return;         /* Identity operation! */
    logn = int_log2(n);


    if (N == 0) return;

    /* bit reverse the input data for decimation in time algorithm */
    FFT_bitreverse(N, data) ;

    /* apply fft recursion */
    /* this loop executed int_log2(N) times */
    for (bit = 0; bit < logn; bit++, dual *= 2) {
      double w_real = 1.0;
      double w_imag = 0.0;
      int a;
      int b;

      double theta = 2.0 * direction * PI / (2.0 * (double) dual);
      double s = sin(theta);
      double t = sin(theta / 2.0);
      double s2 = 2.0 * t * t;

      for (a=0, b = 0; b < n; b += 2 * dual) {
        int i = 2*b ;
        int j = 2*(b + dual);

        double wd_real = data[j] ;
        double wd_imag = data[j+1] ;

        data[j]   = data[i]   - wd_real;
        data[j+1] = data[i+1] - wd_imag;
        data[i]  += wd_real;
        data[i+1]+= wd_imag;
      }

      /* a = 1 .. (dual-1) */
      for (a = 1; a < dual; a++) {
        /* trignometric recurrence for w-> exp(i theta) w */
        {
          double tmp_real = w_real - s * w_imag - s2 * w_real;
          double tmp_imag = w_imag + s * w_real - s2 * w_imag;
          w_real = tmp_real;
          w_imag = tmp_imag;
        }
        for (b = 0; b < n; b += 2 * dual) {
          int i = 2*(b + a);
          int j = 2*(b + a + dual);

          double z1_real = data[j];
          double z1_imag = data[j+1];

          double wd_real = w_real * z1_real - w_imag * z1_imag;
          double wd_imag = w_real * z1_imag + w_imag * z1_real;

          data[j]   = data[i]   - wd_real;
          data[j+1] = data[i+1] - wd_imag;
          data[i]  += wd_real;
          data[i+1]+= wd_imag;
        }
      }
    }
  }
// bitflip

static void FFT_transform_internal_bitflip (int N, double *data, int direction) {
    int n = N/2;
    int bit = 0;
    int logn;
    int dual = 1;

    if (n == 1) return;         /* Identity operation! */
    logn = int_log2(n);


    if (N == 0) return;

    /* bit reverse the input data for decimation in time algorithm */
    FFT_bitreverse_bitflip(N, data) ;

    /* apply fft recursion */
    /* this loop executed int_log2(N) times */
    for (bit = 0; bit < logn; bit++, dual *= 2) {
      double w_real = 1.0;
      double w_imag = 0.0;
      int a;
      int b;

      double theta = 2.0 * direction * PI / (2.0 * (double) dual);
      double s = sin(theta);
      double t = sin(theta / 2.0);
      double s2 = 2.0 * t * t;

      for (a=0, b = 0; b < n; b += 2 * dual) {
        int i = 2*b ;
        int j = 2*(b + dual);

        double wd_real = data[j] ;

#ifdef AGGRESSIVE
        double wd_imag = data[j+1] ;
#else
        double wd_imag = bitflip_float(data[j+1]) ;
#endif
        data[j]   = bitflip_float(data[i]   - wd_real);
        data[j+1] = bitflip_float(data[i+1] - wd_imag);

#ifdef AGGRESSIVE
        data[i]  += bitflip_float(wd_real);
#else
          data[i]  += wd_real;
#endif
        data[i+1]+= (wd_imag);

      }

      /* a = 1 .. (dual-1) */
      for (a = 1; a < dual; a++) {
        /* trignometric recurrence for w-> exp(i theta) w */
        {
          double tmp_real = w_real - s * w_imag - s2 * w_real;
          double tmp_imag = w_imag + s * w_real - s2 * w_imag;
          w_real = tmp_real;
          w_imag = tmp_imag;
        }
        for (b = 0; b < n; b += 2 * dual) {
          int i = 2*(b + a);
          int j = 2*(b + a + dual);

#ifdef AGGRESSIVE
          double z1_real = bitflip_float(data[j]);
#else
          double z1_real = data[j];
#endif
          double z1_imag = data[j+1];

          double wd_real = w_real * z1_real - w_imag * z1_imag;
#ifdef AGGRESSIVE
          double wd_imag = w_real * z1_imag + w_imag * z1_real;
#else
          double wd_imag = bitflip_float(w_real * z1_imag + w_imag * z1_real);
#endif
          data[j]   = bitflip_float(data[i]   - wd_real);
          data[j+1] = bitflip_float(data[i+1] - wd_imag);

#ifdef AGGRESSIVE
          data[i]  += bitflip_float(wd_real);
#else
          data[i]  += (wd_real);
#endif
          data[i+1]+= (wd_imag);
        }
      }
    }
  }

  static void FFT_transform_internal_mpfr (int N, double *data, int direction) {
      int n = N/2;
      int bit = 0;
      int logn;
      int dual = 1;

      if (n == 1) return;         /* Identity operation! */
      logn = int_log2(n);

      mpfr_t w_real_mpfr;
      mpfr_t w_imag_mpfr;
      mpfr_t theta_mpfr;
      mpfr_t s_mpfr;
      mpfr_t s2_mpfr;
      mpfr_t t_mpfr;
      mpfr_init2(w_real_mpfr, FP_APPROX_FRACTION_BIT);
      mpfr_init2(w_imag_mpfr, FP_APPROX_FRACTION_BIT);
      mpfr_init2(theta_mpfr, FP_APPROX_FRACTION_BIT);
      mpfr_init2(s_mpfr, FP_APPROX_FRACTION_BIT);
      mpfr_init2(s2_mpfr, FP_FULL_FRACTION_BIT); //False negative
      mpfr_init2(t_mpfr, FP_APPROX_FRACTION_BIT);

      if (N == 0) return;

      /* bit reverse the input data for decimation in time algorithm */
      FFT_bitreverse_mpfr(N, data) ;

      /* apply fft recursion */
      /* this loop executed int_log2(N) times */
      for (bit = 0; bit < logn; bit++, dual *= 2) {
        double w_real = 1.0;
        double w_imag = 0.0;
        int a;
        int b;

        double theta = 2.0 * direction * PI / (2.0 * (double) dual);

        mpfr_set_d(theta_mpfr, theta, MPFR_RNDZ);
        theta = mpfr_get_d(theta_mpfr,MPFR_RNDZ);

        double s = sin(theta);
        mpfr_set_d(s_mpfr, s, MPFR_RNDZ);
        s = mpfr_get_d(s_mpfr,MPFR_RNDZ);

        double t = sin(theta / 2.0);
        mpfr_set_d(t_mpfr, t, MPFR_RNDZ);
        t = mpfr_get_d(t_mpfr,MPFR_RNDZ);

        double s2 = 2.0 * t * t;
        mpfr_set_d(s2_mpfr, s2, MPFR_RNDZ);
        s2 = mpfr_get_d(s2_mpfr,MPFR_RNDZ);

        for (a=0, b = 0; b < n; b += 2 * dual) {
          int i = 2*b ;
          int j = 2*(b + dual);

          double wd_real = data[j] ;
          double wd_imag = data[j+1] ;

          data[j]   = data[i]   - wd_real;
          data[j+1] = data[i+1] - wd_imag;
          data[i]  += wd_real;
          data[i+1]+= wd_imag;
        }

        /* a = 1 .. (dual-1) */
        for (a = 1; a < dual; a++) {
          /* trignometric recurrence for w-> exp(i theta) w */
          {
            double tmp_real = w_real - s * w_imag - s2 * w_real;
            double tmp_imag = w_imag + s * w_real - s2 * w_imag;
            w_real = tmp_real;
            mpfr_set_d(w_real_mpfr, w_real, MPFR_RNDZ);
            w_real = mpfr_get_d(w_real_mpfr, MPFR_RNDZ);

            w_imag = tmp_imag;
            mpfr_set_d(w_imag_mpfr, w_imag, MPFR_RNDZ);
            w_imag = mpfr_get_d(w_imag_mpfr, MPFR_RNDZ);

          }
          for (b = 0; b < n; b += 2 * dual) {
            int i = 2*(b + a);
            int j = 2*(b + a + dual);

            double z1_real = data[j];
            double z1_imag = data[j+1];

            double wd_real = w_real * z1_real - w_imag * z1_imag;
            double wd_imag = w_real * z1_imag + w_imag * z1_real;

            data[j]   = data[i]   - wd_real;
            data[j+1] = data[i+1] - wd_imag;
            data[i]  += wd_real;
            data[i+1]+= wd_imag;
          }
        }
      }
    }

void FFT_bitreverse(int N, double *data) {
    /* This is the Goldrader bit-reversal algorithm */
    int n=N/2;
    int nm1 = n-1;
    int i=0;
    int j=0;
    for (; i < nm1; i++) {

      /*int ii = 2*i; */
      int ii = i << 1;

      /*int jj = 2*j; */
      int jj = j << 1;

      /* int k = n / 2 ; */
      int k = n >> 1;

      if (i < j) {
        double tmp_real    = data[ii];
        double tmp_imag    = data[ii+1];
        data[ii]   = data[jj];
        data[ii+1] = data[jj+1];
        data[jj]   = tmp_real;
        data[jj+1] = tmp_imag; }

      while (k <= j)
      {
        /*j = j - k ; */
        j -= k;

        /*k = k / 2 ;  */
        k >>= 1 ;
      }
      j += k ;
    }
  }



  void FFT_bitreverse_mpfr(int N, double *data) {
      /* This is the Goldrader bit-reversal algorithm */

      mpfr_t tmp_real_mpfr;
      mpfr_t tmp_imag_mpfr;
      mpfr_t data_mpfr;
      mpfr_init2(tmp_real_mpfr, FP_APPROX_FRACTION_BIT);
      mpfr_init2(tmp_imag_mpfr, FP_APPROX_FRACTION_BIT);
      mpfr_init2(data_mpfr, FP_APPROX_FRACTION_BIT);

      int n=N/2;
      int nm1 = n-1;
      int i=0;
      int j=0;
      for (; i < nm1; i++) {

        /*int ii = 2*i; */
        int ii = i << 1;

        /*int jj = 2*j; */
        int jj = j << 1;

        /* int k = n / 2 ; */
        int k = n >> 1;

        if (i < j) {
          mpfr_set_d(tmp_real_mpfr,  data[ii], MPFR_RNDZ);
          mpfr_set_d(tmp_imag_mpfr,  data[ii+1], MPFR_RNDZ);
          //double tmp_real    = data[ii];
          //double tmp_imag    = data[ii+1];

          mpfr_set_d(data_mpfr,  data[jj], MPFR_RNDZ);
          data[ii] = mpfr_get_d(data_mpfr,MPFR_RNDZ);
        //  data[ii]   = data[jj];

          mpfr_set_d(data_mpfr,  data[jj+1], MPFR_RNDZ);
          data[ii+1] = mpfr_get_d(data_mpfr,MPFR_RNDZ);

        //  data[ii+1] = data[jj+1];

          data[jj] = mpfr_get_d(tmp_real_mpfr,MPFR_RNDZ);
          data[jj+1] = mpfr_get_d(tmp_imag_mpfr,MPFR_RNDZ);
        //data[jj]   = tmp_real;
        //  data[jj+1] = tmp_imag;
       }

        while (k <= j)
        {
          /*j = j - k ; */
          j -= k;

          /*k = k / 2 ;  */
          k >>= 1 ;
        }
        j += k ;
      }
    }



  void FFT_bitreverse_bitflip(int N, double *data) {
      /* This is the Goldrader bit-reversal algorithm */
      int n=N/2;
      int nm1 = n-1;
      int i=0;
      int j=0;
      for (; i < nm1; i++) {

        /*int ii = 2*i; */
        int ii = i << 1;

        /*int jj = 2*j; */
        int jj = j << 1;

        /* int k = n / 2 ; */
        int k = n >> 1;

        if (i < j) {
#ifdef AGGRESSIVE
          double tmp_real    = bitflip_float(data[ii]);
#else
          double tmp_real    = (data[ii]);
#endif
          double tmp_imag    = bitflip_float(data[ii+1]);
#ifdef AGGRESSIVE
          data[ii]   = bitflip_float(data[jj]);
#else
          data[ii]   = (data[jj]);
#endif
          data[ii+1] = (data[jj+1]);
          data[jj]   = (tmp_real);
          data[jj+1] = (tmp_imag); }

        while (k <= j)
        {
          /*j = j - k ; */
          j -= k;

          /*k = k / 2 ;  */
          k >>= 1 ;
        }
        j += k ;
      }
    }


void FFT_transform(int N, double *data)
{
    FFT_transform_internal(N, data, -1);
}


void FFT_inverse(int N, double *data)
{
    int n = N/2;
    double norm = 0.0;
    int i=0;
    FFT_transform_internal(N, data, +1);

    /* Normalize */


    norm=1/((double) n);
    for(i=0; i<N; i++)
      data[i] *= norm;

}

void FFT_transform_mpfr(int N, double *data)
{
    FFT_transform_internal_mpfr(N, data, -1);
}


void FFT_inverse_mpfr(int N, double *data)
{
    int n = N/2;
    double norm = 0.0;

    int i=0;
    FFT_transform_internal_mpfr(N, data, +1);

    /* Normalize */


    norm=1/((double) n);
    mpfr_t norm_mpfr;
    mpfr_init2(norm_mpfr, FP_FULL_FRACTION_BIT);
    mpfr_set_d(norm_mpfr,  1/((double) n), MPFR_RNDZ);
    norm = mpfr_get_d(norm_mpfr,MPFR_RNDZ);

    for(i=0; i<N; i++)
      data[i] *= norm;

}

void FFT_transform_bitflip(int N, double *data)
{
    FFT_transform_internal_bitflip(N, data, -1);
}


void FFT_inverse_bitflip(int N, double *data)
{
    int n = N/2;
    double norm = 0.0;
    int i=0;
    FFT_transform_internal_bitflip(N, data, +1);

    /* Normalize */


    norm=1/((double) n);
    for(i=0; i<N; i++)
      data[i] *= norm;

}





int main()
{
    /* initialize FFT data as complex (N real/img pairs) */
double mintime = RESOLUTION_DEFAULT;


mpfr_t pi_mpfr;

mpfr_init2(pi_mpfr, FP_APPROX_FRACTION_BIT);
mpfr_set_d(pi_mpfr,   3.1415926535897932, MPFR_RNDZ);
  PI = mpfr_get_d(pi_mpfr, MPFR_RNDZ);

 int N = FFT_SIZE;

Random R = new_Random_seed(RANDOM_SEED);


    int twoN = 2*N;
    double *x = RandomVector(twoN, R);
    Random R_bitflip = new_Random_seed(RANDOM_SEED);

    Random R_mpfr = new_Random_seed(RANDOM_SEED);

    double *x_bitflip = RandomVector(twoN, R_bitflip);
    double *x_mpfr = RandomVector(twoN, R_mpfr);
    int i=0;
    mpfr_t data_mpfr;
    mpfr_init2(data_mpfr, FP_APPROX_FRACTION_BIT);
    for(i=0;i<twoN;i++){

  		mpfr_set_d(data_mpfr,x_mpfr[i],MPFR_RNDZ);
      x_mpfr[i] = mpfr_get_d(data_mpfr,MPFR_RNDZ);

  	}

    long cycles = 1;
    Stopwatch Q = new_Stopwatch();

    while(1)
    {
        Stopwatch_start(Q);
        for (i=0; i<cycles; i++)
        {
            FFT_transform(twoN, x);     /* forward transform */
            FFT_inverse(twoN, x);       /* backward transform */

            FFT_transform_bitflip(twoN, x_bitflip);     /* forward transform */
            FFT_inverse_bitflip(twoN, x_bitflip);       /* backward transform */

            FFT_transform_mpfr(twoN, x_mpfr);     /* forward transform */
            FFT_inverse_mpfr(twoN, x_mpfr);       /* backward transform */
        }
        Stopwatch_stop(Q);
        if (Stopwatch_read(Q) >= mintime)
            break;

        cycles *= 2;
    }

    double err =0 ,sum_bitflip = 0, sum_mpfr =0;
    /* approx Mflops */
    for(i=0;i<twoN;i++){
  		err = x[i] - x_bitflip[i];
  		sum_bitflip += (err*err);
      //reuse err
      err = x[i] - x_mpfr[i];
      sum_mpfr += (err*err);
  	}
  	printf("bitflip %lf\n",sqrt(sum_bitflip/twoN));
    printf("mpfr %lf\n",sqrt(sum_mpfr/twoN));



    free(x);

}
