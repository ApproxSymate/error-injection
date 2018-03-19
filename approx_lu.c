#include <math.h>
//#include "LU.h"
#include "Random.h"
#include "Stopwatch.h"
#include "SparseCompRow.h"
#include "array.h"
#include "constants.h"
#include <mpfr.h>
#include <stdio.h>
#include <stdlib.h>
#include "bitflip.h"


double LU_num_flops(int N)
{
        /* rougly 2/3*N^3 */

    double Nd = (double) N;

    return (2.0 * Nd *Nd *Nd/ 3.0);
}


double LU_diff_matrix(int M, int N, double **lu, double **A)
{
            int i;
                int j;
                  double sum=0.0;
                    for (i=0; i<M; i++){
                            for (j=0; j<N; j++){
                                    double err = lu[i][j]-A[i][j];
                                    sum+=err*err;
                            }
                    }
                    return sqrt(sum/(M*N));
}

void LU_copy_matrix(int M, int N, double **lu, double **A)
{
    int i;
    int j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            lu[i][j] = A[i][j];
}


int LU_factor(int M, int N, double **A,  int *pivot)
{


    int minMN =  M < N ? M : N;
    int j=0;

    for (j=0; j<minMN; j++)
    {
        /* find pivot in column j and  test for singularity. */

        int jp=j;
        int i;

        double t = fabs(A[j][j]);
        for (i=j+1; i<M; i++)
        {
            double ab = fabs(A[i][j]);
            if ( ab > t)
            {
                jp = i;
                t = ab;
            }
        }

        pivot[j] = jp;

        /* jp now has the index of maximum element  */
        /* of column j, below the diagonal          */

        if ( A[jp][j] == 0 )
            return 1;       /* factorization failed because of zero pivot */


        if (jp != j)
        {
            /* swap rows j and jp */
            double *tA = A[j];
            A[j] = A[jp];
            A[jp] = tA;
        }

        if (j<M-1)                /* compute elements j+1:M of jth column  */
        {
            /* note A(j,j), was A(jp,p) previously which was */
            /* guarranteed not to be zero (Label #1)         */

            double recp =  1.0 / A[j][j];
            int k;
            for (k=j+1; k<M; k++)
                A[k][j] *= recp;
        }


        if (j < minMN-1)
        {
            /* rank-1 update to trailing submatrix:   E = E - x*y; */
            /* E is the region A(j+1:M, j+1:N) */
            /* x is the column vector A(j+1:M,j) */
            /* y is row vector A(j,j+1:N)        */

            int ii;
            for (ii=j+1; ii<M; ii++)
            {
                double *Aii = A[ii];
                double *Aj = A[j];
                double AiiJ = Aii[j];
                int jj;
                for (jj=j+1; jj<N; jj++)
                  Aii[jj] -= AiiJ * Aj[jj];

            }
        }
    }

    return 0;
}

int LU_factor_bitflip(int M, int N, double **A,  int *pivot)
{


    int minMN =  M < N ? M : N;
    int j=0;

    for (j=0; j<minMN; j++)
    {
        /* find pivot in column j and  test for singularity. */

        int jp=j;
        int i;

        double t = fabs(bitflip_float(A[j][j]));
        for (i=j+1; i<M; i++)
        {
            double ab = fabs(bitflip_float(A[i][j]));
            if ( ab > t)
            {
                jp = i;
                t = ab;
            }
        }

        pivot[j] = bitflip_int(jp);

        /* jp now has the index of maximum element  */
        /* of column j, below the diagonal          */
        A[jp][j] = bitflip_float(A[jp][j]);
        if ( A[jp][j] == 0 )
            return 1;       /* factorization failed because of zero pivot */


        if (jp != j)
        {
            /* swap rows j and jp */
            double *tA = A[j];
            A[j] = A[jp];
            A[jp] = tA;
        }

        if (j<M-1)                /* compute elements j+1:M of jth column  */
        {
            /* note A(j,j), was A(jp,p) previously which was */
            /* guarranteed not to be zero (Label #1)         */

            double recp =  bitflip_float(1.0 / A[j][j]);
            int k;
            for (k=j+1; k<M; k++)
                A[k][j] *= bitflip_float(recp);
        }


        if (j < minMN-1)
        {
            /* rank-1 update to trailing submatrix:   E = E - x*y; */
            /* E is the region A(j+1:M, j+1:N) */
            /* x is the column vector A(j+1:M,j) */
            /* y is row vector A(j,j+1:N)        */

            int ii;
            for (ii=j+1; ii<M; ii++)
            {
                double *Aii = A[ii];
                double *Aj = A[j];
                double AiiJ = bitflip_float(Aii[j]);
                int jj;
                for (jj=j+1; jj<N; jj++)
                  Aii[jj] -= bitflip_float(AiiJ * Aj[jj]);

            }
        }
    }

    return 0;
}

int LU_factor_mpfr(int M, int N, double **A,  int *pivot)
{

    mpfr_t t_mpfr;
    mpfr_t ab_mpfr;
    mpfr_t A_mpfr; //for all vars in double ** A
    mpfr_t recp_mpfr;
    mpfr_t  AiiJ_mpfr;
    mpfr_init2(t_mpfr, FP_APPROX_FRACTION_BIT);
    mpfr_init2(ab_mpfr, FP_APPROX_FRACTION_BIT);
    mpfr_init2(A_mpfr, FP_APPROX_FRACTION_BIT);
    mpfr_init2(recp_mpfr, FP_APPROX_FRACTION_BIT);
    mpfr_init2(AiiJ_mpfr, FP_APPROX_FRACTION_BIT);

    int minMN =  M < N ? M : N;
    int j=0;

    for (j=0; j<minMN; j++)
    {
        /* find pivot in column j and  test for singularity. */

        int jp=j;
        int i;

      //  double t = fabs(bitflip_float(A[j][j]));
        double t;
        mpfr_set_d(t_mpfr,fabs(A[j][j]),MPFR_RNDZ);
        t = mpfr_get_d(t_mpfr,MPFR_RNDZ);

        for (i=j+1; i<M; i++)
        {
            //double ab = fabs(bitflip_float(A[i][j]));
            double ab;
            mpfr_set_d(ab_mpfr,fabs(A[i][j]),MPFR_RNDZ);
            ab = mpfr_get_d(ab_mpfr,MPFR_RNDZ);

            if ( ab > t)
            {
                jp = i;
                t = ab;
            }
        }

        pivot[j] = bitflip_int(jp);

        /* jp now has the index of maximum element  */
        /* of column j, below the diagonal          */
        mpfr_set_d(A_mpfr,A[jp][j],MPFR_RNDZ);
        A[jp][j] = mpfr_get_d(A_mpfr,MPFR_RNDZ);

        if ( A[jp][j] == 0 )
            return 1;       /* factorization failed because of zero pivot */


        if (jp != j)
        {
            /* swap rows j and jp */
            double *tA = A[j];
            A[j] = A[jp];
            A[jp] = tA;
        }

        if (j<M-1)                /* compute elements j+1:M of jth column  */
        {
            /* note A(j,j), was A(jp,p) previously which was */
            /* guarranteed not to be zero (Label #1)         */

            //double recp =  bitflip_float(1.0 / A[j][j]);
            double recp;
            mpfr_set_d(recp_mpfr,1.0 / A[j][j],MPFR_RNDZ);
            recp = mpfr_get_d(recp_mpfr,MPFR_RNDZ);

            int k;
            for (k=j+1; k<M; k++)
                A[k][j] *= recp;
        }


        if (j < minMN-1)
        {
            /* rank-1 update to trailing submatrix:   E = E - x*y; */
            /* E is the region A(j+1:M, j+1:N) */
            /* x is the column vector A(j+1:M,j) */
            /* y is row vector A(j,j+1:N)        */

            int ii;
            for (ii=j+1; ii<M; ii++)
            {
                double *Aii = A[ii];
                double *Aj = A[j];
                double AiiJ ;//= bitflip_float(Aii[j]);
                mpfr_set_d(AiiJ_mpfr,Aii[j],MPFR_RNDZ);
                AiiJ = mpfr_get_d(AiiJ_mpfr,MPFR_RNDZ);

                int jj;
                for (jj=j+1; jj<N; jj++){
                  mpfr_set_d(A_mpfr,AiiJ * Aj[jj],MPFR_RNDZ);
                  Aii[jj] -= mpfr_get_d(A_mpfr,MPFR_RNDZ);
                  //Aii[jj] -= AiiJ * Aj[jj];
                }
            }
        }
    }

    return 0;
}


int main(){
    FILE* args;
    char line[256];
    double **A = NULL;
    double **lu = NULL;
    double **lu_bitflip = NULL;
    double **lu_mpfr = NULL;
    int *pivot = NULL;
    int N = LU_SIZE;
    double min_time = RESOLUTION_DEFAULT;
    Random R = new_Random_seed(RANDOM_SEED);
    Stopwatch Q = new_Stopwatch();
    int i=0;
    int cycles=1;

    if ((A = RandomMatrix(N, N,  R)) == NULL) exit(1);
    if ((lu = new_Array2D_double(N, N)) == NULL) exit(1);
    if ((lu_bitflip = new_Array2D_double(N, N)) == NULL) exit(1);
    if ((lu_mpfr = new_Array2D_double(N, N)) == NULL) exit(1);
    if ((pivot = (int *) malloc(N * sizeof(int))) == NULL) exit(1);


    while(1)
    {
        Stopwatch_start(Q);
        for (i=0; i<cycles; i++)
        {
            Array2D_double_copy(N, N, lu, A);
            Array2D_double_copy(N, N, lu_bitflip, A);
            Array2D_double_copy(N, N, lu_mpfr, A);
            LU_factor(N, N, lu, pivot);
            LU_factor_bitflip(N, N, lu_bitflip, pivot);
            LU_factor_mpfr(N, N, lu_mpfr, pivot);
        }
        Stopwatch_stop(Q);
        if (Stopwatch_read(Q) >= min_time) break;

        cycles *= 2;
    }

    printf("bitflip %lf\n",LU_diff_matrix(N,N,lu,lu_bitflip));
    printf("mpfr %lf\n",LU_diff_matrix(N,N,lu,lu_mpfr));


    Stopwatch_delete(Q);
    free(pivot);
    Array2D_double_delete(N, N, lu);
    Array2D_double_delete(N, N, A);

  return 0;
}
