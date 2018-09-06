/*
 * Stupid simple Raytracer.
 *
 * The C version of the original Java program of EnerJ benchmark.
 * https://sampa.cs.washington.edu/new/research/approximation/enerj.html
 *
 * Portions copyright 2017 National University of Singapore.
 */
#include <stdio.h>
#include <stdlib.h>
#include "bitflip.h"
#include <mpfr.h>
//#include <klee/klee.h>
// mpfr stuff , dont touch
#define LEN 24
 int config_vals[LEN];
  mpfr_t k_init_mpfr;
  mpfr_t lcoff_init_mpfr;
  mpfr_t sng_init_mpfr;
  mpfr_t xe_init_mpfr;
  mpfr_t ye_init_mpfr;
  mpfr_t ze_init_mpfr;
  mpfr_t xd_init_mpfr;
  mpfr_t yd_init_mpfr;
  mpfr_t zd_init_mpfr;
  mpfr_t ix_init_mpfr;
  mpfr_t iy_init_mpfr;
  mpfr_t iz_init_mpfr;
  mpfr_t nx_init_mpfr;
  mpfr_t ny_init_mpfr;
  mpfr_t nz_init_mpfr;
  mpfr_t lx_init_mpfr;
  mpfr_t ly_init_mpfr;
  mpfr_t lz_init_mpfr;
  mpfr_t lly_init_mpfr;
  mpfr_t t_init_mpfr;
  mpfr_t l_init_mpfr;
  mpfr_t w1_init_mpfr;
  mpfr_t h1_init_mpfr;
    mpfr_t temp_var_1;
    mpfr_t temp_var_2;
            mpfr_t temp_var_3;
            mpfr_t temp_var_4;
            mpfr_t temp_var_5;
            mpfr_t temp_var_6;
            mpfr_t temp_var_7;
            mpfr_t temp_var_8;
            mpfr_t temp_var_9;
            mpfr_t temp_var_10;
                mpfr_t temp_var_11;
            mpfr_t temp_var_12;
                mpfr_t temp_var_13;
                mpfr_t temp_var_14;
                mpfr_t temp_var_15;
                mpfr_t temp_var_16;
                mpfr_t temp_var_17;
                mpfr_t temp_var_18;
                mpfr_t temp_var_19;
                mpfr_t temp_var_20;
                mpfr_t temp_var_21;
                mpfr_t temp_var_22;
                mpfr_t temp_var_23;
int my_init_mpfr() {
  mpfr_init2(k_init_mpfr, config_vals[1]);
  //approx
  mpfr_init2(lcoff_init_mpfr, config_vals[2]);
  mpfr_init2(sng_init_mpfr, config_vals[3]);
  mpfr_init2(xe_init_mpfr, config_vals[4]);
  //approx
  mpfr_init2(ye_init_mpfr, config_vals[5]);
  mpfr_init2(ze_init_mpfr, config_vals[6]);
  mpfr_init2(xd_init_mpfr, config_vals[7]);
  mpfr_init2(yd_init_mpfr, config_vals[8]);
  mpfr_init2(zd_init_mpfr, config_vals[9]);
  mpfr_init2(ix_init_mpfr, config_vals[10]);
  //approx
  mpfr_init2(iy_init_mpfr, config_vals[11]);
  //approx
  mpfr_init2(iz_init_mpfr, config_vals[12]);
  mpfr_init2(nx_init_mpfr, config_vals[13]);
  mpfr_init2(ny_init_mpfr, config_vals[14]);
  mpfr_init2(nz_init_mpfr, config_vals[15]);
  mpfr_init2(lx_init_mpfr, config_vals[16]);
  //approx
  mpfr_init2(ly_init_mpfr, config_vals[17]);
  //approx
  mpfr_init2(lz_init_mpfr, config_vals[18]);
  //approx
  mpfr_init2(lly_init_mpfr, config_vals[19]);
  //approx
  mpfr_init2(t_init_mpfr, config_vals[20]);
  mpfr_init2(l_init_mpfr, config_vals[21]);
  mpfr_init2(w1_init_mpfr, config_vals[22]);
  mpfr_init2(h1_init_mpfr, config_vals[23]);

    mpfr_init2 (temp_var_1, config_vals[0]);
    mpfr_init2 (temp_var_2, config_vals[0]);
            mpfr_init2 (temp_var_3, config_vals[0]);
            mpfr_init2 (temp_var_4, config_vals[0]);
            mpfr_init2 (temp_var_5, config_vals[0]);
            mpfr_init2 (temp_var_6, config_vals[0]);
            mpfr_init2 (temp_var_7, config_vals[0]);
            mpfr_init2 (temp_var_8, config_vals[0]);
            mpfr_init2 (temp_var_9, config_vals[0]);
            mpfr_init2 (temp_var_10, config_vals[0]);
                mpfr_init2 (temp_var_11, config_vals[0]);
            mpfr_init2 (temp_var_12, config_vals[0]);
                mpfr_init2 (temp_var_13, config_vals[0]);
                mpfr_init2 (temp_var_14, config_vals[0]);
                mpfr_init2 (temp_var_15, config_vals[0]);
                mpfr_init2 (temp_var_16, config_vals[0]);
                mpfr_init2 (temp_var_17, config_vals[0]);
                mpfr_init2 (temp_var_18, config_vals[0]);
                mpfr_init2 (temp_var_19, config_vals[0]);
                mpfr_init2 (temp_var_20, config_vals[0]);
                mpfr_init2 (temp_var_21, config_vals[0]);
                mpfr_init2 (temp_var_22, config_vals[0]);
                mpfr_init2 (temp_var_23, config_vals[0]);
}

int init_readconfig() {

  // For reading precision contents of config_file into array
   FILE *myFile;
     myFile = fopen("plane_config_file.txt", "r");

        if (myFile == NULL) {
					printf("Error Reading File\n");
                exit (0);
                }

        int s;
        for (s = 0; s < LEN; s++) {
            fscanf(myFile, "%d,", &config_vals[s]);
                          }
        config_vals[0]= 53;
        fclose(myFile);
        my_init_mpfr();
        return 0;
}



// end mpfr stuff


 int  w = 400;
 int  h = 256;

int abs(int value) {
  if (value < 0)
    return -value;
  return value;
}

// From
// https://stackoverflow.com/questions/3581528/how-is-the-square-root-function-implemented
double sqrt(double n) {
  double lo = 0, hi = n, mid;
  int i;
  for ( i = 0; i < 1000; i++) {
    mid = (lo + hi) / 2;
    if (mid * mid == n)
      return mid;
    if (mid * mid > n) {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  return mid;
}

// From
// https://stackoverflow.com/questions/4572556/concise-way-to-implement-round-in-c
int my_round(double x) {
  if (x < 0.0)
    return (int)(x - 0.5);
  else
    return (int)(x + 0.5);
}

/* Approx */
int texture1(
    /* Approx */
    float x,
    /* Approx */
    float y,
    /* Approx */
    float z,
    float lcoff,
    int texture,
    int light
    ) {
  int v;
  /* Approx */
  int col;
  /* Approx */
  int r, g, b;
  r = 255;
  b = 0;
  col = 0;
  if (light != 0) {
    r = (
        /* Approx */
        int)(255 * lcoff);
  }
  b = r;

  if (texture == 1) {
    col = (255 << 24) | (255 << 16);
  } else if (texture == 2) {
    v = (my_round(x) + my_round(z)) % 2;
    if (v == 0) {
      col = (255 << 24) | b;
    } else {
      col = (255 << 24) | (r << 16);
    }
  }

  return col;
}

/* Approx */
int texture1_bitflip(
    /* Approx */
    float x,
    /* Approx */
    float y,
    /* Approx */
    float z,
    float lcoff,
    int texture,
    int light
    ) {
  int v;
  /* Approx */
  int col;
  /* Approx */
  int r, g, b;
  r = 255;
  b = 0;
  col = 0;
  if (bitflip_int(light) != 0) {
    r = (
        /* Approx */
        int)(255 * lcoff);
  }

  b = r;
  // inject error to r abd b
#ifdef AGGRESSIVE
  b= bitflip_int(b);
  r= bitflip_int(r);
#endif
  if (bitflip_int(texture) == 1) {
    col = (255 << 24) | (255 << 16);
  } else if (bitflip_int(texture) == 2) {

#ifdef AGRESSIVE
    v = (my_round(x) + my_round(bitflip_float(z))) % 2;
#else
    v = (my_round(x) + my_round(z)) % 2;
#endif

    if (v == 0) {
      col = (255 << 24) | b;
    } else {
      col = (255 << 24) | (r << 16);
    }
  }

  return col;
}

void init(int arg0, int arg1, int arg2, int arg3, int* pixels) {
  w = 400;
  h = 256;
  //~ w = 40;
  //~ h = 25;
 float k; // what the hell is this variable for?
/* Approx */
int texture, light;
float lcoff;
/* Approx */

float sng; // could maybe make approximate
int numIterations = 0;


  texture = arg0; // getParameter("texture"));
  light = arg1;   // getParameter("light"));
  /* Approx */
  //~ int pixels[w * h];

  int index, x, y; // not approx --> for loops and array indexing.
  /* Approx */
  float xe, ye, ze, xd, yd, zd;
  /* Approx */
  float ix, iy, iz;
  /* Approx */
  float nx, ny, nz;
  /* Approx */
  float lx, ly, lz;
  float lly;
  lly = arg2; // getParameter("lighty"));
  ye = arg3;  // getParameter("viewy"));

  nx = 0;
  ny = 1;
  nz = 0;
  int bl = (255 << 24); // this stands for black, a constant, maybe?
  /* Approx */
  float t; // who knows
  /* Approx */
  float l;      // who knows
  float w1, h1; // positioning in image? so don't make approx?
  w1 = w / 2;
  h1 = h / 2;

  xe = 0;

  ze = 0;
  k = -1;

  for (y = 0; y < h; y++) {
    for (x = 0; x < w; x++) {
      t = -1;
      xd = (x - w1) / w1;
      yd = (h1 - y) / h1;
      zd = -1;
      l = xd * xd + yd * yd + zd * zd;
      xd /= l;
      yd /= l;
      zd /= l;

      if ((k - ye) * yd <= 0) {
        t = -1;
      } else {
        t = (k - ye) / yd;
      }

      index = y * w + x;
      if (t >= 0) {
        ix = xe + t * xd;
        iy = ye + t * yd;
        iz = ze + t * zd;
        lx = 0;
        ly = lly;
        lz = 0;
        lx = lx - ix;
        ly = ly - iy;
        lz = lz - iz;
        sng = (float)sqrt(lx * lx + ly * ly + lz * lz); // sng=1.7f/sng;
        sng = 1.0f / sng;
        lcoff = (lx * nx + ly * ny + lz * nz) * sng;
        pixels[index] = texture1(ix, iy, iz, lcoff, texture,light);
      } else {
        pixels[index] = bl;
      }
      numIterations++;
    }
  }

}


void init_mpfr(int arg0, int arg1, int arg2, int arg3, int *pixels)
{
  w = 400;
  h = 256;
  int texture;
  int light;
  int numIterations = 0;
  texture = (arg0);
  light = (arg1);
  int index;
  int x;
  int y;
  mpfr_set_d(lly_init_mpfr, (arg2), MPFR_RNDN);
  mpfr_set_d(ye_init_mpfr, (arg3), MPFR_RNDN);
  mpfr_set_d(nx_init_mpfr, 0, MPFR_RNDN);
  mpfr_set_d(ny_init_mpfr, 1, MPFR_RNDN);
  mpfr_set_d(nz_init_mpfr, 0, MPFR_RNDN);
  int bl = 255 << 24;
  mpfr_set_d(w1_init_mpfr, w / 2, MPFR_RNDN);
  mpfr_set_d(h1_init_mpfr, h / 2, MPFR_RNDN);
  mpfr_set_d(xe_init_mpfr, 0, MPFR_RNDN);
  mpfr_set_d(ze_init_mpfr, 0, MPFR_RNDN);
  mpfr_set_d(k_init_mpfr, -1, MPFR_RNDN);
  for (y = 0; y < h; y++){
  {
    for (x = 0; x < w; x++){
    {
      mpfr_set_d(t_init_mpfr, -1, MPFR_RNDN);



            mpfr_d_sub(temp_var_3, x, w1_init_mpfr, MPFR_RNDN);
            mpfr_div(xd_init_mpfr, (temp_var_3), w1_init_mpfr, MPFR_RNDN);

            mpfr_sub_d(temp_var_4, h1_init_mpfr, y, MPFR_RNDN);
            mpfr_div(yd_init_mpfr, (temp_var_4), h1_init_mpfr, MPFR_RNDN);
      mpfr_set_d(zd_init_mpfr, -1, MPFR_RNDN);

            mpfr_mul(temp_var_5, xd_init_mpfr, xd_init_mpfr, MPFR_RNDN);

            mpfr_mul(temp_var_6, yd_init_mpfr, yd_init_mpfr, MPFR_RNDN);

            mpfr_add(temp_var_7, (temp_var_5), (temp_var_6), MPFR_RNDN);

            mpfr_mul(temp_var_8, zd_init_mpfr, zd_init_mpfr, MPFR_RNDN);
            mpfr_add(l_init_mpfr, (temp_var_7), (temp_var_8), MPFR_RNDN);
            mpfr_div(xd_init_mpfr, xd_init_mpfr, l_init_mpfr, MPFR_RNDN);
            mpfr_div(yd_init_mpfr, yd_init_mpfr, l_init_mpfr, MPFR_RNDN);
            mpfr_div(zd_init_mpfr, zd_init_mpfr, l_init_mpfr, MPFR_RNDN);

            mpfr_sub(temp_var_9, k_init_mpfr, ye_init_mpfr, MPFR_RNDN);

            mpfr_mul(temp_var_10, (temp_var_9), yd_init_mpfr, MPFR_RNDN);
if ((mpfr_cmp_d((temp_var_10),0) <= 0)
//original  (temp_var_10) <= 0
)
      {
        mpfr_set_d(t_init_mpfr, -1, MPFR_RNDN);
      }
      else
      {

                mpfr_sub(temp_var_11, k_init_mpfr, ye_init_mpfr, MPFR_RNDN);
                mpfr_div(t_init_mpfr, (temp_var_11), yd_init_mpfr, MPFR_RNDN);
      }

      index = (y * w) + x;

if ((mpfr_cmp_d(t_init_mpfr,0) >= 0)
//original  t_init_mpfr >= 0
)
      {

                mpfr_mul(temp_var_13, t_init_mpfr, xd_init_mpfr, MPFR_RNDN);
                mpfr_add(ix_init_mpfr, xe_init_mpfr, (temp_var_13), MPFR_RNDN);

                mpfr_mul(temp_var_14, t_init_mpfr, yd_init_mpfr, MPFR_RNDN);
                mpfr_add(iy_init_mpfr, ye_init_mpfr, (temp_var_14), MPFR_RNDN);

                mpfr_mul(temp_var_15, t_init_mpfr, zd_init_mpfr, MPFR_RNDN);
                mpfr_add(iz_init_mpfr, ze_init_mpfr, (temp_var_15), MPFR_RNDN);
        mpfr_set_d(lx_init_mpfr, 0, MPFR_RNDN);
        mpfr_set(ly_init_mpfr, lly_init_mpfr, MPFR_RNDN);
        mpfr_set_d(lz_init_mpfr, 0, MPFR_RNDN);
                        mpfr_sub(lx_init_mpfr, lx_init_mpfr, ix_init_mpfr, MPFR_RNDN);
                        mpfr_sub(ly_init_mpfr, ly_init_mpfr, iy_init_mpfr, MPFR_RNDN);
                        mpfr_sub(lz_init_mpfr, lz_init_mpfr, iz_init_mpfr, MPFR_RNDN);
        mpfr_set_d(sng_init_mpfr, (float) sqrt(((mpfr_get_d(lx_init_mpfr, MPFR_RNDN) * mpfr_get_d(lx_init_mpfr, MPFR_RNDN)) + (mpfr_get_d(ly_init_mpfr, MPFR_RNDN) * mpfr_get_d(ly_init_mpfr, MPFR_RNDN))) + (mpfr_get_d(lz_init_mpfr, MPFR_RNDN) * mpfr_get_d(lz_init_mpfr, MPFR_RNDN))), MPFR_RNDN);



                mpfr_d_div(sng_init_mpfr, 1.0f, sng_init_mpfr, MPFR_RNDN);

                mpfr_mul(temp_var_19, lx_init_mpfr, nx_init_mpfr, MPFR_RNDN);

                mpfr_mul(temp_var_20, ly_init_mpfr, ny_init_mpfr, MPFR_RNDN);

                mpfr_add(temp_var_21, (temp_var_19), (temp_var_20), MPFR_RNDN);

                mpfr_mul(temp_var_22, lz_init_mpfr, nz_init_mpfr, MPFR_RNDN);

                mpfr_add(temp_var_23, (temp_var_21), (temp_var_22), MPFR_RNDN);
                mpfr_mul(lcoff_init_mpfr, (temp_var_23), sng_init_mpfr, MPFR_RNDN);
        pixels[index] = texture1(mpfr_get_d(ix_init_mpfr, MPFR_RNDN), mpfr_get_d(iy_init_mpfr, MPFR_RNDN), mpfr_get_d(iz_init_mpfr, MPFR_RNDN), mpfr_get_d(lcoff_init_mpfr, MPFR_RNDN), bitflip_int(texture), bitflip_int(light));
      }
      else
      {
        pixels[index] = bl;
      }

      numIterations++;
    }

        }

  }

    }

}


void init_bitflip(int arg0, int arg1, int arg2, int arg3, int* pixels) {
  w = 400;
  h = 256;
  //~ w = 40;
  //~ h = 25;
	float k; // what the hell is this variable for?
	/* Approx */
	int texture, light;
	/* Approx */
	float lcoff;
	float sng; // could maybe make approximate
	int numIterations = 0;


  texture = arg0; // getParameter("texture"));
  light = bitflip_int(arg1);   // getParameter("light"));
  /* Approx */
  //~ int pixels[w * h];

  int index, x, y; // not approx --> for loops and array indexing.
  /* Approx */
  float xe, ye, ze, xd, yd, zd;
  /* Approx */
  float ix, iy, iz;
  /* Approx */
  float nx, ny, nz;
  /* Approx */
  float lx, ly, lz;
  float lly;
  lly = bitflip_int(arg2); // getParameter("lighty"));
#ifdef AGRESSIVE
    ye = bitflip_int(arg3);  // getParameter("viewy"));
#else
    ye = (arg3);
#endif

  nx = 0;
  ny = 1;
  nz = 0;
  int bl = (255 << 24); // this stands for black, a constant, maybe?
  /* Approx */
  float t; // who knows
  /* Approx */
  float l;      // who knows
  float w1, h1; // positioning in image? so don't make approx?
  w1 = w / 2;
  h1 = h / 2;

  xe = 0;

  ze = 0;
  k = -1;

  for (y = 0; y < h; y++) {
    for (x = 0; x < w; x++) {
      t = -1;
      xd = (x - w1) / w1;
      yd = (h1 - y) / h1;
      zd = -1;
      l = xd * xd + yd * yd + zd * zd;
      xd /= l;
      yd /= l;
      zd /= l;

      if ((k - ye) * yd <= 0) {
        t = -1;
      } else {
        t = (k - ye) / yd;
      }
#ifdef AGGRESSIVE
		t = bitflip_float(t);
#endif
      index = y * w + x;
      if (t >= 0) {
        ix = xe + t * xd;
#ifdef AGGRESSIVE
        iy = bitflip_float(ye + t * yd);
        iz = bitflip_float(ze + t * zd);
#else
        iy = (ye + t * yd);
        iz = (ze + t * zd);
#endif
        lx = 0;
        ly = lly;
        lz = 0;
        lx = lx - ix;
#ifdef AGGRESSIVE
        ly = bitflip_float(ly - iy);
#else
        ly = (ly - iy);
#endif
        lz = lz - iz;

#ifdef AGGRESSIVE
        sng = (float)sqrt(bitflip_float(lx * lx + ly * ly + lz * lz)); // sng=1.7f/sng;
#else
        sng = (float)sqrt(lx * lx + ly * ly + lz * lz); // sng=1.7f/sng;
#endif
        sng = 1.0f / sng;
        lcoff = (lx * nx + ly * ny + lz * nz) * sng;
        pixels[index] = texture1_bitflip(ix, iy, iz, lcoff, texture,light);
      } else {
        pixels[index] = bl;
      }
      numIterations++;
    }
  }

}

int main(int argc, char **argv) {
  int arg0 =2  , arg1 = 1, arg2 = 30, arg3= 10;
//  int arg0 =1  , arg1 = 16843009, arg2 = 30, arg3= 10;

	init_readconfig(); //for mpfr

	int * pixels_ref = malloc(w*h*sizeof(int));
	int * pixels_bitflip = malloc(w*h*sizeof(int));
	int * pixels_mpfr = malloc(w*h*sizeof(int));

  init(arg0, arg1, arg2, arg3, pixels_ref);
//  init_bitflip(arg0, arg1, arg2, arg3, pixels_bitflip);
  init_mpfr(arg0, arg1, arg2, arg3, pixels_mpfr);

  int i;
  double diff_bitflip = 0.0;
  double diff_mpfr = 0.0;
  for ( i = 0; i < w * h; i++) {
    int ref = pixels_ref[i] & 0xff;
     int approx_bitflip = pixels_bitflip[i] & 0xff;
    int approx_mpfr = pixels_mpfr[i] & 0xff;
    diff_bitflip+=abs(approx_bitflip - ref);
    diff_mpfr+=abs(approx_mpfr - ref);

  }

  for ( i = 0; i < w * h; i++) {
    int ref = (pixels_ref[i]>>8) & 0xff;
     int approx_bitflip = (pixels_bitflip[i] >>8)& 0xff;
    int approx_mpfr = (pixels_mpfr[i] >>8)& 0xff;
    diff_bitflip+=abs(approx_bitflip - ref);
    diff_mpfr+=abs(approx_mpfr - ref);

  }

 for ( i = 0; i < w * h; i++) {
    int ref = (pixels_ref[i]>>16) & 0xff;
     int approx_bitflip = (pixels_bitflip[i] >>16)& 0xff;
    int approx_mpfr = (pixels_mpfr[i] >>16)& 0xff;
    diff_bitflip+=abs(approx_bitflip - ref);
    diff_mpfr+=abs(approx_mpfr - ref);

  }


  printf("%lf ", diff_mpfr/(3*w*h));
///  printf("mean diff  bitflip %lf \n", diff_bitflip/(3*w*h));
}
