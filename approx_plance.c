/*
 * Stupid simple Raytracer.
 *
 * The C version of the original Java program of EnerJ benchmark.
 * https://sampa.cs.washington.edu/new/research/approximation/enerj.html
 *
 * Portions copyright 2017 National University of Singapore.
 */
#include <stdio.h>

//#include <klee/klee.h>

 int  w = 400;
 int  h = 256;
int texture, light;
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
    float lcoff
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
        pixels[index] = texture1(ix, iy, iz, lcoff);
      } else {
        pixels[index] = bl;
      }
      numIterations++;
    }
  }
  int i;
  //~ for ( i = 0; i < w * h; i++) {
    //~ printf("%d, ", pixels[i] & 0xff);
  //~ }
  //~ printf("\n");

  //~ for ( i = 0; i < w * h; i+s+) {
    //~ printf("%d, ", (pixels[i] >> 8) & 0xff);
  //~ }
  //~ printf("\n");
  
   //~ for ( i = 0; i < w * h; i++) {
    //~ printf("%d, ", (pixels[i] >> 16) & 0xff);
  //~ }
  //~ printf("\n");
  
  
   //~ checkErrors(pixels);
}

void init_mpfr(int arg0, int arg1, int arg2, int arg3, int* pixels) {
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
        pixels[index] = texture1(ix, iy, iz, lcoff);
      } else {
        pixels[index] = bl;
      }
      numIterations++;
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
        pixels[index] = texture1(ix, iy, iz, lcoff);
      } else {
        pixels[index] = bl;
      }
      numIterations++;
    }
  }

}

int main(int argc, char **argv) {
  int arg0 =2  , arg1 = 1, arg2 = 30, arg3= 10;
  texture = arg0; // getParameter("texture"));
  light = arg1;   // getParameter("light"));
	int * pixels_ref = malloc(w*h*sizeof(int));
	int * pixels_bitflip = malloc(w*h*sizeof(int));
	int * pixels_mpfr = malloc(w*h*sizeof(int));
	
  init(arg0, arg1, arg2, arg3, pixels_ref);
  init_bitflip(arg0, arg1, arg2, arg3, pixels_bitflip);
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
    int ref = pixels_ref[i] & 0xff;
     int approx_bitflip = (pixels_bitflip[i] >>8)& 0xff;
    int approx_mpfr = (pixels_mpfr[i] >>8)& 0xff;
    diff_bitflip+=abs(approx_bitflip - ref);
    diff_mpfr+=abs(approx_mpfr - ref);
    
  }

 for ( i = 0; i < w * h; i++) {
    int ref = pixels_ref[i] & 0xff;
     int approx_bitflip = (pixels_bitflip[i] >>16)& 0xff;
    int approx_mpfr = (pixels_mpfr[i] >>16)& 0xff;
    diff_bitflip+=abs(approx_bitflip - ref);
    diff_mpfr+=abs(approx_mpfr - ref);
    
  } 
  
  
  printf("mean diff  mpfr %lf \n", diff_mpfr/(3*w*h));
  printf("mean diff  bitflip %lf \n", diff_bitflip/(3*w*h));
}
