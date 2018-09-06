#ifndef APPROXFFT
#define APPROXFFT
void FFT_transform(int N, double *data);
void FFT_inverse(int N, double *data);
void FFT_bitreverse(int N, double *data);
double FFT_num_flops(int N);

void FFT_transform_bitflip(int N, double *data);
void FFT_inverse_bitflip(int N, double *data);
void FFT_bitreverse_bitflip(int N, double *data);


void FFT_transform_mpfr(int N, double *data);
void FFT_inverse_mpfr(int N, double *data);
void FFT_bitreverse_mpfr(int N, double *data);

#endif
