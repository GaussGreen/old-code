#ifndef NUM_H_FFT_H
#define NUM_H_FFT_H

void FFT1D(double *data, unsigned long nn, int isign);

void FFTboundary_out(long int n, double *yin, double *yy, int iconv);

void FFTboundary_in(long int n, double *y, double *ycompl);

void Conv1D(double *in_v, unsigned long n, unsigned long n_conv, double delta,
            int iconv);

void ConvolveMarginals(double **AbscMarg, double **OrdMarg, unsigned long n_pts,
                       unsigned long n_assets, unsigned long n_conv,
                       int conv_flag // this is 1 if you want to convolve
                                     // and -1 if you want to deconvolve
);

void twofft(double data1[], double data2[], double fft1[], double fft2[],
            unsigned long n);

#endif
