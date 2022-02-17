#include "UTALLHDR.H>
#include "math.h"
#include "num_h_allhdr.h"

// S. Galluccio: 21 November 2000
// Created Source File

#define FFTSWAP(a, b)                                                          \
  tempr = (a);                                                                 \
  (a) = (b);                                                                   \
  (b) = tempr

void twofft(double data1[], double data2[], double fft1[], double fft2[],
            unsigned long n) {
  void four1(double data[], unsigned long nn, int isign);
  unsigned long nn3, nn2, jj, j;
  double rep, rem, aip, aim;

  nn3 = 1 + (nn2 = 2 + n + n);
  for (j = 1, jj = 2; j <= n; j++, jj += 2) {
    fft1[jj - 1] = data1[j];
    fft1[jj] = data2[j];
  }
  four1(fft1, n, 1);
  fft2[1] = fft1[2];
  fft1[2] = fft2[2] = 0.0;
  for (j = 3; j <= n + 1; j += 2) {
    rep = 0.5 * (fft1[j] + fft1[nn2 - j]);
    rem = 0.5 * (fft1[j] - fft1[nn2 - j]);
    aip = 0.5 * (fft1[j + 1] + fft1[nn3 - j]);
    aim = 0.5 * (fft1[j + 1] - fft1[nn3 - j]);
    fft1[j] = rep;
    fft1[j + 1] = aim;
    fft1[nn2 - j] = rep;
    fft1[nn3 - j] = -aim;
    fft2[j] = aip;
    fft2[j + 1] = -rem;
    fft2[nn2 - j] = aip;
    fft2[nn3 - j] = rem;
  }
}

/*****************************************************************************************/

void FFT1D(double *data, unsigned long nn, int isign) {
  unsigned long n, mmax, m, j, istep, i;
  double wtemp, wr, wpr, wpi, wi, theta;
  double tempr, tempi;

  n = nn << 1;
  j = 1;
  for (i = 1; i < n; i += 2) {
    if (j > i) {
      FFTSWAP(data[j], data[i]);
      FFTSWAP(data[j + 1], data[i + 1]);
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax = 2;
  while (n > mmax) {
    istep = mmax << 1;
    theta = isign * (6.28318530717959 / mmax);
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m = 1; m < mmax; m += 2) {
      for (i = m; i <= n; i += istep) {
        j = i + mmax;
        tempr = wr * data[j] - wi * data[j + 1];
        tempi = wr * data[j + 1] + wi * data[j];
        data[j] = data[i] - tempr;
        data[j + 1] = data[i + 1] - tempi;
        data[i] += tempr;
        data[i + 1] += tempi;
      }
      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
    }
    mmax = istep;
  }
}

/*************************************************************************************/

void FFTboundary_out(long int n, double *yin, double *yy, int iconv) {
  long int i, idx;

  switch (iconv) {

  case 1:

    idx = 1;
    for (i = 1; i <= n; i++) {
      yy[i] = sqrt(yin[idx] * yin[idx] + yin[idx + 1] * yin[idx + 1]);
      idx += 2;
    }

    break;

  case -1:

    idx = n + 1;
    for (i = 1; i <= n / 2; i++) {
      yy[i] = sqrt(yin[idx] * yin[idx] + yin[idx + 1] * yin[idx + 1]);
      idx += 2;
    }

    idx = 1;
    for (i = n / 2 + 1; i <= n; i++) {
      yy[i] = sqrt(yin[idx] * yin[idx] + yin[idx + 1] * yin[idx + 1]);
      idx += 2;
    }

    break;
  }
}

/************************************************************************************/
void FFTboundary_in(long int n, double *y, double *ycompl) {
  long int i;

  for (i = 1; i <= n; i++) {
    ycompl[2 * i - 1] = y[i];
    ycompl[2 * i] = 0.0;
  }
}

/************************************************************************************/

void Conv1D(double *in_v, unsigned long n, unsigned long n_conv, double delta,
            int iconv)

{
  unsigned long i, j;
  double *storage, stock1, stock2, corr_f = 1. / sqrt(8. * n),
                                   corr_f2 = 1. / sqrt(2. * n), rho;
  int k = 0; // this choses the root by varying the integer k in [0 ,n_conv-1]
             // as there are n_conv distinct roots
  storage = dvector(1, 2 * n);

  for (i = 1; i <= 2 * n; i++) {
    storage[i] = in_v[i];
  }

  switch (iconv) {

  case 1:

    for (i = 1; i <= n_conv; i++) {

      for (j = 1; j <= (2 * n) / 2; j++) {

        //					rho=pow(sqrt(storage[2*j-1]*storage[2*j-1]+storage[2*j]*storage[2*j])
        //      ,(n_conv+1)); 					stock1=
        //      rho*cos((n_conv+1)
        //*atan(storage[2*j]/storage[2*j-1]))*corr_f;
        //stock2= rho*sin((n_conv+1) *atan(storage[2*j]/storage[2*j-1]))*corr_f;

        stock1 = (storage[2 * j - 1] * in_v[2 * j - 1] -
                  storage[2 * j] * in_v[2 * j]);
        stock2 = (storage[2 * j - 1] * in_v[2 * j] +
                  storage[2 * j] * in_v[2 * j - 1]);

        storage[2 * j - 1] = stock1;
        storage[2 * j] = stock2;
      }
    }

    break;

  case -1:

    for (j = 1; j <= n; j++) {

      rho = sqrt(storage[2 * j - 1] * storage[2 * j - 1] +
                 storage[2 * j] * storage[2 * j]);
      rho = pow(rho, 1. / (n_conv + 1));

      stock1 = atan(storage[2 * j] / storage[2 * j - 1]);
      stock2 = atan(storage[2 * j] / storage[2 * j - 1]);

      stock1 = rho * cos(1. / (n_conv + 1.) * (2 * SRT_PI * k + stock1));
      stock2 = rho * sin(1. / (n_conv + 1.) * (2 * SRT_PI * k + stock2));

      storage[2 * j - 1] = stock1;
      storage[2 * j] = stock2;
    }

    break;
  }

  for (i = 1; i <= 2 * n; i++) {
    in_v[i] = storage[i];
  }

  free_dvector(storage, 1, 2 * n);
}

void ConvolveMarginals(double **AbscMarg, double **OrdMarg, unsigned long n_pts,
                       unsigned long n_assets, unsigned long n_conv,
                       int conv_flag // this is 1 if you want to convolve
                                     // and -1 if you want to deconvolve
)

{
  double *avg, *PrDens, *ycompl, *FT_PrDens, *xPrDens, Norm, delta, Mean,
      StdDev, MidPointX, MidPointY, DeltaX, InfBound, SupBound, dx, q,
      *xPrDens_new, *PrDens_new, *Xdata, *Ydata, low_thresh, sigma2;

  unsigned long i, j, k;
  unsigned int n_lagr_pts = 4, n_kill_pts = 4;

  Xdata = dvector(1, n_lagr_pts);
  Ydata = dvector(1, n_lagr_pts);

  avg = dvector(0, n_assets - 1);
  xPrDens = dvector(0, n_pts - 1);
  PrDens =
      dvector(0, n_pts - 1); // NB. n_pts must be a  1 + a power of 2 to be used
                             //     in the FFT algorithm later on
  FT_PrDens = dvector(0, n_pts - 1);
  PrDens_new = dvector(0, n_pts - 1);
  xPrDens_new = dvector(0, n_pts - 1);
  ycompl = dvector(0, 2 * (n_pts - 1));

  switch (n_conv) {
  case 0:

    break;

  default:

    for (i = 0; i < n_assets; i++) {
      avg[i] = 0.0;

      for (j = 1; j < n_pts; j++) {

        ////////////////////////////////////////////////////////////////////////////////////
        /////////// step 1: computes the average for each cumulative marginal
        ////////////////////////////////////////////////////////////////////////////////////

        avg[i] += (OrdMarg[i][j] - OrdMarg[i][j - 1]) *
                  (AbscMarg[i][j] + AbscMarg[i][j - 1]) / 2.; ///(n_conv+1.);

        ////////////////////////////////////////////////////////////////////////////////////
        /////////// step 2: computes the centered Pr. density from each
        /// cumulative marginal
        ////////////////////////////////////////////////////////////////////////////////////

        PrDens[j] = (OrdMarg[i][j] - OrdMarg[i][j - 1]) /
                    (AbscMarg[i][j] - AbscMarg[i][j - 1]);

        xPrDens[j] = (AbscMarg[i][j] + AbscMarg[i][j - 1]) / 2.;
      }

      for (j = 1; j < n_pts; j++) {
        xPrDens[j] -= avg[i];
      }

      ////////////////////////////////////////////////////////////////////////////////////
      /////////// step 3: deconvolves n_conv times each Prob. density
      ////////////////////////////////////////////////////////////////////////////////////

      for (j = 1; j <= n_conv; j++) {

        FFTboundary_in(n_pts - 1, PrDens,
                       ycompl); //  n_pts-1 must be a power a 2
        FFT1D(ycompl, n_pts - 1, 1);

        delta = 1.;
        Conv1D(ycompl, n_pts - 1, 1, delta, conv_flag);

        FFT1D(ycompl, n_pts - 1, -1);
        FFTboundary_out(n_pts - 1, ycompl, PrDens, conv_flag);

        low_thresh = PrDens[(n_pts - 1) / 2] / 40.;
        sigma2 = 1. / pow(2., j);
        for (k = 1; k < n_pts; k++) {
          if (PrDens[k] < low_thresh)
            PrDens[k] = 0.0; // exp(-xPrDens[k]*xPrDens[k]/2./sigma2);
        }

        /////////////////////////////////////////////////////
        /////  step 3b: smooths the Pr. density function
        /////
        /////////////////////////////////////////////////////

        FunctionSmooth(PrDens, n_pts - 1, 8);

        /////////////////////////////////////////////////////
        /////  step 3a: Interpolates the middle points
        /////           using the Lagrange polinomial
        /////////////////////////////////////////////////////

        for (k = 1; k <= n_lagr_pts / 2; k++) {
          Xdata[k] =
              xPrDens[(n_pts - 1) / 2 - n_kill_pts - n_lagr_pts / 2 + (k - 1)];
          Ydata[k] =
              PrDens[(n_pts - 1) / 2 - n_kill_pts - n_lagr_pts / 2 + (k - 1)];
        }
        for (k = 1; k <= n_lagr_pts / 2; k++) {
          Xdata[k + n_lagr_pts / 2] =
              xPrDens[((n_pts - 1) / 2 + 1) + (n_kill_pts + 1) + k];
          Ydata[k + n_lagr_pts / 2] =
              PrDens[((n_pts - 1) / 2 + 1) + (n_kill_pts + 1) + k];
        }

        for (k = 1; k <= 2 * n_kill_pts + 3; k++)
          PrDens[(n_pts - 1) / 2 - n_kill_pts + (k - 1)] =
              InterpLagrange1D(Xdata, Ydata, n_lagr_pts,
                               xPrDens[(n_pts - 1) / 2 - n_kill_pts + (k - 1)]);

        /////////////////////////////////////////////////////
        /////  step 3c: normalizes the deconvoluted density
        /////           and computes the standard deviation
        /////////////////////////////////////////////////////

        Norm = 0.0;
        Mean = 0.0;
        StdDev = 0.0;
        for (k = 2; k < n_pts; k++)
          Norm +=
              (xPrDens[k] - xPrDens[k - 1]) * (PrDens[k] + PrDens[k - 1]) / 2.;

        for (k = 1; k < n_pts; k++) {
          PrDens[k] = PrDens[k] / Norm;
        }

        for (k = 2; k < n_pts; k++) {
          MidPointX = (xPrDens[k] + xPrDens[k - 1]) / 2.;
          MidPointY = (PrDens[k] + PrDens[k - 1]) / 2.;
          DeltaX = (xPrDens[k] - xPrDens[k - 1]);

          Mean += MidPointX * MidPointY * DeltaX;
          StdDev += MidPointX * MidPointX * MidPointY * DeltaX;
        }

        StdDev = sqrt(StdDev - Mean * Mean);

        /////////////////////////////////////////////////////
        /////  step 3d: resets the bounds for the Pr. support
        /////           for the deconvolved density according
        /////           to the new Standard deviation
        /////////////////////////////////////////////////////

        InfBound = -4. * StdDev;
        SupBound = -InfBound;
        dx = (SupBound - InfBound) / (n_pts - 2);

        xPrDens[0] = xPrDens[1] - dx;
        PrDens[0] = PrDens[n_pts - 1];

        for (k = 1; k < n_pts; k++) {

          xPrDens_new[k] = InfBound + (k - 1) * dx;

          PrDens_new[k] = interp(xPrDens, PrDens, n_pts, xPrDens_new[k], 1, &q);
        }
        xPrDens_new[0] = xPrDens_new[1] - dx;
        PrDens_new[0] = PrDens_new[n_pts - 1];

        Norm = 0.0;
        Mean = 0.0;
        StdDev = 0.0;
        for (k = 1; k < n_pts; k++)
          Norm += (xPrDens_new[k] - xPrDens_new[k - 1]) *
                  (PrDens_new[k] + PrDens_new[k - 1]) / 2.;

        for (k = 0; k < n_pts; k++) {
          xPrDens[k] = xPrDens_new[k];
          PrDens[k] = PrDens_new[k] / Norm;
        }

        for (k = 2; k < n_pts; k++) {
          MidPointX = (xPrDens[k] + xPrDens[k - 1]) / 2.;
          MidPointY = (PrDens[k] + PrDens[k - 1]) / 2.;
          DeltaX = (xPrDens[k] - xPrDens[k - 1]);

          Mean += MidPointX * MidPointY * DeltaX;
          StdDev += MidPointX * MidPointX * MidPointY * DeltaX;
        }

        for (k = 0; k < n_pts; k++) {
          xPrDens[k] -= Mean;
        }
      }

      ////////////////////////////////////////////////////////////////////////////////////
      /////////// step 4: computes the marginal cumulative distributions
      ////////////////////////////////////////////////////////////////////////////////////

      OrdMarg[i][0] = 0.0;
      AbscMarg[i][0] = xPrDens[0] + avg[i] / (n_conv + 1.);
      for (j = 1; j < n_pts; j++) {
        AbscMarg[i][j] = xPrDens[j] + avg[i] / (n_conv + 1.);

        MidPointX = xPrDens[j] - xPrDens[j - 1];
        MidPointY = (PrDens[j] + PrDens[j - 1]) / 2.;
        OrdMarg[i][j] = OrdMarg[i][j - 1] + MidPointX * MidPointY;
      }
    }

    //		fprintf(fp2        ,"{ %f        ,  %f}        , \n"        , Mean
    //,StdDev); 	for (i=0;i<n_pts;i++) { 		fprintf(fp2        ,"{ %f        ,  %f}
    //, \n"        , AbscMarg[1][i]
    //      ,OrdMarg[1][i]);
    //	}

    free_dvector(Xdata, 1, n_lagr_pts);
    free_dvector(Ydata, 1, n_lagr_pts);
    free_dvector(avg, 0, n_assets - 1);
    free_dvector(xPrDens, 0, n_pts - 1);
    free_dvector(PrDens, 0, n_pts - 1);
    free_dvector(FT_PrDens, 0, n_pts - 1);
    free_dvector(PrDens_new, 0, n_pts - 1);
    free_dvector(xPrDens_new, 0, n_pts - 1);
    free_dvector(ycompl, 0, 2 * (n_pts - 1));

    break;
  }
}

#undef FFTSWAP
