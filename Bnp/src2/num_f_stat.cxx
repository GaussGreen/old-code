/* ===============================================================
   FILE_NAME:	num_f_stat.cxx

   PURPOSE:     A few functions to do minimal statistics
   =============================================================== */
#include "math.h"
#include "num_h_allhdr.h"

static double avef(double *x, int len) {
  double sum = 0.0;
  int i;
  for (i = 0; i < len; i++) {
    sum += x[i];
  }
  return (sum / ((double)len));
}

static double inner(double *x, double *y, int len) {
  int i;
  double inner = 0.0;
  for (i = 0; i < len; i++) {
    inner += x[i] * y[i];
  }
  return inner;
}
static void dif(double *x, int len, int flg) {
  int i;
  if (flg == SRT_LOGNORMAL) {
    for (i = 0; i < len - 1; i++) {
      x[i] = log(x[i + 1] / x[i]);
    }
  } else {
    for (i = 0; i < len - 1; i++) {
      x[i] = x[i + 1] - x[i];
    }
  }
}

/* public functions in this module */

/*
 * returns stdev of data        , or of its returns or differences depending:
 * if flg=0 returns stdev of returns
 * if flg=1 returns stdev of first differences
 * if flg=2 returns stdev
 * if aveptr is not NULL it is set equal to
 * the average/avg difference/avg return. Almost no error checking
 * is performed.
 */

double srt_f_stdev(double *data1, int len, int flg, double *aveptr) {
  double vol, ave;

  if (flg == SRT_LOGNORMAL || flg == SRT_NORMAL) {
    dif(data1, len, flg);
    len--;
  }
  vol = inner(data1, data1, len);
  ave = avef(data1, len);
  vol = sqrt(vol / (len - 1.0) - ave * ave * len / (len - 1.0));
  if (aveptr) {
    *aveptr = ave;
  }
  return vol;
}

/*
 * computes correlation between data1 and data2.  These are first transformed
 * by calling srt_f_stdev.  It is assumed that the lengths are the same;
 * e.g. this will break if you choose flg1 = 0 or 1 and flg2 = 2. (for example)
 */
double srt_f_corr(double *data1, double *data2, int len, int flg1, int flg2) {
  double ave1, ave2, volx, voly, corr;
  volx = srt_f_stdev(data1, len, flg1, &ave1);
  voly = srt_f_stdev(data2, len, flg2, &ave2);
  corr = inner(data1, data2, len - 1);
  corr = corr / (len - 2.0) - ave1 * ave2 * (len - 1.0) / (len - 2.0);
  corr = corr / (volx * voly);
  return corr;
}
