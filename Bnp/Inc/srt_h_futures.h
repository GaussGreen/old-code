/* =======================================================================================
   FILENAME :       srt_h_futures.h
   ======================================================================================= */
   
#ifndef SRT_H_FUTURES_H
#define SRT_H_FUTURES_H

extern double*  dDates;
extern double*  dVolStruct;
extern double*  dVolStruct2;
extern double** ppdFwdVol;
extern double** ppdOptFutVol;
extern double*  dB0T;
extern double*  dL0T;
extern double*  dFwd;
extern double*  dCoverage;
extern int      iNumFut;
extern int      nbTenors2;

Err srt_f_futures(long     lToday,
                  int      iType,
                  long     *lFutDates,
                  int      numFut,
                  double   *dParams,
                  double   *futuresPrice,
									double   *convexity,
                  double   *rate);

Err srt_f_calibfutures(long     lToday,
                       long     *lFutDates,
                       int      iNumFut,
                       int      iType,
                       double   *alpha,
                       double   *beta,
												double	*convexities,
                       double   *error);
#endif //SRT_H_FUTURES_H
