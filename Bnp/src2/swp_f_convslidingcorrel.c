/* ===================================================================================
   FILENAME:      swp_f_convslidingcorrel.c

   PURPOSE:       Computes the correlation between two fwd swap rates
                  (converging sliding method)
                                  Inputs are spot swaps correlation matrix and
   tenor dates associated
   ===================================================================================
 */

#pragma warning(disable : 4786) // Disable long name warnings

#include "math.h"
#include "swp_h_convslidingcorrel.h"
#include "swp_h_vol.h"

Err Compute_CoInitalSwaps_Correl2(
    long lToday, char *CorrelCubeId,
    Err (*get_correl)(char *correl_cube_name, double start_date,
                      double end_date, double strike, double *vol),
    long lStartDate, long lEndDate, SrtCompounding srtFreq,
    SrtBasisCode srtBasis, int NDimCorrel, double ***dCorrelMatrix) {
  Err err = NULL;
  int i, j, iNDates;

  double rho, strike;
  long start, end;

  iNDates = NDimCorrel;
  (*dCorrelMatrix) = dmatrix(0, iNDates - 1, 0, iNDates - 1);
  strike = coverage(lToday, lStartDate, BASIS_ACT_365) / 100.0;
  for (i = 0; i < iNDates; ++i) {
    (*dCorrelMatrix)[i][i] = 1.0;
    for (j = i + 1; j < iNDates; ++j) {
      start =
          add_unit(lToday, i * 12 / srtFreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
      end = add_unit((long)(start), (long)(fabs(i - j) * 12 / srtFreq),
                     SRT_MONTH, NO_BUSDAY_CONVENTION);
      err = get_correl(CorrelCubeId, start, end, strike, &rho);
      (*dCorrelMatrix)[i][j] = DMIN(0.99999, DMAX(-0.99999, rho));
      (*dCorrelMatrix)[j][i] = (*dCorrelMatrix)[i][j];
    }
  }

  return err;
}

Err Compute_CoInitalSwaps_Correl(
    long lToday, char *CorrelCubeId,
    Err (*get_correl)(char *correl_cube_name, double start_date,
                      double end_date, double strike, double *vol),
    long lStartDate, long lEndDate, SrtCompounding srtFreq,
    SrtBasisCode srtBasis, int *NDimCorrel, double ***dCorrelMatrix) {
  SwapDP Swap;
  Err err = NULL;
  long *lPayDates = NULL;
  long *lStartDates = NULL;
  long *lEndDates = NULL;
  double *dCoverages = NULL;
  int i, j, iNPayDates, iNDates;

  double rho, strike;
  long start, end;

  err = swp_f_setSwapDP(lStartDate, lEndDate, srtFreq, srtBasis, &Swap);
  if (err) {
    goto FREE_RETURN;
  }

  err = swp_f_make_FixedLegDatesAndCoverages(&Swap, lToday, &lPayDates,
                                             &iNPayDates, &lStartDates,
                                             &lEndDates, &dCoverages, &iNDates);
  if (err) {
    goto FREE_RETURN;
  }

  *NDimCorrel = iNDates;
  (*dCorrelMatrix) = dmatrix(0, iNDates - 1, 0, iNDates - 1);
  strike = coverage(lToday, lStartDate, BASIS_ACT_365) / 100.0;
  for (i = 0; i < iNDates; ++i) {
    (*dCorrelMatrix)[i][i] = 1.0;
    for (j = i + 1; j < iNDates; ++j) {
      start =
          add_unit(lToday, i * 12 / srtFreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
      end = add_unit((long)(start), (long)(fabs(i - j) * 12 / srtFreq),
                     SRT_MONTH, NO_BUSDAY_CONVENTION);
      //			err = swp_f_vol( CorrelCubeId  , start  , end  ,
      //strike  , &rho  , &pow);
      err = get_correl(CorrelCubeId, start, end, strike, &rho);
      (*dCorrelMatrix)[i][j] = DMIN(0.99999, DMAX(-0.99999, rho));
      (*dCorrelMatrix)[j][i] = (*dCorrelMatrix)[i][j];
    }
  }

FREE_RETURN:

  if (lPayDates)
    free(lPayDates);
  if (lStartDates)
    free(lStartDates);
  if (lEndDates)
    free(lEndDates);
  if (dCoverages)
    free(dCoverages);

  return err;
}
