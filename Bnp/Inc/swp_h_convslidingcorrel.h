/* ===================================================================================
   FILENAME:      swp_h_convslidingcorrel.h

   PURPOSE:       Computes correlation (converging sliding method)
   ===================================================================================
 */

#ifndef SWP_H_CONVSLIDINGCORREL_H
#define SWP_H_CONVSLIDINGCORREL_H

#include "swp_h_all.h"

double converging_sliding_correlation(long lToday, SrtCompounding srtFreqCorrel,
                                      SrtBasisCode srtBasis, int nTenorDates,
                                      long *lTenorDates, double **CorrelMatrix,
                                      long lStartDateCorrel,
                                      long lEndDateCorrel, long lStart1,
                                      long lEnd1, long lStart2, long lEnd2);

Err Compute_CoInitalSwaps_Correl(
    long lToday, char *CorrelCubeId,
    Err (*get_correl)(char *correl_cube_name, double start_date,
                      double end_date, double strike, double *vol),
    long lStartDate, long lEndDate, SrtCompounding srtFreq,
    SrtBasisCode srtBasis, int *NDimCorrel, double ***dCorrelMatrix);

Err Compute_CoInitalSwaps_Correl2(
    long lToday, char *CorrelCubeId,
    Err (*get_correl)(char *correl_cube_name, double start_date,
                      double end_date, double strike, double *vol),
    long lStartDate, long lEndDate, SrtCompounding srtFreq,
    SrtBasisCode srtBasis, int NDimCorrel, double ***dCorrelMatrix);

#endif