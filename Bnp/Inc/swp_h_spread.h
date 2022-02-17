/* ===========================================================================

         FILENAME:     swp_h_spread.h

     PURPOSE:      The function that will have to be used everywhere to compute
                       a spread for a given reference rate wrt Cash

   ===========================================================================
 */

#ifndef SWP_H_SPREAD_H
#define SWP_H_SPREAD_H

#ifndef SRT_SPREAD_ERROR
#define SRT_SPREAD_ERROR DBL_MAX
#endif

double swp_f_spread(Ddate start, Ddate end, String ref_rate_name);

Err swp_f_get_ref_rate_details(String ref_rate_name, SrtBasisCode *float_basis,
                               SrtCompounding *float_compounding);

Err srt_f_get_spot_lag_from_refrate(String ref_rate_name, int *spot_lag);

double swp_f_fra(Ddate start, Ddate end, BasisCode basis, String yc_name,
                 String ref_rate);
#endif