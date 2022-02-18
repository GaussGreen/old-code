/* ===========================================================================

         FILENAME:     swp_h_df.h

     PURPOSE:      The function that will have to be used everywhere to compute
                       a discount factor or a FwdCash
                                   They use  a pointer to a function DiscFunc
                   passed to the SrtExternalFunctions static
                                   ( see srt_f_external_fct.c )

   =========================================================================== */

#ifndef SWP_H_DISC_CURVE_H
#define SWP_H_DISC_CURVE_H

#ifndef SRT_DF_ERROR
#define SRT_DF_ERROR DBL_MAX
#endif

double swp_f_disc_or_zrate(Ddate dstart, Ddate dend, char* crv_name, int rate_or_disc);

double swp_f_fwdcash(Ddate start, Ddate end, BasisCode basis, char* crv_name);

double swp_f_fwdcashrate(Ddate fixing, Ddate start, Ddate end, BasisCode b, char* crvname);

double swp_f_swapcashrate(
    long           start,
    long           end,
    SrtBasisCode   Basis,
    SrtCompounding Comp,
    String         ycName,
    String         refRateCode);

#define swp_f_df(start, end, name) swp_f_disc_or_zrate((Ddate)start, (Ddate)end, (char*)name, 1)

#define swp_f_zr(start, end, name) swp_f_disc_or_zrate((Ddate)start, (Ddate)end, (char*)name, 0)

#endif