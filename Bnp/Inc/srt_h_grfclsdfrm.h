#ifndef GRF_H_CLSDFRM_H
#define GRF_H_CLSDFRM_H

#include "grf_h_all.h"
#include "srt_h_all.h"

/* srt_f_grfclsdfrm.c */

/* -------------------------------------------------------------------------- */

Err grfn_SwapDParray_to_GrfnCells(long *numeventdatesptr, Date **eventdatesptr,
                                  long *nrowsptr, long *ncolsptr,
                                  GrfnCell ***sprdshtptr, Date today,
                                  int lNumInstr, SwapDP *sdp, double *strike,
                                  double *bk, SrtReceiverType *rec_pay,
                                  StructType *type, String und_name,
                                  String *ref_rate_code);

Err new_grfn_SwapDParray_to_GrfnCells(
    long *numeventdatesptr, Date **eventdatesptr, long *nrowsptr,
    long *ncolsptr, GrfnCell ***sprdshtptr, Date today, int lNumInstr,
    SwapDP *sdp, double *strike, double *bk, SrtReceiverType *rec_pay,
    StructType *type, String und_name, String *ref_rate_code);

/* -------------------------------------------------------------------------- */

Err grfn_make_swaption_string(String swptn_str, SwapDP *sdp, double strike,
                              SrtReceiverType rec_pay, String und_name,
                              String ref_rate);

/* -------------------------------------------------------------------------- */

Err grfn_make_caplet_string(String caplet_str, SwapDP *sdp, double strike,
                            SrtReceiverType rec_pay, String und_name,
                            String ref_rate);

/* -------------------------------------------------------------------------- */

Err grfn_SwapDP_to_GrfnCells(long *numeventdatesptr, Date **eventdatesptr,
                             long *nrowsptr, long *ncolsptr,
                             GrfnCell ***sprdshtptr, Date today, SwapDP *sdp,
                             double strike, double bk, SrtReceiverType rec_pay,
                             StructType type, String und_name,
                             String ref_rate_code);

/* -------------------------------------------------------------------------- */

Err grfn_SwapDP_to_BDT_GrfnCells(long *numeventdatesptr, Date **eventdatesptr,
                                 long *nrowsptr, long *ncolsptr,
                                 GrfnCell ***sprdshtptr, SrtUndPtr und,
                                 SwapDP *sdp, double strike, double bk,
                                 SrtReceiverType rec_pay, StructType t);

/* -------------------------------------------------------------------------- */

/*****************************************************************************
   FUNCTION		:		GRFN_CLSDFRM_VEGA
   DESCRIPTION	:		COMPUTES VEGAS FOR GRFN_CLSDFRM
******************************************************************************/

Err srt_f_grfn_clsdfrm_vega(SrtUndPtr und, SrtGrfnParam *grfnparam, SwapDP *sdp,
                            double strike, double bk, SrtReceiverType rec_pay,
                            StructType type, String ref_rate_code,
                            double *initial_price, double **sigma_vega,
                            long *num_sig, double **tau_vega, long *num_tau);

/* ---------------------------------------------------------------------------
 */

/* srt_f_grf_f_clsdfrm.c */
Err srt_f_grfn_clsdfrm(SrtUndPtr und, SrtGrfnParam *grfparam, SwapDP *sdp,
                       double strike, double bk, SrtReceiverType rec_pay,
                       StructType type, String ref_rate_code, double *answer);

/* srt_f_grf_f_newclsdfrm.c */
Err srt_f_grfn_newclsdfrm(SrtUndPtr und, SrtGrfnParam *grfnparam, int lNumInstr,
                          SwapDP *sdp, double *strike, double *bk,
                          SrtReceiverType *rec_pay, StructType *type,
                          String *ref_rate_code, double *answer);

/* This Function is in srt_f_grfnfutureclosedform */
Err srt_f_grfn_future_closedform(long future_date, SrtUndPtr und,
                                 SrtGrfnParam *grfnparam, SwapDP *sdp,
                                 double strike, double bondstrike,
                                 SrtReceiverType rec_pay, StructType type,
                                 String ref_rate_code, double *price);
/* ---------------------------------------------------------------------------
 */

/*** grf_f_midat_clsdfrm.c ***/
Err grfn_midat_clsdfrm(
    long num_exercise_dates,    /* len of next two arrays */
    Date *exercise_dates,       /* dates when option can be exercise_d */
    Date *exercise_start_dates, /* dates when bonds start      ,>= correponding
                          exercise dates. */
    double *
        exercise_premiums, /* amounts that must be payed to exercise_ options */
    long num_prod_dates, Date *prod_dates, double *prod_cfs,
    SrtReceiverType rec_pay, SrtUndPtr und, /* name of und to use */
    SrtGrfnParam *grfnparam,                /* model and implementation
                        details */
    double *answer                          /* value returned here*/
);

/*** srt_f_grfcntcap ***/
Err grf_cntcap_clsdfrm(int compd, int basis, long num_dates,
                       /* len of next two arrays */
                       Date *dates,
                       /* dates when option can be exercise_d */
                       double *strikes, double *contingent_premiums,
                       /* amounts to be paid if strike is reached */
                       SrtReceiverType cap_floor,
                       /* cap or floor */
                       SrtUndPtr und,
                       /* name of und to use */
                       SrtGrfnParam *grfnparam,
                       /* model and implementation details */
                       double *answer
                       /* value returned */
);

/*** srt_f_grfameswp ***/
Err grf_ameswp_clsdfrm(double start, double end, int nfp, int delay, int compd,
                       int basis, double *strikes, SrtReceiverType rec_pay,
                       /* receiver or payer */
                       SrtUndPtr und,
                       /* name of und to use */
                       SrtGrfnParam *grfnparam,
                       /* model and implementation details */
                       double *answer
                       /* value returned */
);

#endif
