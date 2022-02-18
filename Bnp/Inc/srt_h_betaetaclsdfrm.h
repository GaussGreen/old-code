#ifndef SRT_H_BETAETACLSDFRM_H
#define SRT_H_BETAETACLSDFRM_H

double srt_f_etabeta_capfloor(
    TermStruct*     ts,
    double*         fixing_time,
    double*         pay_time, /* pay_time[0] is the strike payment time */
    double*         df,
    double*         coupon,
    int             num_caplets,
    SrtReceiverType rec_pay,
    SrtMdlDim       mdl_dim);

double srt_f_etabeta_coupon_bond_option(
    int             n,
    double          bond_strike,
    TermStruct*     ts,
    double          fixing_time,
    double*         coupon,
    double*         pay_time,
    double*         df,
    SrtReceiverType rec_pay,
    SrtMdlDim       mdl_dim);

#endif
