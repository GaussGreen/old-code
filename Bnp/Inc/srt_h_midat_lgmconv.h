#ifndef SRT_MIDAT_LGMCONV_H
#define SRT_MIDAT_LGMCONV_H

Err srt_f_lgm_midat_conv_eval(
    Date            tnow,
    long            n,
    SrtUndPtr       und,
    long            nex,
    Date            t_ex[],
    Date            t_start[],
    long            npay,
    Date            tpay[],
    double          pvpay[],
    double          redpay[],
    long            ifirst[],
    Date            dates[],
    double          prem[],
    SrtReceiverType rec_pay,
    double*         LGMvalue,
    double          phi_cr[]);

Err Findexerbdry(
    Date           today,
    long           nex,
    Date           t_ex[],
    Date           tend,
    Date           dates[],
    double         g[],
    double         zeta[],
    double         phi_cr[],
    long           n,
    SrtUndPtr      und,
    SrtCompounding compounding,
    BasisCode      natural_basis,
    BusDayConv     natural_conv,
    int            natural_spot_lag,
    long*          n_exer_bdryptr,
    double*        r_exer_bdryptr[],
    Date*          t_exer_bdryptr[]);

#endif
