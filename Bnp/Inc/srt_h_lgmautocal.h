#ifndef SRT_LGMAUTOCAL_H
#define SRT_LGMAUTOCAL_H

#include "srt_h_lgmtypes.h"

/* This holds data for a single exercise boundary swap */
/* see srt_h_lgmtypes.h */

/* This holds all exercise boundary swap data */
/* see srt_h_lgmtypes.h */

/* This holds data for a single LGM reference swaption */
/* see srt_h_lgmtypes.h */

/* This holds all LGM reference swaption data */
/* see srt_h_lgmtypes.h */

/* This holds LGM tau/sigma data for a single point */
/* see srt_h_lgmtypes.h */

/* This holds all LGM tau/sigma data */
/* see srt_h_lgmtypes.h */

String LGMamerican(
    long   nfix,           /* number of fixed periods                                  */
    Date   tfix_start[],   /* [0,1,...,nfix-1]  start dates for fixed coupons          */
    Date   tfix_end[],     /* [0,1,...,nfix-1]  end dates for fixed coupons            */
    Date   tfix_pay[],     /* [0,1,...,nfix-1]  pay dates for fixed coupons            */
    double full_pay_fix[], /* [0,1,...,nfix-1]  fixed coupons (with notional at end)   */
    double Inprem[],       /* [0,1,...,nfix-1]  premium for exer. at j (with notional) */
    char*  fix_basis,      /* basis for fixed coupons                                  */
    int    early_flag_fix, /* 0=subtract accrual from frst pymnt; 1=add accrual to fee */
    long   nflt,           /* number of floating dates                                 */
    Date   t_flt[],        /* [0,1,...,nflt-1] all flting dtes (all strt dtes plus final end dte) */
                           /* BFA: remove this line and uncomment the next line to implement */
    /* double   rflt_current */ /* current floating rate (if fixed); ignored if non-positive   */
    char*      flt_basis,       /* basis for floating coupons                                   */
    int        early_flag_flt,  /* 0 = subtrct accrual from frst pymnt; 1 = add accrual to prem */
    Date       t_first_ex,      /* first exerercise date                                        */
    int        lag_exer_start,  /* days between exercise and start                              */
    int        cal_or_bus,      /* 0 = cal. days, 1 = bus. days for lag_exer_start              */
    BusDayConv conv_start,      /* business day convention for start                            */
    char*      pay_rec_str,     /* SRT_RECEIVER or SRT_PAYER                                */
    String     yc_name,         /* Yield curve name    */
    int        end_of_day_flag, /*      */
    int        fix_tau_flag,    /* 1 means calibrate with fixed tau */
    int        usecaps, /* if fixed tau, then 1=use caplets to calibrate, 0=use long swptns */
    double     In_tau,  /* if fixed tau, use this value for tau (in years) */
    String (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* function to get swaption vols */
    char*   char_vol_type,                                      /*       */
    int     update_ts_flag, /* 1=cmpute new sigmas and taus & store in yldcrv; 0=don't bother */
    int     find_exer_bdry, /* 1=find swap rates at exercise bocrvary; 0=don't bother       */
    String  outfile,        /* output log file (unused)                                     */
    double* LGMamervalue,   /* answer                                                       */
    SrtLgmExerBdryData* lgmExerBdryData); /* ptr to exercise boundary data structure         */

String LGMautocal(
    long   In_ex,           /* J = In_ex-1; In_ex is # of exercises                         */
    Date   It_ex[],         /* notification (exercise) dates   T[0,1,...,J]                 */
    Date   It_start[],      /* start dates for exercise Tj                     s[0,1,...,J] */
    double Iprem[],         /* premium (inc. notional) paid to exercise at Tj  f[0,1,...,J] */
    long   npay,            /* I = npay-1; npay is # of fixed leg coupons                   */
    Date   tfix_start[],    /* start dates for coupon i                   tst[0,1, ..., I]  */
    Date   tfix_end[],      /* end dates for coupon i                     tend[0,1, ..., I] */
    Date   tpay[],          /* pay dates for coupon i                     tpay[0,1, ..., I] */
    double fixpymnt[],      /* coupons (incl. notional in last)           fixpymnt[0,...,I] */
    double Iredpay[],       /* reduction in 1rst payment after exercise   Iredpay[0,..J]    */
    char*  pay_rec_str,     /* SRT_RECEIVER or SRT_PAYER                                    */
    String yc_name,         /* pointer to market structures                                 */
    int    end_of_day_flag, /* 0: intra day 1: end of day */
    int    fix_tau_flag,    /* 1 means calibrate with fixed tau */
    int    usecaps,         /* if fixed tau, then 1=use caplets to calibrate, 0=use long swptns */
    double In_tau,          /* if fixed tau, use this value for tau (in years) */
    String (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* function to get swaption vols */
    char*   char_vol_type,
    int     update_ts_flag, /* 1=compte new sigs and taus & store in yldcrv; 0=don't bother */
    int     find_exer_bdry, /* 1 = find swap rates at exercise bocrvary; 0 = don't bother  */
    String  outfile,        /* output file name for log file (unused)                      */
    double* LGMvalue,       /* LGM value of mid-atlantic             */
    SrtLgmExerBdryData*
        lgmExerBdryData, /* ptr to exercise boundary data structure (NULL => not req'd) */
    SrtLgmRefSwptnData*
        lgmRefSwptnData,          /* ptr to reference swaption data structure (NULL => not req'd) */
    SrtLgmTSData* lgmTSData,      /* ptr to tau/sigma data (NULL => not req'd) */
    double*       intrinsicValue, /* ptr to intrinsic value (NULL => not req'd) */
    int skip_deal_evaluation);    /* 0 for a noraml valuation, 1 for a simple volatility points
                                     extraction */

#endif /* if SRT_LGMAUTOCAL_H not defined */
