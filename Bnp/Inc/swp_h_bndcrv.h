#ifndef SRT_H_BNDCRV_H
#define SRT_H_BNDCRV_H
/*
 * HEADER FOR SRT_F_BNDCRV.C
 * AUTHOR E.AULD
 * MATH R.Roy
 */

/*
 * data structure to represent Bond Curve (interest rate curve estimated
 * from bond prices)
 */

#define MAXNAM 32

typedef struct srtbndcrv
{
    int      algo_flag;      /* algo flag 0 constant, 1 linear, 2 spline */
    char     ccy[4];         /* currency */
    Ddate    today;          /* today */
    Ddate    spot;           /* spot */
    Ddate    lastdate;       /* last date that the stripper can interpolate */
    int      num_knot_point; /* number of knot points */
    int      lastknot;       /* index of the last optimised knot */
    Ddate*   knot_date;      /* date of knot points, indexed from 1 */
    double*  knot_point;     /* time in years from today of knot points, indexed from 1 */
    double*  fwd_rate_coef;  /* coef of fwd rate, indexed from 0 */
    double** zcpn_yld_coef;  /* coef of the integral of fwd rate indexed from 0 */
                             /* for each knot intervals indexed from 0,0 */
} SrtBndCrv, *SrtBndCrvPtr;

typedef struct srtbndcrvlist
{
    SrtBndCrv* crv;
    char       crv_name[MAXNAM];
    long       ticker;

} SrtBndCrvList, *SrtBndCrvListPtr;

/*
 * externally supported functions using srt bnd curves
 * (disc_mkt will be extended to call these)
 */

/* discount factor */
SrtErr srt_f_bndcrv_df(void* crv, Ddate start, Ddate end, double* answer);

/* level rate */
SrtErr srt_f_bndcrv_level(SrtBndCrv* crv, SwapDP* sdp, double* answer);

/* zero rate */
SrtErr srt_f_bndcrv_zr(SrtBndCrv* crv, Ddate start, Ddate end, SrtBasisCode basis, double* answer);

/* instantaneous forward rate */
SrtErr srt_f_bndcrv_fwdr(SrtBndCrv* crv, Ddate start, double* answer);

/* swap rate/fwd treasury rate */
SrtErr srt_f_bndcrv_fwdtr(SrtBndCrv* crv, Date s, Date e, char* f, char* b, double* answer);

/* fwd bond price */
SrtErr srt_f_bndcrv_fwdpr(
    SrtBndCrv* crv, SwapDP* sdp, double cpn, double redemption, double* answer);

/* fwd bond price or yield */
SrtErr srt_f_bndcrv_fwdyld(
    SrtBndCrv* crv, SwapDP* sdp, double cpn, double redemption, double* answer);

/* return today */
Date srt_f_bndcrv_today(SrtBndCrv* crv);

/* return lastdate */
Date srt_f_bndcrv_lastdate(SrtBndCrv* crv);

/* free */
void srt_f_bndcrv_free(SrtBndCrv** crv);

/* allocate and estimate curve */
SrtErr srt_f_bndcrv_init(
    SrtBndCrv**     crvptr,        /* Return of the bond curve */
    SrtBasisCode*   basis,         /* array containing basis of bonds */
    SrtCompounding* compounding,   /* array containing compounding of bonds */
    Ddate*          prev_cpn_date, /* array containing previous coupon date of bonds */
    int*            num_coupons,   /* array containing number of coupon in bonds */
    Ddate**         coupon_dates,  /* array containing coupon dates of bonds*/
    double**        coupons,       /* array containing coupons of bonds*/
    double*         weight,        /* array containing weight (standart deviation) of bonds */
    double*         fwd_price,     /* array containing the forward clean price */
    double*         redemption,    /* array containing the redemption */
    int             numbond,       /* number of bonds */
    double*         knot_point,    /* array containing the knot points (in times) */
    int             numknot,       /* number of knot points */
    Ddate           today,         /* today's date */
    double          coeff_line,    /* coefficient for smoothness */
    int             num_iter,      /* Number of iterations for levenberg */
    int             algo_flag);                /* 0 = LINEAR, 1 = SPLINE */

#endif
