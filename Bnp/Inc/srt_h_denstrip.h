#ifndef SRT_H_DENSTRIP_H
#define SRT_H_DENSTRIP_H
/*
 * HEADER FOR SRT_F_DENSTRIP.C
 * AUTHOR D.DESBIEZ
 */

/*
 * data structure to represent different products
 * (Futures  , Cashs  , Swaps)
 */

typedef struct deninst {
  InstType type;
  Ddate theo_start;
  Ddate theo_end;
  double rate;
  SrtBasisCode basis;
  SrtCompounding freq;

} DenInst, *DenInstPtr;

/*
 * data structure to represent Strip Curve (interest rate curve estimated
 * from futures  , cashs and swaps)
 */

#define MAXNAM 32

typedef struct denstrpcrv {
  char ccy[4];       /* currency */
  Ddate today;       /* today */
  double coeff_line; /* coefficient for smoothness */
  int num_iter;      /* Number of iterations for levenberg */
  int algo_flag;     /* 0 = LINEAR  , 1 = SPLINE */
  int calib_flag;    /* 0 = LEVENBERG  , 1 = ANNEALING */

  Ddate lastdate;         /* last date that the stripper can interpolate */
  int numknot;            /* number of knot points */
  int lastknot;           /* index of the last optimised knot */
  Ddate *knot_date;       /* date of knot points  , indexed from 1 */
  double *knot_point;     /* time in years from today of knot points */
                          /* indexed from 1 */
  DenInst *inst;          /* Instrument used to optimize over */
  int numinst;            /* Number of instruments indexed from 1 */
  double *fwd_rate_coef;  /* coef of fwd rate  , indexed from 1 */
  double **zcpn_yld_coef; /* coef of the integral of fwd rate indexed from 0 */
                          /* for each knot intervals indexed from [0  ,0] */
} DenStrpCrv, *DenStrpCrvPtr;

typedef struct denstrpcrvlist {
  DenStrpCrv *crv;
  char crv_name[MAXNAM];
  long ticker;

} DenStrpCrvList, *DenStrpCrvListPtr;

/***************/
/** FUNCTIONS **/
/***************/

/* discount factor */
SrtErr srt_f_dencrv_df(DenStrpCrv *crv, Ddate start, Ddate end, double *answer);

/* zero rate */
SrtErr srt_f_dencrv_zr(DenStrpCrv *crv, Ddate start, Ddate end,
                       SrtBasisCode basis, double *answer);

/* instantaneous forward rate */
SrtErr srt_f_dencrv_fwdr(DenStrpCrv *crv, Ddate start, double *answer);

/* To price a swap */
double srt_f_den_swappr(Ddate start, Ddate end, SrtCompounding freq,
                        SrtBasisCode basis, DenStrpCrv *crv);

/* To price a deposit */
double srt_f_den_depositpr(Ddate start, Ddate end, SrtBasisCode basis,
                           DenStrpCrv *crv);

/* return today */
Ddate srt_f_dencrv_today(DenStrpCrv *crv);

/* return lastdate */
Ddate srt_f_dencrv_lastdate(DenStrpCrv *crv);

/* free */
void srt_f_dencrv_free(DenStrpCrv **crv);

/* allocate and estimate curve */
SrtErr srt_f_dencrv_init(
    DenStrpCrv **crvptr, /* Return of the strip curve */
    DenInst *inst,       /* array containing information on instruments */
    int numinst,         /* number of instruments */
    Ddate today,         /* today's date */
    char *ccy,           /* currency */
    double coeff_line,   /* coefficient for smoothness */
    int num_iter,        /* Number of iterations for levenberg */
    int algo_flag,       /* 0 = LINEAR  , 1 = SPLINE */
    int calib_flag);     /* 0 = LEVENBERG  , 1 = ANNEALING */

#endif
