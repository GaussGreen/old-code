/**********************************************************************
 *      Name: srt_f_rngfltr.c                                         *
 *  Function: Functions to compute fraction of days accrued for       *
 *            accrual swaps and range floaters and resettable 		  *
 *            floaters.                                               *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Adam Litke  , Rishin Roy  , Eric Auld                       *
 *      Date: 10/02/95                                                *
 *--------------------------------------------------------------------*
 *    Inputs:                                                         *
 *   Returns:                                                         *
 *   Globals:                                                         *
 *                                                                    *
 * FYI
 ** (1)accrual swap:you get the fixed coupon times the number of	  * days
 *libor is within a range
 ** (2)range floater: you get the floating coupon fixed at t1 times	  * the
 *number of days libor is within a range
 ** (3)resettable: you get the above where you get to * fix the band at t1.
 **
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 10/02/95 AL      Started from acc_addin.c
 **
 **********************************************************************/

#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "num_h_interp.h"
#include "num_h_proba.h"
#include "srt_h_rngfltr.h"
#include "swp_h_vol_interpol.h"
#include "utallhdr.h"

/*
  suppose y = f(x)  , fit a quadratic (there is only one) through these
  points  , and return the value x that min/maximizes that quadratic
  we find the quadratic using lagranges formula
*/
static double quadratic_minmax(double x[3], double y[3]) {
  double xminmax;
  double x12, x10, y12, y10;

  x12 = x[1] - x[2];
  x10 = x[1] - x[0];
  y12 = y[1] - y[2];
  y10 = y[1] - y[0];

  if (fabs(y10) + fabs(y12) < EPS)
    return x[1];

  xminmax = 2.0 * (x10 * y12 - x12 * y10);
  xminmax = x[1] - (x10 * x10 * y12 - x12 * x12 * y10) / xminmax;
  return xminmax;
}

/* WE MUST ROUND KA TO THE NEAREST 32ND IN THE RIGHT DIRECTION */
static double nearest_32nd(double strike, double cf, double wl) {
  double ka;
  /* FIRST COMPUTE  THE NUMBER OF 16THS THAT IS <=K*/
  ka = (double)(DTOL(strike * 1600.0));
  /* NOW WE CHECK THAT KA IS REALLY LESS THAN NOT <= */
  if (ka < strike * 1600.0 - EPS) {
    ka = 2.0 * ka + 1.0;
    ka /= 3200.0;
  } else /* KA==strike */
  {
    ka = 2.0 * ka + cf * wl;
    ka /= 3200.0;
  }
  return ka;
}

static double get_vol(double vol_date, double strike, double *mat_dates_vec,
                      double *mkt_vol_1darray, double *strikes_vec,
                      double **mkt_vol_2darray, long m_mat_dates,
                      long m_strikes, int method, int smile_flg) {

  double vol, q;

#ifdef RFLTR_DBG_PRINT
  printf("getvol: vol_date %lf strike %lf\n", vol_date, strike);
#endif

  if (smile_flg) {
    vol = vol_interpol_function(vol_date, strike, mat_dates_vec, m_mat_dates,
                                strikes_vec, m_strikes, mkt_vol_2darray,
                                m_mat_dates, m_strikes);
  } else {
    vol = interp(mat_dates_vec, mkt_vol_1darray, m_mat_dates, vol_date, method,
                 &q);
    vol = q;
  }
#ifdef RFLTR_DBG_PRINT
  printf("getvol: vol %lf strike %lf\n", vol);
#endif
  return vol;
}

/*
  CALCULATE D2 FOR THE MULTITUDE OF DIFFERENT CASES WE NOW ALLOW:
*/

/* x is libor we are going to receive  , y is libor we are contingent on*/
/* float_flg determines if we are going to receive x or k (if we are going
  to receive x  , we have to make a change of numeraire so there will be an
  adjustment); reset_flg determines whether or not the band is resettable;
  if it is then receiving will be contingent on y>x+k (for a cap); we value
  this by making a margrabe type assumption that x+k or y-k is lognormal.
*/

/*
get_d2(fwdy  ,fwdx  ,k  ,fvol  ,vol  ,float_set  ,y_set  ,today  ,correlation
,float_tenor  , float_flg  ,reset_flg);
*/
static double get_d2(double fwdy, double fwdx, double k, double volx,
                     double voly, Date x_set, Date y_set, Date today, double p,
                     double float_tenor, int float_flg, int reset_flg) {
  double d2;
  double adj = 0;
  double tx, ty;
  double var;

  tx = (x_set - today) * YEARS_IN_DAY;
  ty = (y_set - today) * YEARS_IN_DAY;
  p = DMAX(1.0 - (1.0 - p) * DMAX(ty - tx, 0.0) * 12.0 / float_tenor, p);

#ifdef RFLTR_DBG_PRINT
  printf("getd2: fwdy %lf fwdx %lf k %lf volx %lf voly %lf\n", fwdy, fwdx, k,
         volx, voly);
#endif

  if (!float_flg && !reset_flg) {
    d2 = (log(fwdy / k) - .5 * voly * voly * ty) / (voly * sqrt(ty));
  } else if (!reset_flg) {
    d2 = (voly * sqrt(ty));
    d2 = (log(fwdy / k) - .5 * voly * voly * ty) / d2;
    adj = p * volx * tx / sqrt(ty);
    d2 += adj;
  } else /* contingent on Y(t2) > X(t1)+k */
  {
    if (k >= 0) /* say x + k is lognormal */
    {
      volx *= fwdx / (fwdx + k);
      fwdx += k;
    } else /* say y - k is lognormal */
    {
      voly *= fwdy / (fwdy - k);
      fwdy -= k;
    }
    var = sqrt(volx * volx * tx + voly * voly * ty - 2 * p * volx * voly * tx);
    if (float_flg) {
      d2 = (log(fwdy / fwdx) - 0.5 * var * var) / var;
    } else {
      d2 =
          (log(fwdy / fwdx) - 0.5 * voly * voly * ty + 0.5 * volx * volx * tx) /
          var;
    }
  }
#ifdef RFLTR_DBG_PRINT
  printf("getd2: d2 %lf\n", d2);
#endif

  return d2;
}

/*
   THIS FUNCTION VALUES THE FRACTION OF DAYS ACCRUED FOR THE FLOATING PIECE OF
   RANGE FLOATERS AND THE FIXED DAYS ACCRUED FOR ACCRUAL SWAPS AND THE
   SPREAD PIECE OF RANGE FLOATERS.  THIS FUNCTION ONLY VALUES ONE "SIDE"
   --- RETURNS NUMBER OF DAYS LIKELY TO BE ABOVE STRIKE FOR CAP  , OR
   BELOW STRIKE FOR FLOOR.  THEN NUM OF DAYS LIKELY TO BE INSIDE BAND =
   TOTAL DAYS - NUM OF DAYS LIKELY TO BE BELOW BAND - NUMBER OF DAYS LIKELY TO
   BE ABOVE BAND = TOTAL DAYS -
     @ACCRUAL("CAP"  ,CAPSTRIKE)-@ACCRUAL("FLOOR"  ,FLOORSTRIKE)
*/

/************************************************************************
INPUT KEY:

  Date today
  Date dcnt_start_date--first day for daycount purposes
  Date start_date--first accrual day
  Date end_date--last accrual day
  char *cap_floor_code--"CAP" or "FLOOR"
  char *win_lose--"WIN" or "LOSE"
  int lag--number of days to spot
  int end_lag--number of days at the end which are lumped together into one
    fixing
  double strike--level for determining if you accrue or not
  Date float_set--if you are receiving the floating rate  , when it is fixed
  double correlation--correlation between floating rate you receive and libor
  double float_tenor--length in months of floating rate you receive (ifyoudo)
  Date *bus_day_array--array of all business days
  double *rate_array--fixings and forward rates corresponding to bus_day_array
  int ra_len--length of prev two arrays
  double *vol_date_array--maturities for which we know vols
  double *vol_array--1d array of vols corresp. to vol_date_array(maybe NULL)
  double *strikes_vec--strikes for which we know vols(maybe NULL)
  double **mkt_vol_2darray--matrix of vols across mat and strike(maybe NULL)
  long m_mat_dates--num of vol dates
  long m_strikes--num of vol strikes
  int float_flg--1 if a floating rate is accruing  , 0 if a fixed rate.
    Note that if a floating rate is accruing  , we make a change of
    numeraire to that in which that rate is a martingale  , that is why
    (1)we don't need to know the forward rate for that rate  , and
    (2)we add adj to d2 below.
  int smile_flg--1 if we use smile effect of vol  , 0 if not.  If 1  , this
means we are going to look at mkt_vol_2darray; if 0 we are going to look at
    vol_array.
  int FHLB_flg--1 this routine knows two rules for business days.
    1:{Friday  , Saturday   ,Sunday} are all considered business days accruing
    off of Friday's Rate
    0:we don't count Saturday and Sunday.
  double *answer--answer = expected fraction of accruing days/business days
************************************************************************/

#define NUMD(DAY2, DAY1) (FHLB_flg ? (DAY2 - DAY1) : 1.0)

SrtErr srt_f_rngfltr(Date today, Date dcnt_start_date, Date start_date,
                     Date end_date, char *cap_floor_code, char *win_lose,
                     int lag, int end_lag, double strike, double float_fwd,
                     Date float_set, double correlation, double float_tenor,
                     Date *bus_day_array, double *rate_array, int ra_len,
                     double *vol_date_array, double *vol_array,
                     double *strikes_vec, double **mkt_vol_2darray,
                     long m_mat_dates, long m_strikes, int float_flg,
                     int smile_flg, int FHLB_flg, int reset_flg,
                     double *answer) {
  double cf, wl, ka, kao, optval, fraction_accruing, d2, vol, fvol;
  int method = 0; /* LINEAR INTERPOLATION */
  int i;
  int total_bus_days;

  /*** INPUT PROCESSING AND TESTING ***/

  strupper(cap_floor_code);
  if ((cap_floor_code[0] != 'C') && (cap_floor_code[0] != 'F'))
    return serror("Option must be a cap or floor");
  strupper(win_lose);
  if ((win_lose[0] != 'W') && (win_lose[0] != 'L'))
    return serror("The option owner must win or lose on the boundary");

  /**	CHECK 	DATES IN INCREASING ORDER
          THIS SECTION COMMENTED OUT FOR SPEED
    for(i=0;i<ra_len-1;i++)
    {
      if(bus_day_array[i] >= bus_day_array[i+1])
        return serror("@interp: bus days for rates must be increasing.");
    }
  ***/

  /**	CHECK 	DATES IN INCREASING ORDER
          THIS SECTION COMMENTED OUT FOR SPEED
    for(i=0;i<va_len-1;i++)
    {
      if(vol_date_array[i] >= vol_date_array[i+1])
        return serror("@interp: vol dates must be increasing.");
    }
  ***/

  /* NOW VALUE THE OPTIONS */
  /* COMPUTE THE APPROPRIATE STRIKE ADJUSTMENTS */
  cf = 1.0;
  if (cap_floor_code[0] == 'C')
    cf = (-1.0);
  wl = 1.0;
  if (win_lose[0] == 'L')
    wl = (-1.0);
  kao = strike + cf * wl * (1.0e-9);

  /* WE MUST ROUND KA TO THE NEAREST 32ND IN THE RIGHT DIRECTION */
  if (!reset_flg)
    ka = nearest_32nd(strike, cf, wl);
  else
    ka = strike;

  /* FIRST FIND THE INITIAL DATE TO VALUE OPTIONS FOR */
  i = 0;
  optval = 0.0;
  /* BUS_DAY_ARRAY IS THE RATE SET DATE */
  while (bus_day_array[i + lag] < dcnt_start_date)
    i++;
  total_bus_days = i;
  while (bus_day_array[i + lag] < start_date)
    i++;

  /* IF THIS IS A RESETTABLE THEN WE STIPULATE THAT THE BAND MUST BE SET
     BEFORE WE START TO ACCRUE */
  if (reset_flg && bus_day_array[i] < float_set) {
    return serror("band must be reset before start of accrual");
  }

  /* NOW DEAL WITH DAYS ALREADY SET BUT BEFORE LAST DAY */
  /* THE START AND END DATES ARE SPECIAL BECAUSE THEY MAY BE ON WEEKENDS
     THANKS TO THE MEAN OLD FHLB */
  if ((bus_day_array[i] <= today) && (bus_day_array[i + lag] > start_date)) {
    if (cap_floor_code[0] == 'C') {
      if (rate_array[i - 1] > kao)
        optval += 1.0 * NUMD(bus_day_array[i + lag], start_date);
    } else {
      if (rate_array[i - 1] < kao)
        optval += 1.0 * NUMD(bus_day_array[i + lag], start_date);
    }
  }
  while ((bus_day_array[i] <= today) &&
         (bus_day_array[i + lag + end_lag] < end_date)) {
    if (cap_floor_code[0] == 'C') {
      if (rate_array[i] > kao)
        optval +=
            1.0 * NUMD(bus_day_array[i + 1 + lag], bus_day_array[i + lag]);
    } else {
      if (rate_array[i] < kao)
        optval +=
            1.0 * NUMD(bus_day_array[i + 1 + lag], bus_day_array[i + lag]);
    }
    i++;
  }
  /* NOW DEAL WITH DAYS ALREADY SET AND ON LAST DAY */
  /*  WE NEED >= END_DATE TO ACCOUNT FOR FHLB WEEKEND RULES */
  if ((bus_day_array[i] <= today) &&
      (bus_day_array[i + lag + end_lag] >= end_date)) {
    if (cap_floor_code[0] == 'C') {
      if (rate_array[i] > kao)
        optval += 1.0 * NUMD(end_date, bus_day_array[i + lag]);
    } else {
      if (rate_array[i] < kao)
        optval += 1.0 * NUMD(end_date, bus_day_array[i + lag]);
    }
    /*    i++;*/
  }

  fvol = 0.0;
  if (float_set > today) {
    fvol = get_vol((double)float_set, (reset_flg ? float_fwd : strike),
                   vol_date_array, vol_array, strikes_vec, mkt_vol_2darray,
                   m_mat_dates, m_strikes, method, smile_flg);
  }
  /* NOW DEAL WITH DAYS TO VALUE OPTION BEFORE LAST DAY*/
  /* THE START AND END DATES ARE SPECIAL BECAUSE THEY MAY BE ON WEEKENDS
     THANKS TO THE MEAN OLD FHLB */
  if ((bus_day_array[i + lag] > start_date) &&
      (bus_day_array[i + lag - 1] < start_date) && (bus_day_array[i] > today)) {
    /* THIS IS THE BUG FIX */
    if (bus_day_array[i - 1] == today) {
      if (cap_floor_code[0] == 'C') {
        if (rate_array[i - 1] > kao)
          optval += 1.0 * NUMD(bus_day_array[i + lag], start_date);
      } else {
        if (rate_array[i - 1] < kao)
          optval += 1.0 * NUMD(bus_day_array[i + lag], start_date);
      }
    } else {
      vol =
          get_vol((double)bus_day_array[i - 1], strike + reset_flg * float_fwd,
                  vol_date_array, vol_array, strikes_vec, mkt_vol_2darray,
                  m_mat_dates, m_strikes, method, smile_flg);
      d2 = get_d2(rate_array[i - 1], float_fwd, ka, fvol, vol, float_set,
                  bus_day_array[i - 1], today, correlation, float_tenor,
                  float_flg, reset_flg);

      if (cap_floor_code[0] == 'C') {
        optval += norm(d2) * NUMD(bus_day_array[i + lag], start_date);
      } else {
        optval += (1.0 - norm(d2)) * NUMD(bus_day_array[i + lag], start_date);
      }
    }
    /* END BUG FIX */
  }
  while (bus_day_array[i + lag + end_lag] < end_date) {
    vol = get_vol((double)bus_day_array[i], strike + reset_flg * float_fwd,
                  vol_date_array, vol_array, strikes_vec, mkt_vol_2darray,
                  m_mat_dates, m_strikes, method, smile_flg);

    d2 = get_d2(rate_array[i], float_fwd, ka, fvol, vol, float_set,
                bus_day_array[i], today, correlation, float_tenor, float_flg,
                reset_flg);

    if (cap_floor_code[0] == 'C') {
      optval +=
          norm(d2) * NUMD(bus_day_array[i + 1 + lag], bus_day_array[i + lag]);
    } else {
      optval += (1.0 - norm(d2)) *
                NUMD(bus_day_array[i + 1 + lag], bus_day_array[i + lag]);
    }
    i++;
  }
  /* NOW DEAL WITH DAYS TO VALUE OPTION ON LAST DAY*/
  /*  WE NEED >= END_DATE TO ACCOUNT FOR FHLB WEEKEND RULES */

  if (FHLB_flg && (bus_day_array[i] > today) &&
      (bus_day_array[i + lag + end_lag] >= end_date)) {
    vol = get_vol((double)bus_day_array[i], strike + reset_flg * float_fwd,
                  vol_date_array, vol_array, strikes_vec, mkt_vol_2darray,
                  m_mat_dates, m_strikes, method, smile_flg);

    d2 = get_d2(rate_array[i], float_fwd, ka, fvol, vol, float_set,
                bus_day_array[i], today, correlation, float_tenor, float_flg,
                reset_flg);

    if (cap_floor_code[0] == 'C') {
      optval += norm(d2) * NUMD(end_date, bus_day_array[i + lag]);
    } else {
      optval += (1.0 - norm(d2)) * NUMD(end_date, bus_day_array[i + lag]);
    }
    i++;
  }

  fraction_accruing = optval;
  total_bus_days = i - total_bus_days;

  if (FHLB_flg) {
    fraction_accruing /= (double)(end_date - dcnt_start_date);
  } else {
    fraction_accruing /= (double)total_bus_days;
  }

  *answer = fraction_accruing;
  return (NULL);
}

/*
  RETURNS EXPECTED COUPON LOST FROM BEING OUTSIDE OF BAND
  ASSUMING OWNER CAN ADJUST POSITION OF BAND TO MINIMIZE THIS
  NUMBER.  WE ASSUME THERE IS AN OPTIMAL STRATEGY OF THE FORM:
  "WHEN IT IS TIME TO RESET THE BAND I AM GOING TO CENTER THE BAND AROUND
  FWD LIBOR + K" (K IS UNKNOWN); WE ALSO ASSUME THAT THE NO. DAYS
  NOT ACCRUING IS A QUADRATIC
  FUNCTION IN K  , WHOSE ARGMIN CAN THEREFORE BE COMPUTED AFTER TRYING
  THREE DIFFERENT VALUES FOR K.
*/
SrtErr srt_f_rngfltr_resettable(
    Date today, Date dcnt_start_date, Date start_date, Date end_date,
    char *win_lose_cap, char *win_lose_floor, int lag, int end_lag,
    double strike, double band_width, double float_fwd_fwd, double float_fwd,
    double spr_or_cpn, Date float_set, double correlation, double float_tenor,
    Date *bus_day_array, double *rate_array, int ra_len, double *vol_date_array,
    double *vol_array, double *strikes_vec, double **mkt_vol_2darray,
    long m_mat_dates, long m_strikes, int float_flg, int smile_flg,
    int FHLB_flg, double *answer) {
  double val, val2;
  SrtErr err;
  double x[3];
  double y[3];
  int i, j;

  /* check if band has already been set  , in which case we take strike
  to be the bottom of the band */
  if (today > float_set) {
    val2 = 0;
    for (j = 0; j < float_flg + 1; j++) {

      err = srt_f_rngfltr(today, dcnt_start_date, start_date, end_date, "CAP",
                          win_lose_cap, lag, end_lag, strike + band_width,
                          float_fwd, float_set, correlation, float_tenor,
                          bus_day_array, rate_array, ra_len, vol_date_array,
                          vol_array, strikes_vec, mkt_vol_2darray, m_mat_dates,
                          m_strikes, j, smile_flg, FHLB_flg, 0, &val);
      if (err)
        return err;
      val2 += val * (j ? float_fwd : spr_or_cpn);

      err = srt_f_rngfltr(today, dcnt_start_date, start_date, end_date, "FLOOR",
                          win_lose_floor, lag, end_lag, strike, float_fwd,
                          float_set, correlation, float_tenor, bus_day_array,
                          rate_array, ra_len, vol_date_array, vol_array,
                          strikes_vec, mkt_vol_2darray, m_mat_dates, m_strikes,
                          j, smile_flg, FHLB_flg, 0, &val);
      if (err)
        return err;
      val2 += val * (j ? float_fwd : spr_or_cpn);
    }
    *answer = val2;
    return NULL;
  } else {
    x[0] = -band_width * 0.5;
    x[1] = 0.5 * (float_fwd_fwd - float_fwd) + x[0];
    x[2] = float_fwd_fwd - float_fwd + x[0];
    memset(y, 0, 3 * sizeof(double));

    for (i = 0; i < 3; i++) {
      for (j = 0; j < float_flg + 1; j++) {
        err = srt_f_rngfltr(
            today, dcnt_start_date, start_date, end_date, "CAP", win_lose_cap,
            lag, end_lag, x[i] + band_width, float_fwd, float_set, correlation,
            float_tenor, bus_day_array, rate_array, ra_len, vol_date_array,
            vol_array, strikes_vec, mkt_vol_2darray, m_mat_dates, m_strikes, j,
            smile_flg, FHLB_flg, 1, &val);
        y[i] += val * (j ? float_fwd : spr_or_cpn);

        err = srt_f_rngfltr(
            today, dcnt_start_date, start_date, end_date, "FLOOR",
            win_lose_floor, lag, end_lag, x[i], float_fwd, float_set,
            correlation, float_tenor, bus_day_array, rate_array, ra_len,
            vol_date_array, vol_array, strikes_vec, mkt_vol_2darray,
            m_mat_dates, m_strikes, j, smile_flg, FHLB_flg, 1, &val);
        y[i] += val * (j ? float_fwd : spr_or_cpn);
      }
    }
    /*NOTE THAT WE WANT TO MAXIMIZE THE NUMBER OF DAYS INSIDE THE BAND--->
    MINIMIZE THE NUMBER OF DAYS OUTSIDE THE BAND  , WHICH IS WHAT srt_f_rngfltr
    COMPUTES.*/

    strike = quadratic_minmax(x, y);
    val2 = 0;
    for (j = 0; j < float_flg + 1; j++) {
      err = srt_f_rngfltr(today, dcnt_start_date, start_date, end_date, "CAP",
                          win_lose_cap, lag, end_lag, strike + band_width,
                          float_fwd, float_set, correlation, float_tenor,
                          bus_day_array, rate_array, ra_len, vol_date_array,
                          vol_array, strikes_vec, mkt_vol_2darray, m_mat_dates,
                          m_strikes, j, smile_flg, FHLB_flg, 1, &val);
      val2 += val * (j ? float_fwd : spr_or_cpn);

      err = srt_f_rngfltr(today, dcnt_start_date, start_date, end_date, "FLOOR",
                          win_lose_floor, lag, end_lag, strike, float_fwd,
                          float_set, correlation, float_tenor, bus_day_array,
                          rate_array, ra_len, vol_date_array, vol_array,
                          strikes_vec, mkt_vol_2darray, m_mat_dates, m_strikes,
                          j, smile_flg, FHLB_flg, 1, &val);
      val2 += val * (j ? float_fwd : spr_or_cpn);
    }
    for (i = 0; i < 3; i++) {
      val2 = DMIN(val2, y[i]);
    }
    *answer = val2;
    return NULL;
  }
}
