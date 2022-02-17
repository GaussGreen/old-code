/* =======================================================================================

   FUNCNAME        :swp_f_cms

   DESCRIPTION     :computes value of cms through integration with swaptions

   =======================================================================================
 */

#include "math.h"
#include "opfnctns.h"
#include "swp_h_all.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"

#define JMAX 30
/* ----------------------------------------------------------------------------------
 */

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :swp_f_cmswap
  AUTHOR          :E.Auld from code by P.Varnish
  DESCRIPTION     :computes value of cms option        ,
                                        computing the forward rate itself
  MODIFIES        :
  CALL            :

<%%END>---------------------------------------------------------------------*/

Err swp_f_cmswp(SwapDP *swapdp, String refRatecode, double strike,
                double volatility, SrtReceiverType pay_or_receive,
                Date pay_date, double delta, int nswaps, String ycName,
                SrtDiffusionType lognrm_nrm, double *ans) {
  double forward, maturity, cms, delay;
  Err err;
  Date fix_date;
  SrtCurvePtr yldcrv;

  err = swp_f_ForwardRate_SwapDP(swapdp, ycName, refRatecode, &forward);
  if (err)
    return err;

  yldcrv = lookup_curve(ycName);
  if (!yldcrv)
    return serror("Could not find %s yield curve", ycName);

  fix_date = add_unit(swapdp->start, -swapdp->spot_lag, SRT_BDAY, SUCCEEDING);
  maturity =
      coverage(get_clcndate_from_yldcrv(yldcrv), fix_date, BASIS_ACT_365);

  if (swapdp->basis_code == BASIS_ACT_360) {
    forward = forward * 365.0 / 360.0;
    strike = strike * 365.0 / 360.0;
  }
  if (pay_date < swapdp->start)
    return serror("pay date < swap start cms");

  delay = coverage(swapdp->start, pay_date, BASIS_ACT_365);

  /* so delay is number of years */

  err = swp_f_cmsoption(forward, swapdp->nfp, swapdp->compd, strike, volatility,
                        maturity, pay_or_receive, delay, delta, nswaps,
                        lognrm_nrm, &cms);
  if (err)
    return err;

  if (swapdp->basis_code == BASIS_ACT_360) {
    cms *= 360.0 / 365.0;
  }

  (*ans) = cms;

  return NULL;
}

/* -------------------------------------------------------------------------------
 */
/* Returns the value of a CMS rate        , as Call(F) - Put(F) + F */
Err swp_f_cmsrate(double forward, double num_periods, SrtCompounding frequency,
                  double volatility, double maturity, double delay,
                  double delta, int num_swaps, SrtDiffusionType lognrm_nrm,
                  double *ans) {
  Err err = NULL;
  double cms_receiver;
  double cms_payer;

  /* Computes the value of a Receiver CMS Option        , strike at Forward */
  err = swp_f_cmsoption(forward, num_periods, frequency, forward, volatility,
                        maturity, SRT_RECEIVER, delay, delta, num_swaps,
                        lognrm_nrm, &cms_receiver);

  /* Computes the value of a Payer CMS Option        , strike at Forward */
  err = swp_f_cmsoption(forward, num_periods, frequency, forward, volatility,
                        maturity, SRT_PAYER, delay, delta, num_swaps,
                        lognrm_nrm, &cms_payer);

  /* The Final value: Payer - Receiver + Forward */
  *ans = cms_payer - cms_receiver + forward;

  /* Return a success message */
  return NULL;
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :swp_f_cmswap
  DESCRIPTION     :Returns the value of a CMS rate        , as Call(F) - Put(F)
+ F computing the forward rate itself MODIFIES        : CALL            :
<%%END>---------------------------------------------------------------------*/

Err swp_f_cmsratefwd(SwapDP *swapdp, String refRatecode, double volatility,
                     Date pay_date, double delta, int nswaps, String ycName,
                     SrtDiffusionType lognrm_nrm, double *ans) {
  double forward;
  Err err = NULL;
  Date today;
  SrtCurvePtr yldcrv;

  err = swp_f_ForwardRate_SwapDP(swapdp, ycName, refRatecode, &forward);
  if (err)
    return err;

  yldcrv = lookup_curve(ycName);
  if (!yldcrv)
    return serror("Could not find %s yield curve", ycName);

  today = get_clcndate_from_yldcrv(yldcrv);

  err = swp_f_cmsratedp(swapdp, today, forward, volatility, pay_date, delta,
                        nswaps, lognrm_nrm, ans);

  return err;
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :swp_f_cmsratedp
  DESCRIPTION     :Returns the value of a CMS rate        , as Call(F) - Put(F)
+ F from a swapdp MODIFIES        : CALL            :
<%%END>---------------------------------------------------------------------*/

Err swp_f_cmsratedp(SwapDP *swapdp, Date today, double forward,
                    double volatility, Date pay_date, double delta, int nswaps,
                    SrtDiffusionType lognrm_nrm, double *ans) {
  Date fix_date;
  double maturity, delay;
  double cms;
  double nfp;
  Err err;

  fix_date = add_unit(swapdp->start, -swapdp->spot_lag, SRT_BDAY, SUCCEEDING);
  maturity = coverage(today, fix_date, BASIS_ACT_365);

  if (swapdp->basis_code == BASIS_ACT_360) {
    forward = forward * 365.0 / 360.0;
  }

  if (pay_date < swapdp->start)
    return serror("pay date < swap start cms");

  delay = coverage(swapdp->start, pay_date, BASIS_ACT_365);

  /*	Calculate number of full periods	*/
  if (swapdp->nfp < 400) {
    nfp = (double)swapdp->nfp;
  } else {
    nfp = (double)(swapdp->end - swapdp->start);
    nfp /= 365;
    nfp *= swapdp->compd;
  }

  /* so delay is number of years */
  err = swp_f_cmsrate(forward, nfp, swapdp->compd, volatility, maturity, delay,
                      delta, nswaps, lognrm_nrm, &cms);
  if (err)
    return err;

  if (swapdp->basis_code == BASIS_ACT_360) {
    cms *= 360.0 / 365.0;
  }

  (*ans) = cms;

  return NULL;
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :swp_f_cmsratedp_clsdfrm
  DESCRIPTION     :Returns the value of a CMS rate from a closed-form
approximation MODIFIES        : CALL            :
<%%END>---------------------------------------------------------------------*/

Err swp_f_cmsratedp_clsdfrm(SwapDP *swapdp, Date today, double forward,
                            double volatility, Date pay_date,
                            SrtDiffusionType lognrm_nrm, double *ans) {
  Date fix_date;
  double maturity, delay;
  double cms;
  double nfp;
  double nvol;
  double lvl, dlvl;
  double cms_adj, del_adj, tot_adj;

  fix_date = add_unit(swapdp->start, -swapdp->spot_lag, SRT_BDAY, SUCCEEDING);
  if (fix_date < today)
    return serror("fix_date < today cms");
  maturity = (fix_date - today) * YEARS_IN_DAY;

  if (pay_date < fix_date)
    return serror("pay date < swap start cms");
  delay = (pay_date - fix_date) * YEARS_IN_DAY;

  if (swapdp->basis_code == BASIS_ACT_360) {
    forward = forward * 365.0 / 360.0;
  }

  /*	Calculate number of full periods	*/
  if (swapdp->nfp < 400) {
    nfp = (double)swapdp->nfp;
  } else {
    nfp = (double)(swapdp->end - swapdp->start);
    nfp /= 365;
    nfp *= swapdp->compd;
  }

  /*	LVL and LVL'	*/
  lvl = (1.0 - pow(1.0 + forward / swapdp->compd, -nfp)) /
        (forward / swapdp->compd);
  dlvl = (nfp / swapdp->compd) * pow(1.0 + forward / swapdp->compd, -nfp - 1) *
             forward / swapdp->compd -
         (1.0 - pow(1.0 + forward / swapdp->compd, -nfp)) / swapdp->compd;
  dlvl /= (forward / swapdp->compd) * (forward / swapdp->compd);

  /*	Adjustments	*/
  cms_adj = -dlvl / lvl;
  del_adj = -delay / (1.0 + forward * delay);
  tot_adj = cms_adj + del_adj;

  switch (lognrm_nrm) {

  case SRT_LOGNORMAL:

    nvol = forward * volatility;
    break;

  case SRT_NORMAL:
    nvol = volatility;

    break;

  default:
    return serror("Vol type unknown");
    break;
  }

  if (swapdp->basis_code == BASIS_ACT_360) {
    forward *= 360.0 / 365.0;
  }

  cms = forward + tot_adj * nvol * nvol * maturity;

  (*ans) = cms;

  return NULL;
}

Err swp_f_Cms_Rate(
    double dFwdSwapRate,                   /* Forward Swap Rate */
    double dMaturity,                      /* Cms Maturity as double */
    double dNumPeriods, double dFrequency, /* Cms Frequency */
    double dDelay, double dRateConv,       /* date adjustment */
    SrtDiffusionType VolType,              /* Vol type Lognormal or Normal */
    double dFlatVol,                       /* Flat Vol if used */
    int iMethod, /* 0: Use Flat Vol        , 1: linear interpolation        , 2:
                    FullSmile*/
    Date dStart, /* The following parameters are used in the GetVol function */
    Date dEnd,   /* and are useless in the rest of the code */
    SRT_Boolean bAdjForSpread, double dSpread, char *szVolCurveName,
    long lNumStrikesInVol, double *pdStrikesVol, double *dCmsRateValue) {
  double dCmsPayVal = 0.0, dCmsRecVal = 0.0;
  Err err;

  /* Get value of Put (i.e. receiver's) option with strike = forward */
  err = swp_f_Cms_Option(dFwdSwapRate, dMaturity, dNumPeriods, dFwdSwapRate,
                         dFrequency, SRT_RECEIVER, dDelay, dRateConv, VolType,
                         dFlatVol, iMethod, dStart, dEnd, bAdjForSpread,
                         dSpread, szVolCurveName, lNumStrikesInVol,
                         pdStrikesVol, &dCmsRecVal);

  /* Get value of Call (i.e. payer's) option with strike = forward */
  err = swp_f_Cms_Option(dFwdSwapRate, dMaturity, dNumPeriods, dFwdSwapRate,
                         dFrequency, SRT_PAYER, dDelay, dRateConv, VolType,
                         dFlatVol, iMethod, dStart, dEnd, bAdjForSpread,
                         dSpread, szVolCurveName, lNumStrikesInVol,
                         pdStrikesVol, &dCmsPayVal);

  /* Use call-put parity: conv(forward) - forward = call - put */
  *dCmsRateValue = dCmsPayVal - dCmsRecVal + dFwdSwapRate;
  return err;
}

/* -------------------------------------------------------------------------------
 */
/* Returns the value of a CMS rate        , as Call(F) - Put(F) + F
 */
Err swp_f_Tec_Rate(
    double dFwdSwapRate, /*Forward Swap Rate */
    double dMaturity, double dNumPeriods, double dAtmStrike, double dMargin,
    double dFrequency,               /* Tec Frequency */
    double dPaymentPower,            /* Tec Power */
    double dDelay, double dRateConv, /* date adjustment */
    SrtDiffusionType VolType, double dFlatVol,
    int iMethod, /* 0: Use Flat Vol        , 1: linear interpolation        , 2:
                    FullSmile*/
    Date dStart, /* The following parameters are used in the GetVol function */
    Date dEnd,   /* and are useless in the rest of the code */
    SRT_Boolean bAdjForSpread, double dSpread, char *szVolCurveName,
    long lNumStrikesInVol, double *pdStrikesVol, double *dTecRateValue) {
  double dTecPayVal = 0.0, dTecRecVal = 0.0, dPaymentEquivalent = 0.0;
  Err err;

  dAtmStrike += dMargin;

  /* Get value of Put (i.e. receiver's) option with strike = ATMStrike */
  err = swp_f_Tec_Option(dFwdSwapRate, dMaturity, dNumPeriods, dAtmStrike,
                         dMargin, dFrequency, dPaymentPower, SRT_RECEIVER,
                         dDelay, dRateConv, VolType, dFlatVol, iMethod, dStart,
                         dEnd, bAdjForSpread, dSpread, szVolCurveName,
                         lNumStrikesInVol, pdStrikesVol, &dTecRecVal);

  /* Get value of Call (i.e. payer's) option with strike = ATMStrike */
  err = swp_f_Tec_Option(dFwdSwapRate, dMaturity, dNumPeriods, dAtmStrike,
                         dMargin, dFrequency, dPaymentPower, SRT_PAYER, dDelay,
                         dRateConv, VolType, dFlatVol, iMethod, dStart, dEnd,
                         bAdjForSpread, dSpread, szVolCurveName,
                         lNumStrikesInVol, pdStrikesVol, &dTecPayVal);

  /* Use call-put parity: conv(forward) - forward = call - put */
  dPaymentEquivalent = dTecPayVal - dTecRecVal -
                       (1 - pow((1 + dAtmStrike), (1 / dPaymentPower)));

  /* Gets the Equivalent Forward (Payment = ( 1+ TEC +m) ^1/4 - 1) */
  *dTecRateValue =
      pow(dPaymentEquivalent + 1.0, dPaymentPower) - (1.0 + dMargin);

  return err;
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :swp_f_cmsrateNewYork
  DESCRIPTION     :Returns the value of a CMS rate using the SABR model through
a one dimensional numerical integration AUTHOR	      :Ezra Nahum LAST MODIFIED
:October 2000
<%%END>---------------------------------------------------------------------*/
double cms_sm_trapzd(double (*function)(double, double, double, double, double,
                                        double, double, double, double),
                     double forward, double mat, double betavol, double alpha,
                     double beta, double rho, double numperiod, double coverage,
                     double a, double b, int n, double accum) {
  double x, tnm, sum, del, ret;
  long it, j;

  if (n == 1) {
    accum = 0.5 * (b - a) *
            (function(a, forward, mat, betavol, alpha, beta, rho, numperiod,
                      coverage) +
             function(b, forward, mat, betavol, alpha, beta, rho, numperiod,
                      coverage));
  } else {
    it = 1;
    for (j = 1; j < n - 1; j++)
      it *= 2;
    tnm = (double)it;
    del = (b - a) / tnm;
    x = a + 0.5 * del;
    sum = 0.0;
    for (j = 1; j <= it; j++) {
      sum += function(x, forward, mat, betavol, alpha, beta, rho, numperiod,
                      coverage);
      x += del;
    }
    accum = 0.5 * (accum + (b - a) * sum / tnm);
  }
  ret = accum;
  return ret;
} /* END sm_trapzd */

double cms_sm_qsimp(double (*func)(double, double, double, double, double,
                                   double, double, double, double),
                    double forward, double mat, double betavol, double alpha,
                    double beta, double rho, double numperiod, double coverage,
                    double a, double b, double precision) {
  int j;
  double sum;
  double oldsum;
  double sum2N;
  double sumN;

  sumN = cms_sm_trapzd(func, forward, mat, betavol, alpha, beta, rho, numperiod,
                       coverage, a, b, 1, 0.0);
  sum2N = oldsum = -1.0e30;

  for (j = 2; j < JMAX; j++) {
    sum2N = cms_sm_trapzd(func, forward, mat, betavol, alpha, beta, rho,
                          numperiod, coverage, a, b, j, sumN);

    sum = ((4.0 * sum2N) - sumN) / 3.0;
    if (fabs(sum - oldsum) <= (precision * fabs(oldsum)))
      return sum;
    oldsum = sum;
    sumN = sum2N;
  }

  return sum;

} /* END sm_qsimp(...) */

double cms_sm_trapzd2(double (*function)(double, double, double, double, double,
                                         double, double, double, double, double,
                                         double, double),
                      double forward, double mat, double betavol, double alpha,
                      double beta, double rho, double numperiod,
                      double coverage, double K1, double K2, double q, double a,
                      double b, int n, double accum) {
  double x, tnm, sum, del, ret;
  long it, j;

  if (n == 1) {
    accum = 0.5 * (b - a) *
            (function(a, forward, mat, betavol, alpha, beta, rho, numperiod,
                      coverage, K1, K2, q) +
             function(b, forward, mat, betavol, alpha, beta, rho, numperiod,
                      coverage, K1, K2, q));
  } else {
    it = 1;
    for (j = 1; j < n - 1; j++)
      it *= 2;
    tnm = (double)it;
    del = (b - a) / tnm;
    x = a + 0.5 * del;
    sum = 0.0;
    for (j = 1; j <= it; j++) {
      sum += function(x, forward, mat, betavol, alpha, beta, rho, numperiod,
                      coverage, K1, K2, q);
      x += del;
    }
    accum = 0.5 * (accum + (b - a) * sum / tnm);
  }
  ret = accum;
  return ret;
} /* END sm_trapzd */

double cms_sm_qsimp2(double (*func)(double, double, double, double, double,
                                    double, double, double, double, double,
                                    double, double),
                     double forward, double mat, double betavol, double alpha,
                     double beta, double rho, double numperiod, double coverage,
                     double K1, double K2, double q, double a, double b,
                     double precision) {
  int j;
  double sum;
  double oldsum;
  double sum2N;
  double sumN;

  sumN = cms_sm_trapzd2(func, forward, mat, betavol, alpha, beta, rho,
                        numperiod, coverage, K1, K2, q, a, b, 1, 0.0);
  sum2N = oldsum = -1.0e30;

  for (j = 2; j < JMAX; j++) {
    sum2N = cms_sm_trapzd2(func, forward, mat, betavol, alpha, beta, rho,
                           numperiod, coverage, K1, K2, q, a, b, j, sumN);

    sum = ((4.0 * sum2N) - sumN) / 3.0;
    if (fabs(sum - oldsum) <= (precision * fabs(oldsum)))
      return sum;
    oldsum = sum;
    sumN = sum2N;
  }

  return sum;

} /* END sm_qsimp(...) */

/*
double ProbaIntegrand(double s        , double S0        , double mat ,double
sigb        , double a        , double b        , double r        , double m ,
double d)
{
double SIGMA;
double RV;
double integ;
double d1        ,d2;
double d2bsd2k        ,d2bsdsigdk        ,dbsdsig        ,d2bsd2sig;
double dsigdk        ,d2sigd2k;
double h;
double sph        ,smh;

h = 0.00001;

RV = s*s/(1-pow(1+d*s        ,-m));


srt_f_optbetastochvoltoblkvol(S0        ,
                                                           s        ,
                                                        sigb        ,
                                                                a        ,
                                                                r        ,
                                                        mat        ,
                                                                b        ,
                                                        &SIGMA);


d1 = log(S0/s)/(SIGMA*sqrt(mat))+SIGMA*sqrt(mat)/2;
d2 = d1-SIGMA*sqrt(mat);

dbsdsig =
S0*(-log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*exp(-d1*d1/2)/sqrt(2*SRT_PI)+
                  s*(log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*exp(-d2*d2/2)/sqrt(2*SRT_PI);

d2bsd2sig =
(2*S0*log(S0/s)/(SIGMA*SIGMA*SIGMA*sqrt(mat))-S0*(-log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*
                         (-log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*d1)*exp(-d1*d1/2)/sqrt(2*SRT_PI)-
                         (2*s*log(S0/s)/(SIGMA*SIGMA*SIGMA*sqrt(mat))+s*(log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*
                         (log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*d2)*exp(-d2*d2/2)/sqrt(2*SRT_PI);


d2bsd2k =
S0*(0.5-log(S0/s)/(SIGMA*SIGMA*mat))*exp(-d1*d1/2)/(sqrt(2*SRT_PI)*SIGMA*sqrt(mat)*s*s)+
                  (0.5+log(S0/s)/(SIGMA*SIGMA*mat))*exp(-d2*d2/2)/(sqrt(2*SRT_PI)*SIGMA*sqrt(mat)*s);


d2bsdsigdk  =
(S0/(s*SIGMA*SIGMA*sqrt(mat))+S0*(-log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*d1/(SIGMA*s*sqrt(mat)))*exp(-d1*d1/2)/sqrt(2*SRT_PI)
                          +(log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2-1/(SIGMA*SIGMA*sqrt(mat))+
                          (log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*d2/(SIGMA*sqrt(mat)))*exp(-d2*d2/2)/sqrt(2*SRT_PI);

srt_f_optbetastochvoltoblkvol(S0        ,
                                                           s+h        ,
                                                        sigb        ,
                                                                a        ,
                                                                r        ,
                                                        mat        ,
                                                                b        ,
                                                        &sph);

srt_f_optbetastochvoltoblkvol(S0        ,
                                                           s-h        ,
                                                        sigb        ,
                                                                a        ,
                                                                r        ,
                                                        mat        ,
                                                                b        ,
                                                        &smh);

dsigdk = (sph-smh)/(2*h);
d2sigd2k = (sph-2*SIGMA+smh)/(h*h);



integ = d2bsd2k+2*d2bsdsigdk*dsigdk+dbsdsig*d2sigd2k+d2bsd2sig*dsigdk*dsigdk;

return integ;
}


double CMSIntegrand(double s        , double S0        , double mat ,double sigb
, double a        , double b        , double r        , double m        , double
d)
{
double SIGMA;
double RV;
double integ;
double d1        ,d2;
double d2bsd2k        ,d2bsdsigdk        ,dbsdsig        ,d2bsd2sig;
double dsigdk        ,d2sigd2k;
double h;
double sph        ,smh;

h = 0.00001;

RV = s*s/(1-pow(1+d*s        ,-m));


srt_f_optbetastochvoltoblkvol(S0        ,
                                                           s        ,
                                                        sigb        ,
                                                                a        ,
                                                                r        ,
                                                        mat        ,
                                                                b        ,
                                                        &SIGMA);


d1 = log(S0/s)/(SIGMA*sqrt(mat))+SIGMA*sqrt(mat)/2;
d2 = d1-SIGMA*sqrt(mat);

dbsdsig =
S0*(-log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*exp(-d1*d1/2)/sqrt(2*SRT_PI)+
                  s*(log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*exp(-d2*d2/2)/sqrt(2*SRT_PI);

d2bsd2sig =
(2*S0*log(S0/s)/(SIGMA*SIGMA*SIGMA*sqrt(mat))-S0*(-log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*
                         (-log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*d1)*exp(-d1*d1/2)/sqrt(2*SRT_PI)-
                         (2*s*log(S0/s)/(SIGMA*SIGMA*SIGMA*sqrt(mat))+s*(log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*
                         (log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*d2)*exp(-d2*d2/2)/sqrt(2*SRT_PI);


d2bsd2k =
S0*(0.5-log(S0/s)/(SIGMA*SIGMA*mat))*exp(-d1*d1/2)/(sqrt(2*SRT_PI)*SIGMA*sqrt(mat)*s*s)+
                  (0.5+log(S0/s)/(SIGMA*SIGMA*mat))*exp(-d2*d2/2)/(sqrt(2*SRT_PI)*SIGMA*sqrt(mat)*s);


d2bsdsigdk  =
(S0/(s*SIGMA*SIGMA*sqrt(mat))+S0*(-log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*d1/(SIGMA*s*sqrt(mat)))*exp(-d1*d1/2)/sqrt(2*SRT_PI)
                          +(log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2-1/(SIGMA*SIGMA*sqrt(mat))+
                          (log(S0/s)/(SIGMA*SIGMA*sqrt(mat))+sqrt(mat)/2)*d2/(SIGMA*sqrt(mat)))*exp(-d2*d2/2)/sqrt(2*SRT_PI);

srt_f_optbetastochvoltoblkvol(S0        ,
                                                           s+h        ,
                                                        sigb        ,
                                                                a        ,
                                                                r        ,
                                                        mat        ,
                                                                b        ,
                                                        &sph);

srt_f_optbetastochvoltoblkvol(S0        ,
                                                           s-h        ,
                                                        sigb        ,
                                                                a        ,
                                                                r        ,
                                                        mat        ,
                                                                b        ,
                                                        &smh);

dsigdk = (sph-smh)/(2*h);
d2sigd2k = (sph-2*SIGMA+smh)/(h*h);



integ =
RV*(d2bsd2k+2*d2bsdsigdk*dsigdk+dbsdsig*d2sigd2k+d2bsd2sig*dsigdk*dsigdk);

return integ;
}
*/

double ProbaIntegrand(double s, double S0, double mat, double sigb, double a,
                      double b, double r, double m, double d) {
  double SIGMA;
  double Expo;
  double normalize;
  double RV;
  double Csi;
  double XofZ;
  double integ;
  double omega;

  omega = sqrt(pow(s, 3 * b) / pow(S0, b));

  RV = s * s / (1 - pow(1 + d * s, -m));
  if (b == 1.0)
    Csi = a * (log(S0) - log(s)) / sigb;
  else
    Csi = a * (pow(S0, 1 - b) - pow(s, 1 - b)) / (sigb * (1 - b));

  XofZ = log((sqrt(1 - 2 * r * Csi + Csi * Csi) + Csi - r) / (1 - r));

  if (b == 1.0) {
    if (a == 0.0) {
      if ((S0 < s + 0.000000001) && (S0 > s - 0.000000001)) {
        SIGMA = s * sigb;
      } else {
        SIGMA = (S0 - s) * sigb / log(S0 / s);
      }
    } else {
      if ((S0 < s + 0.00001) && (S0 > s - 0.00001)) {
        SIGMA = s * sigb;
      } else {
        SIGMA = a * (S0 - s) / XofZ;
      }
    }

  }

  else if (a == 0) {
    if ((S0 < s + 0.00001) && (S0 > s - 0.00001)) {
      SIGMA = sigb * s;
    } else
      SIGMA = (sigb * (1 - b)) * (S0 - s) / (pow(S0, 1 - b) - pow(s, 1 - b));
  } else {
    if ((S0 < s + 0.00001) && (S0 > s - 0.00001)) {
      SIGMA = sigb * s;
    } else
      SIGMA = a * (S0 - s) / XofZ;
  }

  if ((b == 1.0) && (a == 0.0)) {
    normalize = 1 / (sqrt(2 * SRT_PI * mat * s * s * s * sigb * sigb / S0));
    Expo = exp(-(S0 - s) * (S0 - s) / (2 * mat * SIGMA * SIGMA) -
               mat * sigb * sigb / 8);
  } else {

    normalize = pow(Csi * Csi - 2 * r * Csi + 1, -0.75) /
                (sqrt(2 * SRT_PI * mat * sigb * sigb * omega * omega));
    Expo = exp(-(S0 - s) * (S0 - s) / (2 * mat * SIGMA * SIGMA));

    /*	Expo =
     * exp(-(S0-s)*(S0-s)/(2*mat*SIGMA*SIGMA))*(1-(((XofZ*(exp(XofZ)+exp(-XofZ))/(exp(XofZ)-exp(-XofZ)))-1)/(XofZ*XofZ)+1)*mat*a*a/8);
     */
  }

  integ = normalize * Expo;

  return integ;
}

double CMSIntegrand(double s, double S0, double mat, double sigb, double a,
                    double b, double r, double m, double d, double K1,
                    double K2, double q) {
  double SIGMA, SIGMA1, SIGMA2;
  double Expo;
  double normalize;
  double RV;
  double Csi, Csi1, Csi2;
  double XofZ, XofZ1, XofZ2;
  double integ;
  double omega, omega1, omega2;
  double p1, p2;

  omega = sqrt(pow(s, 3 * b) / pow(S0, b));
  omega1 = sqrt(pow(K1, 3 * b) / pow(S0, b));
  omega2 = sqrt(pow(K2, 3 * b) / pow(S0, b));

  RV = s * s / (1 - pow(1 + d * s, -m));

  if (b == 1.0) {
    Csi = a * (log(S0) - log(s)) / sigb;
    Csi1 = a * (log(S0) - log(K1)) / sigb;
    Csi2 = a * (log(S0) - log(K2)) / sigb;
  }

  else {
    Csi = a * (pow(S0, 1 - b) - pow(s, 1 - b)) / (sigb * (1 - b));
    Csi1 = a * (pow(S0, 1 - b) - pow(K1, 1 - b)) / (sigb * (1 - b));
    Csi2 = a * (pow(S0, 1 - b) - pow(K2, 1 - b)) / (sigb * (1 - b));
  }

  XofZ = log((sqrt(1 - 2 * r * Csi + Csi * Csi) + Csi - r) / (1 - r));
  XofZ1 = log((sqrt(1 - 2 * r * Csi1 + Csi1 * Csi1) + Csi1 - r) / (1 - r));
  XofZ2 = log((sqrt(1 - 2 * r * Csi2 + Csi2 * Csi2) + Csi2 - r) / (1 - r));

  if (b == 1.0) {
    if (a == 0.0) {
      if ((S0 < s + 0.000000001) && (S0 > s - 0.000000001)) {
        SIGMA = s * sigb;
        SIGMA1 = K1 * sigb;
        SIGMA2 = K2 * sigb;
      } else {
        SIGMA = (S0 - s) * sigb / log(S0 / s);
        SIGMA1 = (S0 - K1) * sigb / log(S0 / K1);
        SIGMA2 = (S0 - K2) * sigb / log(S0 / K2);
      }
    } else {
      if ((S0 < s + 0.00001) && (S0 > s - 0.00001)) {
        SIGMA = s * sigb;
        SIGMA1 = K1 * sigb;
        SIGMA2 = K2 * sigb;
      } else {
        SIGMA = a * (S0 - s) / XofZ;
        SIGMA1 = a * (S0 - K1) / XofZ1;
        SIGMA2 = a * (S0 - K2) / XofZ2;
      }
    }

  }

  else if (a == 0) {
    if ((S0 < s + 0.00001) && (S0 > s - 0.00001)) {
      SIGMA = sigb * s;
      SIGMA1 = sigb * K1;
      SIGMA2 = sigb * K2;
    } else {
      SIGMA = (sigb * (1 - b)) * (S0 - s) / (pow(S0, 1 - b) - pow(s, 1 - b));
      SIGMA1 = (sigb * (1 - b)) * (S0 - K1) / (pow(S0, 1 - b) - pow(K1, 1 - b));
      SIGMA2 = (sigb * (1 - b)) * (S0 - K2) / (pow(S0, 1 - b) - pow(K2, 1 - b));
    }
  } else {
    if ((S0 < s + 0.00001) && (S0 > s - 0.00001)) {
      SIGMA = sigb * s;
      SIGMA1 = sigb * K1;
      SIGMA2 = sigb * K2;
    } else {
      SIGMA = a * (S0 - s) / XofZ;
      SIGMA1 = a * (S0 - K1) / XofZ1;
      SIGMA2 = a * (S0 - K2) / XofZ2;
    }
  }

  if ((b == 1.0) && (a == 0.0)) {
    normalize = 1 / (sqrt(2 * SRT_PI * mat * s * s * s * sigb * sigb / S0));
    Expo = exp(-(S0 - s) * (S0 - s) / (2 * mat * SIGMA * SIGMA) -
               mat * sigb * sigb / 8);
    p1 = exp(-(S0 - K1) * (S0 - K1) / (2 * mat * SIGMA1 * SIGMA1) -
             mat * sigb * sigb / 8) *
         1 / (sqrt(2 * SRT_PI * mat * K1 * K1 * K1 * sigb * sigb / S0));
    p2 = exp(-(S0 - K2) * (S0 - K2) / (2 * mat * SIGMA2 * SIGMA2) -
             mat * sigb * sigb / 8) *
         1 / (sqrt(2 * SRT_PI * mat * K2 * K2 * K2 * sigb * sigb / S0));
  } else {

    normalize = pow(Csi * Csi - 2 * r * Csi + 1, -0.75) /
                (sqrt(2 * SRT_PI * mat * sigb * sigb * omega * omega));
    Expo = exp(-(S0 - s) * (S0 - s) / (2 * mat * SIGMA * SIGMA));
    p1 = exp(-(S0 - K1) * (S0 - K1) / (2 * mat * SIGMA1 * SIGMA1)) *
         pow(Csi1 * Csi1 - 2 * r * Csi1 + 1, -0.75) /
         (sqrt(2 * SRT_PI * mat * sigb * sigb * omega1 * omega1));
    p2 = exp(-(S0 - K2) * (S0 - K2) / (2 * mat * SIGMA2 * SIGMA2)) *
         pow(Csi2 * Csi2 - 2 * r * Csi2 + 1, -0.75) /
         (sqrt(2 * SRT_PI * mat * sigb * sigb * omega2 * omega2));

    /*		Expo = exp(-(S0-s)*(S0-s)/(2*mat*SIGMA*SIGMA))*
                    (1-(((XofZ*(exp(XofZ)+exp(-XofZ))/(exp(XofZ)-exp(-XofZ)))-1)/(XofZ*XofZ)+1)*mat*a*a/8);
    */
  }

  if (s <= K1) {
    integ = ((3 * p1 / (K1 * K1) - 6 * q / (K1 * K1 * K1)) * s * s +
             (6 * q / (K1 * K1) - 2 * p1 / K1) * s) *
            RV;
  }

  else if (s >= K2) {
    integ = p2 * exp(p2 * (K2 - s) / q) * RV;
  }

  else {
    integ = normalize * Expo * RV;
  }

  return integ;
}

double swp_f_cmsrateNewYork(double maturity, double underlying, double forward,
                            double ATMvol, double alpha, double beta,
                            double rho, double numperiod, double precision) {
  double betavol;
  double a = 0.001;
  double b = 1;
  double rate;
  double coverage;
  double mainpart;
  double K1, K2;
  double nstdev;
  double q;

  nstdev = 2.0;

  coverage = underlying / numperiod;

  srt_f_optblkvolATMtobetavolStochVol(forward, forward, ATMvol, maturity, alpha,
                                      beta, rho, &betavol);

  if (beta == 1) {
    K1 = exp(log(forward) - 0.5 * betavol * betavol * maturity -
             nstdev * betavol * sqrt(maturity));
    K2 = exp(log(forward) - 0.5 * betavol * betavol * maturity +
             nstdev * betavol * sqrt(maturity));
  }

  else {
    K1 = pow((1 - beta) * (pow(forward, 1 - beta) / (1 - beta) -
                           0.5 * beta * betavol * betavol *
                               pow(forward, beta - 1) * maturity -
                           nstdev * betavol * sqrt(maturity)),
             1 / (1 - beta));
    K2 = pow((1 - beta) * (pow(forward, 1 - beta) / (1 - beta) -
                           0.5 * beta * betavol * betavol *
                               pow(forward, beta - 1) * maturity +
                           nstdev * betavol * sqrt(maturity)),
             1 / (1 - beta));
  }

  mainpart = cms_sm_qsimp(&ProbaIntegrand, forward, maturity, betavol, alpha,
                          beta, rho, numperiod, coverage, K1, K2, precision);

  q = (1 - mainpart) / 2;

  rate = cms_sm_qsimp2(&CMSIntegrand, forward, maturity, betavol, alpha, beta,
                       rho, numperiod, coverage, K1, K2, q, a, b, precision);

  return rate;
}
