
#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/*******************************************************************************
 *
 * FUNCTION     : ??
 * PURPOSE      : ??
 * DESCRIPTION  : ??
 * PARAMETERS   : ??
 * RETURNS      : ??
 *****************************************************************************/

/*******************************************************************************
        STAY_INSIDE_FCT
        Calculates the probability of staying within two barrier levels
        between today and maturity        ,
        for either risk neutral (measure_flag = 0) or measure associated
        with the underlying (measure_flag = 1) .

********************************************************************************/

double stay_inside_fct(double spot, double fwd, double up, double low,
                       double mat, double vol, int n_terms, int measure_flag) {
  double b_up;
  double b_do;
  double mu;
  double mult;
  double prob;

  if ((mat <= 0.0) || (vol == 0.0)) {
    if ((fwd > up) || (fwd < low))
      return 0.0;
    else
      return 1.0;
  }

  b_up = log(up / spot);
  b_do = log(low / spot);

  mult = (measure_flag == 1) ? 1 : -1;
  mu = log(fwd / spot) / mat + 0.5 * mult * vol * vol;

  prob = dbarrier_prob(b_do, b_up, b_do, mu, vol, mat, SRT_CALL, n_terms);

  return prob;
} /* END stay_inside_fct() */

/*******************************************************************************
        STAY_INSIDE_FCT_FWD
        Calculates the probability of staying within two barrier levels
        between time T1 (i.e. forward start) and maturity        ,
        for either risk neutral (measure_flag = 0) or measure associated
        with the underlying (measure_flag = 1) .

********************************************************************************/

double stay_inside_fwd_fct(double spot, double fwd1, double fwd2, double up,
                           double low, double t1, double t2, double vol0,
                           double vol, int n, int measure_flag) {

  double price1 = 0;
  double price2 = 0;
  double price = 0;
  double p1, p2, p3, p4;
  double mu0, mu1, mt1, mt2, mv1, mv2, rt, vt1, vt2, sl, su, ul, q1;
  int i, drift_flag;

  if (measure_flag == 1) {
    drift_flag = -1;
  } else {
    drift_flag = 1;
  }
  mu0 = log(fwd1 / spot) / t1 - drift_flag * vol0 * vol0 / 2;
  mu1 = log(fwd2 / fwd1) / (t2 - t1) - drift_flag * vol * vol / 2;
  mt1 = mu0 * t1;
  mt2 = mt1 + mu1 * (t2 - t1);
  rt = vol0 * sqrt(t1 / (vol * vol * (t2 - t1) + vol0 * vol0 * t1));
  vt1 = vol0 * sqrt(t1);
  vt2 = sqrt(vol * vol * (t2 - t1) + vol0 * vol0 * t1);
  sl = log(spot / low);
  su = log(spot / up);
  ul = log(up / low);
  mv1 = (mu0 - 2 * mu1 * vol0 * vol0 / vol / vol) * t1;
  mv2 = mt1 - mu1 * (t2 - t1 + 2 * vol0 * vol0 / vol / vol * t1);

  for (i = -n; i <= n; i++) {
    q1 = exp(log(up / low) * (2 * i * mu1 / vol / vol));
    p1 = -bivar((sl + mt1) / vt1, (su + 2 * i * ul + mt2) / vt2, rt);
    p2 = bivar((su + mt1) / vt1, (su + 2 * i * ul + mt2) / vt2, rt);
    p3 = bivar((sl + mt1) / vt1, (sl + 2 * i * ul + mt2) / vt2, rt);
    p4 = -bivar((su + mt1) / vt1, (sl + 2 * i * ul + mt2) / vt2, rt);
    if (fabs(p1 + p2 + p3 + p4) >= EPS) {
      price1 += q1 * (p1 + p2 + p3 + p4);
    }
    p1 = -bivar((sl + mv1) / vt1, (su - ul + 2 * i * ul + mv2) / vt2, rt);
    p2 = bivar((su + mv1) / vt1, (su - ul + 2 * i * ul + mv2) / vt2, rt);
    p3 = bivar((sl + mv1) / vt1, (su + 2 * i * ul + mv2) / vt2, rt);
    p4 = -bivar((su + mv1) / vt1, (su + 2 * i * ul + mv2) / vt2, rt);
    if (fabs(p1 + p2 + p3 + p4) >= EPS) {
      price2 += (p1 + p2 + p3 + p4) / q1;
    }
  }

  price = price1 - exp(log(up / spot) * (2 * mu1 / vol / vol)) *
                       exp(2 * mu1 * t1 *
                           (vol0 * vol0 / vol / vol * mu1 - mu0) / vol / vol) *
                       price2;

  return (price);
}

/********** Probability to touch b_yes and not b_no : lognormal ***************/

double hit_one_barrier_fct(double spot, double fwd, double b_yes, double b_no,
                           double mat, double vol, int n, int measure_flag)

{
  double a, b, price, yn, coeff_tot, coeff;
  double mu, term, d1, d2;
  int i, mult, mult2;

  if ((mat <= 0.0) || (vol == 0.0)) {
    if ((b_yes > b_no) && (fwd >= b_yes))
      return 1.0;
    else if ((b_yes < b_no) && (fwd <= b_yes))
      return 1.0;
    else
      return 0.0;
  }

  if (b_yes < b_no) {
    if (b_yes >= spot)
      return (1); /*look at this later*/
    if (b_no <= spot)
      return (0);
    mult2 = 1;
  } else {
    if (b_yes <= spot)
      return (1); /*look at this later*/
    if (b_no >= spot)
      return (0);
    mult2 = -1;
  }

  a = log(b_yes / spot);
  b = log(spot / b_no);
  mu = log(fwd / spot) / mat - vol * vol / 2;
  price = 0.0;

  for (i = -n; i <= n; i++) {
    yn = -(2 * i + 1) * (a + b) + b;
    if (yn > 0)
      mult = 1;
    else
      mult = -1;
    coeff = exp(-2 * mu * yn / vol / vol);
    d1 = mult * (-yn - mu * mat) / vol / sqrt(mat);
    d2 = mult * (-yn + mu * mat) / vol / sqrt(mat);
    coeff_tot = exp(-mu * 2 * i * (a + b) / vol / vol);
    if (yn == 0)
      term = -0.5;
    else
      term = mult * (norm(d1) + coeff * norm(d2));
    price += term * coeff_tot;
  }

  return (mult2 * price);
} /* END hit_one_barrier_fct */

double fct_gold_stay_inside(double x, va_list argptr) {
  double spot, fwd, gap, mat, vol;
  double mu, b_up, b_do;
  int n, ff, mult;
  SrtCallPutType call_put = SRT_CALL;
  static double result;

  spot = va_arg(argptr, double);
  fwd = va_arg(argptr, double);
  gap = va_arg(argptr, double);
  mat = va_arg(argptr, double);
  vol = va_arg(argptr, double);
  n = va_arg(argptr, int);
  ff = va_arg(argptr, int);

  b_up = log((x + gap) / spot);
  b_do = log(x / spot);

  mult = (ff == 1) ? 1 : -1;
  mu = log(fwd / spot) / mat + 0.5 * mult * vol * vol;

  /* A trick is used: strike = barrier_down        , with a CALL */
  result = -dbarrier_prob(b_do, b_up, b_do, mu, vol, mat, SRT_CALL, n);

  return (result);
}

double fct_gold_stay_inside_fwd(double x, va_list argptr) {
  double spot, fwd1, fwd2, gap, t1, t2, vol0, vol;
  int n, measure_flag;
  static double result;

  spot = va_arg(argptr, double);
  fwd1 = va_arg(argptr, double);
  fwd2 = va_arg(argptr, double);
  gap = va_arg(argptr, double);
  t1 = va_arg(argptr, double);
  t2 = va_arg(argptr, double);
  vol0 = va_arg(argptr, double);
  vol = va_arg(argptr, double);
  n = va_arg(argptr, int);
  measure_flag = va_arg(argptr, int);

  result = -stay_inside_fwd_fct(spot, fwd1, fwd2, x + gap, x, t1, t2, vol0, vol,
                                n, measure_flag);

  return (result);
}

/* --------------------------------------------------------------------------------------------------------------------------------
 prob_ExpectedDays()
 BBFlag = 0 use Brownian Bridge ( dParameters not used )
                = 1 use no smile model ( dParameters not used )
                = 2 use linear smile model ( dParameters = dSlope )

 dBarrier is a 2 element array with the up barrier as the first element
--------------------------------------------------------------------------------------------------------------------------------
*/
double prob_ExpectedDays_BrownianBridge(double x1, double x0, double dBup,
                                        double dBdown, double dVol,
                                        double dStartTime, double dLength,
                                        int nStep);
double prob_ExpectedDays_NoSmile(double x1, double x0, double dBup,
                                 double dBdown, double dVol, double dStartTime,
                                 double dLength, int nStep, int Flag);
double prob_ExpectedDays_LinearSmile(double dDrift, double dSpot,
                                     double *dBarrier, double dATMVol,
                                     double *dSlope, double dStartTime,
                                     double dLength, int nStep);

double prob_ExpectedDays(double x1, double x0, double *dBarrier, double dVol,
                         double dStartTime, double dLength, int nStep,
                         double *dParameters, int BBflag) {
  if (BBflag == 0)
    return prob_ExpectedDays_BrownianBridge(x1, x0, dBarrier[0], dBarrier[1],
                                            dVol, dStartTime, dLength, nStep);
  else if (BBflag == 1)
    return prob_ExpectedDays_NoSmile(x1, x0, dBarrier[0], dBarrier[1], dVol,
                                     dStartTime, dLength, nStep, BBflag);
  else if (BBflag == 2) {
    return prob_ExpectedDays_LinearSmile(x1, x0, dBarrier, dVol, dParameters,
                                         dStartTime, dLength, nStep);
  } else
    return -1.0;
}

/* --------------------------------------------------------------------------------------------------------------------------------
        prob_ExpectedDays_BB()
        Calculates the expected number of days within a corridor using a
Brownian Bridge
--------------------------------------------------------------------------------------------------------------------------------
*/
double prob_ExpectedDays_BrownianBridge(double x1, double x0, double dBup,
                                        double dBdown, double dVol,
                                        double dStartTime, double dLength,
                                        int nStep) {

  double dLogBup, dLogBdn, dLogF, dLogS, dStepSize, dExpectDays, dTime,
      dTimeRatio, dTotalTime;
  int i;

  dLogBup = log(dBup);
  dLogBdn = log(dBdown);
  dLogF = log(x1);
  dLogS = log(x0);

  dStepSize = dLength / nStep;
  dTotalTime = dLength + dStartTime;

  if (dStartTime == 0.0) {
    dExpectDays = ((x0 < dBup && x0 > dBdown) ? 1.0 : 0.0);
    i = 1;
  } else {
    dExpectDays = 0.0;
    i = 0;
  }

  for (; i < nStep; i++) {
    dTime = dStartTime + i * dStepSize;
    dTimeRatio = dTime / dTotalTime;
    dExpectDays +=
        norm((dLogBup - (1.0 - dTimeRatio) * dLogS - dTimeRatio * dLogF) /
             dVol / sqrt(dTime * (1 - dTimeRatio))) -
        norm((dLogBdn - (1.0 - dTimeRatio) * dLogS - dTimeRatio * dLogF) /
             dVol / sqrt(dTime * (1 - dTimeRatio)));
  }

  return dExpectDays * dLength * 365.0 / (double)nStep;
}

/* --------------------------------------------------------------------------------------------------------------------------------
        prob_ExpectedDays_T0()
        Calculates the expected number of days within a corridor using call
spreads
--------------------------------------------------------------------------------------------------------------------------------
*/

double prob_ExpectedDays_NoSmile(double x1, double x0, double dBup,
                                 double dBdown, double dVol, double dStartTime,
                                 double dLength, int nStep, int Flag) {

  double dLogModBup, dLogModBdn, dStepSize, dExpectDays, dTime, dStdDev, dDrift;
  int i;

  dLogModBup = log(dBup / x0);
  dLogModBdn = log(dBdown / x0);
  if (Flag == -1)
    dDrift = x1;
  else
    dDrift = log(x1 / x0) / dLength - 0.5 * dVol * dVol;

  dStepSize = dLength / nStep;

  if (dStartTime == 0.0) {
    dExpectDays = ((x0 < dBup && x0 > dBdown) ? 1.0 : 0.0);
    i = 1;
  } else {
    dExpectDays = 0.0;
    i = 0;
  }

  for (; i < nStep; i++) {
    dTime = dStartTime + i * dStepSize;
    dStdDev = dVol * sqrt(dTime);
    dExpectDays += norm((dLogModBup - dDrift * dTime) / dStdDev) -
                   norm((dLogModBdn - dDrift * dTime) / dStdDev);
  }

  return dExpectDays * dLength * 365.0 / (double)nStep;
}

/* --------------------------------------------------------------------------------------------------------------------------------
        prob_ExpectedDays_Smile()
        Calculates the expected number of days within a corridor using call
spreads with a simple linear smile model The 0 value is the up and the 1 value
is down
--------------------------------------------------------------------------------------------------------------------------------
*/

double prob_ExpectedDays_LinearSmile(double dDrift, double dSpot,
                                     double *dBarrier, double dATMVol,
                                     double *dSlope, double dStartTime,
                                     double dLength, int nStep) {

  double dStepSize, dExpectDays, dTime, dStdDev, dLog[2], dCallSpread[2],
      dVol[2], dHplusTerm[2], dVegaTerm[2];
  double dHplus, dHminus, dDriftTerm;
  int i, j;

  for (j = 0; j < 2; j++) {
    dLog[j] = log(dSpot / dBarrier[j]);
    dVol[j] = dATMVol + dSlope[j] * dBarrier[j] / dSpot;
    dHplusTerm[j] = dDrift + 0.5 * dVol[j] * dVol[j];
    dVegaTerm[j] = dSpot * INV_SQRT_TWO_PI * dSlope[j];
  }

  dStepSize = dLength / nStep;

  if (dStartTime == 0.0) {
    dExpectDays = ((dSpot < dBarrier[0] && dSpot > dBarrier[1]) ? 1.0 : 0.0);
    i = 1;
  } else {
    dExpectDays = 0.0;
    i = 0;
  }

  for (; i < nStep; i++) {
    dTime = dStartTime + i * dStepSize;
    dDriftTerm = dDrift * dTime;
    for (j = 0; j < 2; j++) {
      dStdDev = dVol[j] * sqrt(dTime);
      dHplus = (dLog[j] + dHplusTerm[j] * dTime) / dStdDev;
      dHminus = dHplus - dStdDev;
      dCallSpread[j] = -norm(dHminus) +
                       dVegaTerm[j] * exp(dDriftTerm - 0.5 * dHplus * dHplus);
    }
    dExpectDays += dCallSpread[1] - dCallSpread[0];
  }

  return dExpectDays * dLength * 365.0 / (double)nStep;
}
