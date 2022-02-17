/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT  , Fixed Income 2020 Addins */
/*      SUB_SYSTEM:     BDT     BOND TOOLS                                    */
/*                                                                            */
/*      MODULE NAME:    SWP_F_BOND_COMPUTE                                    */
/*                                                                            */
/*      PURPOSE:        VARIOUS BOND VALUING FUNCTIONS                        */
/*                                                                            */
/*      AUTHORS:        Julia Matsumoto & ????                         	      */
/*                                                                            */
/*      DATE:                                                                 */
/*                                                                            */
/*      VERSION:        01                                                    */
/*                                                                            */
/*      DESCRIPTION:    XX                                                    */
/*                                                                            */
/*      FUNCTIONS USED: XXX_X_XXXXXXXXX                                       */
/*              Must include all imported function call made by the module    */
/*                                                                            */
/*      PARAMETERS:     <not applicable>                                      */
/*                                                                            */
/*      RETURNS:                                                              */
/*                                                                            */
/*      DATA ACCESSED:                                                        */
/*              Tables  , files  , and global variables accessed. */
/*                                                                            */
/******************************************************************************/
/*                      Amendment History                                     */
/******************************************************************************/
/*                                                                            */
/*      AMENDED BY:     Julia                                                 */
/*                                                                            */
/*      DATE:                                                                 */
/*                                                                            */
/*                                                                            */
/*      REASON:         Merge with Swaptools                                  */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
/*                                                                            */
/*      AMENDED BY:     Jasbir S Malhi                                        */
/*                                                                            */
/*      DATE:           3RD JUNE  , 1994 */
/*                                                                            */
/*                                                                            */
/*      REASON:         Use WITH MAP                                          */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
/*	AMENDED BY:	Alex Benech					      */
/*									      */
/*	DATE:		23RD JUNE  , 1994 */
/*                                                                            */
/*	REASON:		Add new_set_price_fct function			      */
/*									      */
/******************************************************************************/
/*	AMENDED BY:	Julia Matsumoto					      */
/*									      */
/*	DATE:		13 Jan  , 1995					      */
/*                                                                            */
/*	REASON:		faster goalseek_10y     			      */
/*									      */
/***  include files  **********************************************************/

#include "math.h"
#include "num_h_allhdr.h"
#include "swp_h_all.h"

#define ACCURACY 0.000001

/*** declare static variables  ************************************************/

static double mem_rate[50];

/***  functions  **************************************************************/

/* function to calculate the dirty price of a bond     */

double dirty_price_fct(SwapDP p, double coupon, double irr, double redemption,
                       double first_coupon) {
  double price, disc1;
  double dr;
  DateList list;

  /* list of coupon date from maturity to value date */
  list = SwapDP_to_DateList(&p, NO_BUSDAY_CONVENTION);

  disc1 = 1 / (1 + irr / (double)p.compd);

  /* we price the bond from the second coupon (the first 'regular' coupon) */
  /* to maturity */
  if (disc1 != 1) {
    price = (coupon / (double)p.compd) * (pow(disc1, list.len - 2) - 1) /
                (disc1 - 1) +
            redemption * pow(disc1, list.len - 3);
  } else {
    price = (coupon / (double)p.compd) * (list.len - 2) + redemption;
  }

  /* First coupon */
  price = price * disc1 + first_coupon / (double)p.compd;

  /* we have a different discount factor for the first coupon */
  /* if there is an odd first coupon or not */
  if (list.type == BROKEN) {
    if (p.basis_code == BASIS_ACT_USD) {
      dr = day_count_date(list.date[0], list.date[1], p.basis_code);
      dr /= day_count_date(list.prev, list.date[1], p.basis_code);
      dr /= p.compd;
    } else {
      dr = coverage(list.date[0], list.date[1], p.basis_code);
    }

    price *= pow(disc1, p.compd * dr);
  } else {
    price *= disc1;
  }

  srt_free(list.date);

  return (price);
}

/***  Accrued interest function   *********************************************/

double acc_int_fct(SwapDP p, double coupon) {
  DateList list;
  double acc, dr;

  list = SwapDP_to_DateList(&p, NO_BUSDAY_CONVENTION);

  if (list.type == BROKEN) {
    if (p.basis_code == BASIS_ACT_USD) {
      dr = day_count_date(list.prev, list.date[0], p.basis_code);
      dr /= day_count_date(list.prev, list.date[1], p.basis_code);
      dr /= p.compd;
    } else
      dr = coverage(list.prev, list.date[0], p.basis_code);

    acc = coupon * dr;
  } else
    acc = 0;

  srt_free(list.date);

  return acc;
}

/***  Clean Price Function   **************************************************/

double clean_price_fct(SwapDP p, double coupon, double irr, double redemption,
                       double first_coupon) {
  double price;

  price = dirty_price_fct(p, coupon, irr, redemption, first_coupon);
  price -= acc_int_fct(p, first_coupon);

  return (price);
}

/***  Yield Function   **************************************************/

double yield_fct(SwapDP p, double coupon, double clean_price, double redemption,
                 double first_coupon) {
  double dirty_price;
  int i;

  double guess, delta1 = 10, delta2, df, f, marg1 = 0, marg2 = 1, xl, xh, fl,
                fh, dx, dxold, temp, delta = 0.0001;
  DateList list;

  list = SwapDP_to_DateList(&p, NO_BUSDAY_CONVENTION);

  if (p.start >= p.end)
    return (coupon);

  dirty_price = clean_price + acc_int_fct(p, first_coupon);

  fl =
      dirty_price_fct(p, coupon, marg1, redemption, first_coupon) - dirty_price;
  fh =
      dirty_price_fct(p, coupon, marg2, redemption, first_coupon) - dirty_price;

  if (fl * fh >= 0)
    return (0);

  if (fl < 0) {
    xl = marg1;
    xh = marg2;
  } else {
    xl = marg2;
    xh = marg1;
  }

  guess = (xh + xl) / 2;
  dxold = fabs(marg1 - marg2);
  dx = dxold;

  f = delta1 =
      dirty_price_fct(p, coupon, guess, redemption, first_coupon) - dirty_price;
  delta2 = dirty_price_fct(p, coupon, guess + delta, redemption, first_coupon) -
           dirty_price;
  df = (delta2 - delta1) / delta;

  for (i = 0; i < 500; i++) {
    if ((((guess - xh) * df - f) * ((guess - xl) * df - f) >= 0.0) ||
        (fabs(2.0 * f) > fabs(dxold * df))) {
      /* Bisect if Newton out of range  , or not
      decreasing fast enough */
      dxold = dx;
      dx = (xh - xl) / 2;
      guess = xl + dx;
      if (xl == guess)
        break; /* Change in root is neglible */
    } else     /* Newton step acceptable. Take it. */
    {
      dxold = dx;
      dx = f / df;
      temp = guess;
      guess -= dx;
      if (temp == guess)
        break;
    }

    if (fabs(dx) < ACCURACY)
      break; /* Convergence criterion */

    f = delta1 = dirty_price_fct(p, coupon, guess, redemption, first_coupon) -
                 dirty_price;
    delta2 =
        dirty_price_fct(p, coupon, guess + delta, redemption, first_coupon) -
        dirty_price;
    df = (delta2 - delta1) / delta;

    if (f < 0)
      xl = guess;
    else
      xh = guess;
  }

  srt_free(list.date);

  return (guess);
}

/* function to calculate the clean price at a new settlement date */
double new_set_price_fct(Date acc_date, Date new_date, SwapDP p, double coupon,
                         double clean_price, String swap, String repo,
                         double redemption, double first_coupon) {
  double dirty_price, price;
  SwapDP acc_p;
  DateList list;

  acc_p = p;
  acc_p.start = acc_date;
  dirty_price = clean_price + acc_int_fct(acc_p, first_coupon);

  price = dirty_price / swp_f_df((Ddate)(p.start), (Ddate)new_date, repo);

  acc_p.start = new_date;

  /* test if a cash flow occurs between the value date and the new date*/
  list = SwapDP_to_DateList(&p, NO_BUSDAY_CONVENTION);
  if (list.date[1] < new_date) {
    price -=
        first_coupon / swp_f_df((Ddate)(list.date[1]), (Ddate)new_date, swap);
    price -= acc_int_fct(acc_p, coupon);
  } else {
    price -= acc_int_fct(acc_p, first_coupon);
  }

  srt_free(list.date);

  return (price);
}

/* bond sensibility */
double bond_sens_fct(SwapDP p, double coupon, double irr, double redemption,
                     double first_coupon) {
  double price1, price2;
  double delta;

  price1 = dirty_price_fct(p, coupon, irr, redemption, first_coupon);
  price2 = dirty_price_fct(p, coupon, irr + 0.0001, redemption, first_coupon);
  delta = (price1 - price2) * 10000;

  return delta;
}

/* Forward clean price a bond */
double fwd_clean_price_fct(Date fut, SwapDP p, double coupon,
                           double clean_price, String swap, String repo,
                           double first_coupon) {
  int nc, i;
  double dirty_price, fwd_clean_price, acc_int, fwd_acc_int;
  double sum_coup = 0.0, c, df;
  DateList list;
  SwapDP pf;

  list = SwapDP_to_DateList(&p, NO_BUSDAY_CONVENTION);
  acc_int = acc_int_fct(p, first_coupon);
  dirty_price = clean_price + acc_int;

  pf = p; /*Mainly to have pf.end = p.end :
                this alllows us to define a proper BKWD direction
                for the forward bond dates and then get the exact
                forward accrued interests */
  pf.start = fut;
  pf.first_full_fixing = fut;
  pf.direction = BKWD;

  c = first_coupon;

  nc = 0;
  while (list.date[nc + 1] <= fut)
    nc++;

  for (i = 1; i <= nc; i++) {
    sum_coup +=
        (c / p.compd) / swp_f_df((Ddate)(list.date[i]), (Ddate)fut, swap);
    c = coupon;
  }

  fwd_acc_int = acc_int_fct(pf, c);

  df = swp_f_df((Ddate)(list.date[0]), (Ddate)fut, repo);

  fwd_clean_price = dirty_price / df - sum_coup - fwd_acc_int;

  srt_free(list.date);

  return (fwd_clean_price);
}

/* Get today's clean price from a future price at a future date */
double tdy_price_fct(Date fut, SwapDP p, double coupon, double fwd_price,
                     String swap, String repo, double first_coupon)

{
  int i;
  double guess, delta1 = 10, delta2, df, f, marg1 = 0, marg2 = 1, xl, xh, fl,
                fh, dx, dxold, temp, delta = 0.0001;

  marg2 = fwd_price * 3;

  fl = fwd_clean_price_fct(fut, p, coupon, marg1, swap, repo, first_coupon) -
       fwd_price;
  fh = fwd_clean_price_fct(fut, p, coupon, marg2, swap, repo, first_coupon) -
       fwd_price;

  if (fl * fh >= 0)
    return (0);

  if (fl < 0) {
    xl = marg1;
    xh = marg2;
  } else {
    xl = marg2;
    xh = marg1;
  }

  guess = (xh + xl) / 2;
  dxold = fabs(marg1 - marg2);
  dx = dxold;

  f = delta1 =
      fwd_clean_price_fct(fut, p, coupon, guess, swap, repo, first_coupon) -
      fwd_price;
  delta2 = fwd_clean_price_fct(fut, p, coupon, guess + delta, swap, repo,
                               first_coupon) -
           fwd_price;
  df = (delta2 - delta1) / delta;

  for (i = 0; i < 500; i++) {
    if ((((guess - xh) * df - f) * ((guess - xl) * df - f) >= 0.0) ||
        (fabs(2.0 * f) >
         fabs(dxold * df))) { /* Bisect if Newton out of range  , or not
                                    decreasing fast enough */
      dxold = dx;
      dx = (xh - xl) / 2;
      guess = xl + dx;

      if (xl == guess)
        break; /* Change in root is neglible */
    } else     /* Newton step acceptable. Take it. */
    {
      dxold = dx;
      dx = f / df;
      temp = guess;
      guess -= dx;

      if (temp == guess)
        break;
    }

    if (fabs(dx) < ACCURACY)
      break; /* Convergence criterion */

    f = delta1 =
        fwd_clean_price_fct(fut, p, coupon, guess, swap, repo, first_coupon) -
        fwd_price;
    delta2 = fwd_clean_price_fct(fut, p, coupon, guess + delta, swap, repo,
                                 first_coupon) -
             fwd_price;
    df = (delta2 - delta1) / delta;

    if (f < 0)
      xl = guess;
    else
      xh = guess;
  }

  return (guess);
}

/* Forward yield function */
double fwd_irr_fct(Date fut, SwapDP p, double coupon, double clean_price,
                   String swap, String repo, double redemption,
                   double first_coupon) {
  double fwd_clean_price, fwd_yield;
  SwapDP pf;

  pf = p;
  pf.start = fut;
  pf.direction = BKWD;

  fwd_clean_price = fwd_clean_price_fct(fut, p, coupon, clean_price, swap, repo,
                                        first_coupon);

  fwd_yield = yield_fct(pf, coupon, fwd_clean_price, redemption, first_coupon);

  return (fwd_yield);
}

/* Yield function with va_list */
double yield_function(double x, va_list argptr) {
  SwapDP pf;
  double c, r, f;
  pf = va_arg(argptr, SwapDP);
  c = va_arg(argptr, double);
  r = va_arg(argptr, double);
  f = va_arg(argptr, double);
  return (yield_fct(pf, c, x, r, f));
}

/* Price function with va_list */
double price_function(double x, va_list argptr) {

  SwapDP pf;
  double c, r, f, ret;

  pf = va_arg(argptr, SwapDP);
  c = va_arg(argptr, double);
  r = va_arg(argptr, double);
  f = va_arg(argptr, double);

  ret = clean_price_fct(pf, c, x, r, f);
  return (ret);
}

double prim(double (*function)(double x, va_list), double value0, ...)
/*additional parameters: SwapDP p  ,double c  , double r  , double f*/
{
  va_list argptr;
  double delta = 0.0001;
  double ret;
  va_start(argptr, value0);
  ret = (function(value0 + delta / 2, argptr) -
         function(value0 - delta / 2, argptr)) /
        delta;
  va_end(argptr);
  return (ret);
}

double prim2(double (*function)(double, va_list), double value0, ...)
/*additional parameters: SwapDP p  ,double c  , double r  , double f*/
{
  va_list argptr;
  double delta = 0.0001;
  double ret;
  va_start(argptr, value0);
  ret = 4 *
        (function(value0 + delta / 2, argptr) +
         function(value0 - delta / 2, argptr) - 2 * function(value0, argptr)) /
        (delta * delta);
  va_end(argptr);
  return (ret);
}

double fwd_irr_fct_tlr(Date fut, SwapDP p, double coupon, double clean_price,
                       double volatility, String swap, String repo,
                       double redemption, double first_coupon) {
  double yield_prim = 0, yield_prim2 = 0, T, exp_val;
  double fwd_clean_price, fwd_yield, stoch;
  SwapDP pf;

  pf = p;
  pf.start = fut;
  pf.direction = BKWD;

  T = coverage(p.start, fut, p.basis_code);
  fwd_clean_price = fwd_clean_price_fct(fut, p, coupon, clean_price, swap, repo,
                                        first_coupon);
  fwd_yield = yield_fct(pf, coupon, fwd_clean_price, redemption, first_coupon);
  /*lib$signal(SS$_DEBUG);*/
  yield_prim2 = prim2(yield_function, fwd_clean_price, pf, coupon, redemption,
                      first_coupon);
  stoch = fwd_clean_price * fwd_clean_price *
          (exp(volatility * volatility * T) - 1) * yield_prim2 / 2;
  exp_val = fwd_yield + stoch;
  return (exp_val);
}

double fct_to_integrate(double x, va_list argptr)
/*SwapDP pf  ,double coupon  , double clean_price  ,double volatility  ,
                     double myu  , int T  ,
                     double redemption  , double first_coupon*/
{
  SwapDP p;
  double coupon, redemption, first_coupon, myu, clean_price, volatility;
  double T;
  double y;

  p = va_arg(argptr, SwapDP);
  coupon = va_arg(argptr, double);
  clean_price = va_arg(argptr, double);
  volatility = va_arg(argptr, double);
  myu = va_arg(argptr, double);
  T = va_arg(argptr, double);
  redemption = va_arg(argptr, double);
  first_coupon = va_arg(argptr, double);

  y = yield_fct(p, coupon, clean_price * exp(x), redemption, first_coupon);
  return (1 / (volatility * sqrt(2 * SRT_PI * T)) * clean_price * y * exp(x) *
          exp(-pow(x - (myu - volatility * volatility / 2) * T, 2) /
              (2 * volatility * volatility * T)));
}

double fct_to_integrate2(double x, va_list argptr)
/*SwapDP pf  ,double coupon  , double clean_price  ,double volatility  ,
                double myu  , double T  ,
                double redemption  , double first_coupon*/
{
  SwapDP pf;
  double coupon, redemption, first_coupon, myu, clean_price, volatility;
  double T;
  double s, n, d;

  pf = va_arg(argptr, SwapDP);
  coupon = va_arg(argptr, double);
  clean_price = va_arg(argptr, double);
  volatility = va_arg(argptr, double);
  myu = va_arg(argptr, double);
  T = va_arg(argptr, double);
  redemption = va_arg(argptr, double);
  first_coupon = va_arg(argptr, double);

  s = clean_price_fct(pf, coupon, x, redemption, first_coupon);

  d = prim(price_function, x, pf, coupon, redemption, first_coupon);
  n = 1 / (sqrt(2 * SRT_PI * T) * volatility) *
      exp(-pow(log(s / clean_price) - (myu - volatility * volatility / 2) * T,
               2) /
          (2 * volatility * volatility * T));
  return (x * d * n / s);
}
/*****************************************************************************/

double fwd_irr_fct_smp(Date fut, SwapDP p, double coupon, double clean_price,
                       double volatility, String swap, String repo,
                       double redemption, double first_coupon) {
  double exp_val;
  double T;
  double fwd_clean_price, fwd_dirty_price, myu;
  SwapDP pf;

  pf = p;
  pf.start = fut;
  pf.direction = BKWD;

  T = coverage(p.start, fut, p.basis_code);
  fwd_clean_price = fwd_clean_price_fct(fut, p, coupon, clean_price, swap, repo,
                                        first_coupon);
  fwd_dirty_price = fwd_clean_price_fct(fut, p, coupon, clean_price, swap, repo,
                                        first_coupon) +
                    acc_int_fct(pf, coupon);
  myu = log(fwd_clean_price / clean_price) / T;
  /*x0=log(clean_price_fct(p  ,coupon  ,-0.2  ,redemption  ,first_coupon)
          /clean_price);*/
  exp_val = sm_qsimp_list(fct_to_integrate2, -0.5, 1, pf, coupon, clean_price,
                          volatility, myu, T, redemption, first_coupon);

  /*fwd_yield=yield_fct(pf  ,coupon  ,fwd_clean_price  ,redemption  ,
          first_coupon);*/
  return (-exp_val);
}
