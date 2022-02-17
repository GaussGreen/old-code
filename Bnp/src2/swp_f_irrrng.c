/*****************************************************************************
        FILE_NAME:	swp_f_irrrng.c
******************************************************************************/
#include "math.h"
#include <num_h_allhdr.h"
#include <swp_h_all.h"
#include <swp_h_irrrng.h"

#define MAXITER 100

/* -------------------------------------------------------------------------
   Computes the PV of a series of cash-flows using the following formula:
        sum(i=1; i<len) { adj_cf[i]/(1+irr/comp)^i} /
   (1+irr/comp)^first_day_count where:
                - irr/comp is entered as the adjusted irr   ,
                - the cashflows are the coupons adjusted by the compounding
                - first_period_fraction represents the exact fraction of period
   from start to first coupon date computed using the coverage and the basis
   ------------------------------------------------------------------------- */

static double pv_with_irr(double adjusted_irr, double first_period_fraction,
                          double *cf, int len) {
  int i;
  double sum = 0;
  double l;

  l = -log(1.0 + adjusted_irr);

  /* Discount all the cash-flows up to the first coupon date (therefore i-1) */
  for (i = 1; i < len; i++)
    sum += cf[i] * exp(l * (i - 1));

  /* Discount from the first coupon date up to the first date  */
  sum *= exp(l * first_period_fraction);

  return sum;
}

/* ----------------------------------------------------------------------------
 */

Err pv_range_with_irr(Date *d, double *cf, int len, DateListType list_type,
                      Date prev_date, SrtBasisCode basis_code,
                      SrtCompounding compounding, double irr, double *pv) {
  double first_period_fraction;

  if (list_type == BROKEN) {
    if ((basis_code == BASIS_ACT_USD) && (prev_date != 0)) {
      first_period_fraction = day_count_date(d[0], d[1], basis_code);
      first_period_fraction /= day_count_date(prev_date, d[1], basis_code);
    } else {
      first_period_fraction = coverage(d[0], d[1], basis_code);
      first_period_fraction *= (double)compounding;
    }
  } else {
    first_period_fraction = 1.0;
  }

  *pv = pv_with_irr(irr, first_period_fraction, cf, len);

  return NULL;
}

/* ----------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------
   Returns irr for date d[0]  , cashflows occur at dates d[1]..d[len-1]
   and the whole structure is known to be worth pv at date d[0]
   (it can be a forward value)

   irr satisfies the equation:
                pv = sum[i=0;i<len] { cf[i]*cvg[i]/(1 + irr/comp)^(i + alpha) }
   where pv(irr) is computed by the above function
   ------------------------------------------------------------------------ */

Err irr_range(Date *d, double *cf, int len, DateListType list_type,
              Date prev_date, double pv, SrtBasisCode basis_code,
              SrtCompounding compounding, double *irr) {
  int count = 0;
  double first_period_fraction;
  double a[3], b[3];
  double nstop = 0.0, niter = 3.0, yans = 0.0;

  if (list_type == BROKEN) {
    if ((basis_code == BASIS_ACT_USD) && (prev_date != 0)) {
      first_period_fraction = day_count_date(d[0], d[1], basis_code);
      first_period_fraction /= day_count_date(prev_date, d[1], basis_code);
    } else {
      first_period_fraction = coverage(d[0], d[1], basis_code);
      first_period_fraction *= (double)compounding;
    }
  } else {
    first_period_fraction = 1.0;
  }

  a[0] = .06 / (double)compounding;
  b[0] = pv_with_irr(a[0], first_period_fraction, cf, len) - pv;
  a[1] = a[0] + .0001;
  b[1] = pv_with_irr(a[1], first_period_fraction, cf, len) - pv;
  a[2] = a[1] + .0001;

  while ((count < MAXITER) && (nstop < 1.0)) {
    b[2] = pv_with_irr(a[2], first_period_fraction, cf, len) - pv;
    newton(yans, niter, a, b, &nstop);
    count += 1;
  }

  if (count >= MAXITER)
    return serror("irr_range failed to converge.");

  *irr = a[2] * (double)compounding;

  return NULL;
}

/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------
    Computes the duration defined by:
    dur = 1/comp sum[i=0;i<len] { cf[i]*cvg[i]/(1 + irr/comp)^(i + alpha + 1) }
        (i.e. from -dPV/dIRR)
   ----------------------------------------------------------------------------
 */

Err duration_range(Date *d, double *cf, int len, DateListType list_type,
                   Date prev_date, SrtBasisCode basis_code,
                   SrtCompounding compounding, double irr, double *duration) {
  int i;
  double first_period_fraction, l;
  double sum;

  if (list_type == BROKEN) {
    if ((basis_code == BASIS_ACT_USD) && (prev_date != 0)) {
      first_period_fraction = day_count_date(d[0], d[1], basis_code);
      first_period_fraction /= day_count_date(prev_date, d[1], basis_code);
    } else {
      first_period_fraction = coverage(d[0], d[1], basis_code);
      first_period_fraction *= (double)compounding;
    }
  } else {
    first_period_fraction = 1.0;
  }

  l = -log(1.0 + irr / (double)compounding);

  sum = 0.0;
  for (i = 1; i < len; i++)
    sum += cf[i] * ((i - 1) + first_period_fraction) * exp(l * ((i - 1) + 1));

  sum *= exp(l * first_period_fraction);
  sum /= (double)compounding;

  *duration = sum;

  return NULL;
}
/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
    Computes the convexity defined by:
    conv = (1/comp)^2 sum[i=0;i<len] { cf[i]*cvg[i]*(cvg[i]+1)/(1 + irr/comp)^(i
   + alpha + 2) } (i.e. from d2PV/dIRR2)
   ----------------------------------------------------------------------- */

Err convexity_range(Date *d, double *cf, int len, DateListType list_type,
                    Date prev_date, SrtBasisCode basis_code,
                    SrtCompounding compounding, double irr, double *convexity) {
  int i;
  double first_period_fraction, l;
  double sum;

  if (list_type == BROKEN) {
    if ((basis_code == BASIS_ACT_USD) && (prev_date != 0)) {
      first_period_fraction = day_count_date(d[0], d[1], basis_code);
      first_period_fraction /= day_count_date(prev_date, d[1], basis_code);
    } else {
      first_period_fraction = coverage(d[0], d[1], basis_code);
      first_period_fraction *= (double)compounding;
    }
  } else {
    first_period_fraction = 1.0;
  }

  l = -log(1.0 + irr / (double)compounding);

  sum = 0.0;
  for (i = 1; i < len; i++)
    sum += cf[i] * ((i - 1) + first_period_fraction) *
           ((i - 1) + first_period_fraction + 1) * exp(l * ((i - 1) + 2));

  sum *= exp(l * first_period_fraction);
  sum /= (double)(compounding * compounding);

  *convexity = sum;

  return NULL;
}

/* -----------------------------------------------------------------------
    Computes the third moment defined by:
    conv = (1/comp)^3 sum[i=0;i<len] { cf[i]*cvg[i]*(cvg[i]+1)*(cvg[i]+2)/(1 +
   irr/comp)^(i + alpha + 3) } (i.e. from d3PV/dIRR3)
   ----------------------------------------------------------------------- */

Err thirdmoment_range(Date *d, double *cf, int len, DateListType list_type,
                      Date prev_date, SrtBasisCode basis_code,
                      SrtCompounding compounding, double irr,
                      double *thirdmoment) {
  int i;
  double first_period_fraction, l;
  double sum;

  if (list_type == BROKEN) {
    if ((basis_code == BASIS_ACT_USD) && (prev_date != 0)) {
      first_period_fraction = day_count_date(d[0], d[1], basis_code);
      first_period_fraction /= day_count_date(prev_date, d[1], basis_code);
    } else {
      first_period_fraction = coverage(d[0], d[1], basis_code);
      first_period_fraction *= (double)compounding;
    }
  } else {
    first_period_fraction = 1.0;
  }

  l = -log(1.0 + irr / (double)compounding);

  sum = 0.0;
  for (i = 1; i < len; i++)
    sum += cf[i] * ((i - 1) + first_period_fraction) *
           ((i - 1) + first_period_fraction + 1) *
           ((i - 1) + first_period_fraction + 2) * exp(l * ((i - 1) + 2));

  sum *= exp(l * first_period_fraction);
  sum /= (double)(compounding * compounding * compounding);

  *thirdmoment = sum;

  return NULL;
}

/* -------------------------------------------------------------------------------
    Computes a modified duration defined by:
    mod_dur = sum[i=1;i<len] { cf[i]*(i+alpha)/(1 + irr/comp)^(i+alpha) } / comp
   / PV (i.e. from -(1+IRR/comp)/PV * dPV/dIRR)
   -------------------------------------------------------------------------------
 */

Err modified_duration_range(Date *d, double *cf, int len,
                            DateListType list_type, Date prev_date,
                            SrtBasisCode basis_code, SrtCompounding compounding,
                            double irr, double *mod_duration) {
  int i;
  double first_period_fraction, l;
  double sum;
  double pv;

  if (list_type == BROKEN) {
    if ((basis_code == BASIS_ACT_USD) && (prev_date != 0)) {
      first_period_fraction = day_count_date(d[0], d[1], basis_code);
      first_period_fraction /= day_count_date(prev_date, d[1], basis_code);
    } else {
      first_period_fraction = coverage(d[0], d[1], basis_code);
      first_period_fraction *= (double)compounding;
    }
  } else {
    first_period_fraction = 1.0;
  }

  l = -log(1.0 + irr / (double)compounding);

  pv = 0.0;
  for (i = 1; i < len; i++)
    pv += cf[i] * exp(l * (i - 1));
  pv *= exp(l * first_period_fraction);

  sum = 0.0;
  for (i = 1; i < len; i++)
    sum += cf[i] * ((i - 1) + first_period_fraction) * exp(l * (i - 1));

  sum *= exp(l * first_period_fraction);
  sum /= (double)compounding;

  *mod_duration = sum / pv;

  return NULL;
}
/* ----------------------------------------------------------------------- */

/* -----------------------------------------------------------------------
    Computes a modified convexity defined by:
    conv = sum[i=1;i<len] { cf[i]*(i+alpha)*(i+alpha)/(1 + irr)^(i+alpha) } / PV
        (i.e. from (1+IRR)/PV * d[(1+IRR)*dPV/dIRR]/dIRR )
   ----------------------------------------------------------------------- */

Err modified_convexity_range(Date *d, double *cf, int len,
                             DateListType list_type, Date prev_date,
                             SrtBasisCode basis_code,
                             SrtCompounding compounding, double irr,
                             double *mod_convexity) {
  int i;
  double first_period_fraction, l;
  double sum;
  double pv;

  if (list_type == BROKEN) {
    if ((basis_code == BASIS_ACT_USD) && (prev_date != 0)) {
      first_period_fraction = day_count_date(d[0], d[1], basis_code);
      first_period_fraction /= day_count_date(prev_date, d[1], basis_code);
    } else {
      first_period_fraction = coverage(d[0], d[1], basis_code);
      first_period_fraction *= (double)compounding;
    }
  } else {
    first_period_fraction = 1.0;
  }

  l = -log(1.0 + irr / (double)compounding);

  pv = 0.0;
  for (i = 1; i < len; i++)
    pv += cf[i] * exp(l * (i - 1));
  pv *= exp(l * first_period_fraction);

  sum = 0.0;
  for (i = 1; i < len; i++)
    sum += cf[i] * ((i - 1) + first_period_fraction) *
           ((i - 1) + first_period_fraction) * exp(l * (i - 1));

  sum *= exp(l * first_period_fraction);
  sum /= (double)(compounding * compounding);

  *mod_convexity = sum / pv;

  return NULL;
}

#undef MAXITER
