
#include "AmortMidatProdStruct.h"
#include "AmortMidatCalib.h"
#include "EuroAmortSwaption.h"
#include "LGM2Fgrfn.h"
#include "LGMSVUtil.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srt_h_lgmtypes.h"
#include "swp_h_amortswaption.h"
#include "swp_h_convslidingcorrel.h"

#define BIG_BIG_NUMBER 1E40

/*	Functions for the funding leg */

Err am_fill_fund_leg(
    /*	Coupons that started before today are disregarded */
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I  , 1: E */
    double *fund_not, int fund_ncpn, long *fund_fix, long *fund_start,
    long *fund_end, long *fund_pay, char **fund_basis, double *fund_spr,
    double *fund_mrg, AM_FUND_LEG fund_leg) {
  int i, j;
  SrtBasisCode bas;
  AM_FUND_CPN cpn;
  Err err = NULL;

  /*	Initialise pointers to NULL */
  fund_leg->notional = NULL;
  fund_leg->cpn = NULL;
  fund_leg->margin = NULL;
  fund_leg->spread = NULL;

  /*	Skip coupons fixed before today */
  i = 0;
  while (i < fund_ncpn && fund_fix[i] < today + eod_flag) {
    i++;
  }

  /*	Check that at least one coupon is left */
  if (i == fund_ncpn) {
    /*	err = "All funding coupons start before today in am_fill_fund_leg"; */
    fund_leg->num_cpn = 0;
    goto FREE_RETURN;
  }

  /*	Allocate memory */
  fund_leg->num_cpn = fund_ncpn - i;
  fund_leg->cpn = (am_fund_cpn *)calloc(fund_leg->num_cpn, sizeof(am_fund_cpn));
  fund_leg->notional = (double *)calloc(fund_leg->num_cpn, sizeof(double));
  fund_leg->margin = (double *)calloc(fund_leg->num_cpn, sizeof(double));
  fund_leg->spread = (double *)calloc(fund_leg->num_cpn, sizeof(double));
  if ((!fund_leg->cpn) || (!fund_leg->notional) || (!fund_leg->margin) ||
      (!fund_leg->spread)) {
    err = "Allocation error in am_fill_fund_leg";
    goto FREE_RETURN;
  }

  /*	Fill coupons information */
  j = 0;
  while (i < fund_ncpn) {
    cpn = fund_leg->cpn + j;

    fund_leg->notional[j] = fund_not[i];
    fund_leg->margin[j] = fund_mrg[i];
    fund_leg->spread[j] = fund_spr[i];

    /*	Dates */
    cpn->start_date = fund_start[i];
    cpn->end_date = fund_end[i];
    cpn->pay_date = fund_pay[i];

    /*	Times */
    cpn->start_time = (cpn->start_date - today) * YEARS_IN_DAY;
    cpn->end_time = (cpn->end_date - today) * YEARS_IN_DAY;
    cpn->pay_time = (cpn->pay_date - today) * YEARS_IN_DAY;

    /*	Coupon */
    err = interp_basis(fund_basis[i], &bas);
    if (err) {
      goto FREE_RETURN;
    }
    cpn->cvg = coverage(fund_start[i], fund_end[i], bas);
    cpn->cpn = fund_not[i] * cpn->cvg * (fund_spr[i] + fund_mrg[i]);

    if (i < fund_ncpn - 1) {
      cpn->cpn += fund_not[i + 1] - fund_not[i];
    } else {
      cpn->cpn -= fund_not[i];
    }

    i++;
    j++;
  }

  err = am_check_fund_leg(fund_leg);

FREE_RETURN:

  if (err) {
    am_free_fund_leg(fund_leg);
  }

  return err;
}

/*	Check dates consistency */
Err am_check_fund_leg(AM_FUND_LEG fund_leg) {
  int i;

  /*	Check that start and pay dates are increasing */
  for (i = 1; i < fund_leg->num_cpn; i++) {
    if (fund_leg->cpn[i].start_date < fund_leg->cpn[i - 1].start_date) {
      return "Start dates should be increasing in funding leg";
    }

    if (fund_leg->cpn[i].pay_date < fund_leg->cpn[i - 1].pay_date) {
      return "Pay dates should be increasing in funding leg";
    }
  }

  /*	Check that pay dates are after start dates */
  for (i = 0; i < fund_leg->num_cpn; i++) {
    if (fund_leg->cpn[i].pay_date < fund_leg->cpn[i].start_date) {
      return "Pay dates should be after start dates in funding leg";
    }
  }

  /*	OK */
  return NULL;
}

/*	Free */
Err am_free_fund_leg(AM_FUND_LEG fund_leg) {
  if (fund_leg->notional) {
    free(fund_leg->notional);
    fund_leg->notional = NULL;
  }

  if (fund_leg->margin) {
    free(fund_leg->margin);
    fund_leg->margin = NULL;
  }

  if (fund_leg->spread) {
    free(fund_leg->spread);
    fund_leg->spread = NULL;
  }

  if (fund_leg->cpn) {
    free(fund_leg->cpn);
    fund_leg->cpn = NULL;
  }

  return NULL;
}

#define ONE_MONTH 0.083333333
/*	Functions for the fixtic leg */
Err am_fill_fix_leg(
    /*	Coupons that fixed before today are disregarded */
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I  , 1: E */
    double *fix_not, int fix_ncpn,
    // long		*fix_fix  ,
    long *fix_start, long *fix_end, long *fix_pay, char **fix_basis,
    double *fix_rate, double *fix_fee, AM_FIX_LEG fix_leg) {
  int i, j;
  SrtBasisCode bas;
  AM_FIX_CPN cpn;
  Err err = NULL;

  /*	Initialise pointers to NULL */
  fix_leg->notional = NULL;
  fix_leg->fee = NULL;
  fix_leg->cpn = NULL;
  fix_leg->rate = NULL;

  /*	Skip coupons fixed before today */
  i = 0;
  while (i < fix_ncpn && fix_start[i] < today + eod_flag) {
    i++;
  }

  /*	Check that at least one coupon is left */
  if (i == fix_ncpn) {
    /* err = "All funding coupons start before today in am_fill_fix_leg"; */
    fix_leg->num_cpn = 0;
    goto FREE_RETURN;
  }

  /*	Allocate memory */
  fix_leg->num_cpn = fix_ncpn - i;
  fix_leg->cpn = (am_fix_cpn *)calloc(fix_leg->num_cpn, sizeof(am_fix_cpn));
  fix_leg->notional = (double *)calloc(fix_leg->num_cpn, sizeof(double));
  fix_leg->rate = (double *)calloc(fix_leg->num_cpn, sizeof(double));
  fix_leg->fee = (double *)calloc(fix_leg->num_cpn, sizeof(double));
  if ((!fix_leg->cpn) || (!fix_leg->notional) || (!fix_leg->fee)) {
    err = "Allocation error in am_fill_fix_leg";
    goto FREE_RETURN;
  }

  /*	Fill coupons information */
  j = 0;
  while (i < fix_ncpn) {
    cpn = fix_leg->cpn + j;

    /*	Dates */
    cpn->start_date = fix_start[i];
    cpn->end_date = fix_end[i];
    cpn->pay_date = fix_pay[i];

    /*	Times */
    cpn->start_time = (cpn->start_date - today) * YEARS_IN_DAY;
    cpn->end_time = (cpn->end_date - today) * YEARS_IN_DAY;
    cpn->pay_time = (cpn->pay_date - today) * YEARS_IN_DAY;

    /*	Forward Libor and coverage */
    err = interp_basis(fix_basis[i], &bas);
    if (err) {
      goto FREE_RETURN;
    }

    fix_leg->notional[j] = fix_not[i];
    fix_leg->rate[j] = fix_rate[i];
    fix_leg->fee[j] = fix_fee[i];

    cpn->cvg = coverage(cpn->start_date, cpn->end_date, bas);
    cpn->cpn = fix_not[i] * fix_rate[i] * cpn->cvg;

    i++;
    j++;
  }

  err = am_check_fix_leg(fix_leg);

FREE_RETURN:

  if (err) {
    am_free_fix_leg(fix_leg);
  }

  return err;
}

/*	Check dates consistency */
Err am_check_fix_leg(AM_FIX_LEG fix_leg) {
  int i;

  /*	Notional has to be different from 0 */
  //	if (fabs(fix_leg->notional[0]) != 0.0)
  //	{
  //		return "Fix coupon notional has to greater than 0";
  //	}

  /*	Check that start  , pay  , fix and val dates are increasing */
  for (i = 1; i < fix_leg->num_cpn; i++) {
    if (fix_leg->cpn[i].start_date < fix_leg->cpn[i - 1].start_date) {
      return "Start dates should be increasing in fix leg";
    }

    if (fix_leg->cpn[i].pay_date < fix_leg->cpn[i - 1].pay_date) {
      return "Pay dates should be increasing in fix leg";
    }
  }

  /*	Check that pay dates are after start dates and fix dates */
  for (i = 0; i < fix_leg->num_cpn; i++) {
    if (fix_leg->cpn[i].pay_date < fix_leg->cpn[i].start_date) {
      return "Pay dates should be after start dates in fix leg";
    }
  }

  /*	OK */
  return NULL;
}

/*	Free */
Err am_free_fix_leg(AM_FIX_LEG fix_leg) {
  if (fix_leg->notional) {
    free(fix_leg->notional);
    fix_leg->notional = NULL;
  }

  if (fix_leg->fee) {
    free(fix_leg->fee);
    fix_leg->fee = NULL;
  }

  if (fix_leg->rate) {
    free(fix_leg->rate);
    fix_leg->rate = NULL;
  }

  if (fix_leg->cpn) {
    free(fix_leg->cpn);
    fix_leg->cpn = NULL;
  }

  return NULL;
}

/*	Functions for the calls */

Err am_fill_calls(
    /*	Exercises before today are disregarded */
    long today,
    /*	EOD Flag */
    int eod_flag,           /*	0: I  , 1: E */
    int ncall, int pay_rec, /*	0: rec pd  , 1: pay pd */
    long *ex_date, long *set_date, double *fee, AM_STR am) {
  int i, j, k;

  AM_FUND_LEG fund_leg;
  AM_FIX_LEG fix_leg;
  AM_CALL call;
  Err err = NULL;

  /*	Initialise pointers */
  am->call = NULL;
  fund_leg = am->fund_leg;
  fix_leg = am->fix_leg;

  /*	Skip calls to be exercised before today */
  i = 0;
  while (i < ncall && ex_date[i] < today + eod_flag) {
    i++;
  }

  /*	Check that at least one call is left */
  if (i == ncall) {
    /* err = "All calls are to be exercised before today in am_fill_calls"; */
    am->num_calls = 0;
    goto FREE_RETURN;
  }

  /*	Allocate memory */
  am->num_calls = ncall - i;
  am->call = (am_call *)calloc(am->num_calls, sizeof(am_call));
  if (!am->num_calls) {
    err = "Allocation error in am_fill_calls";
    goto FREE_RETURN;
  }

  /*	Fill coupons information */
  j = 0;

  while (i < ncall) {
    call = am->call + j;

    /*	Dates */
    call->ex_date = ex_date[i];
    call->set_date = set_date[i];

    /*	Times */
    call->ex_time = (am->call[j].ex_date - today) * YEARS_IN_DAY;
    call->set_time = (am->call[j].set_date - today) * YEARS_IN_DAY;

    /*	Call on funding leg */
    /*	k = index of the first coupon to be called on funding leg  ,
                    i.e. first coupon with a start date >= ex date */

    k = 0;
    while (k < fund_leg->num_cpn &&
           fund_leg->cpn[k].start_date < call->ex_date) {
      k++;
    }
    if (k == fund_leg->num_cpn) {
      err = serror("Call number %d does not control any coupon in funding leg",
                   i);
      goto FREE_RETURN;
    }
    call->fund_idx = k;
    call->num_fund_cpn = fund_leg->num_cpn - k;

    /*	Call on fix leg */
    /*	k = index of the first coupon to be called on fixtic leg  ,
                    i.e. first coupon with a start date >= ex date */
    k = 0;
    while (k < fix_leg->num_cpn && fix_leg->cpn[k].start_date < call->ex_date) {
      k++;
    }
    if (k == fix_leg->num_cpn) {
      err = serror("Call number %d does not control any coupon in fix leg", i);
      goto FREE_RETURN;
    }
    call->fix_idx = k;
    call->num_fix_cpn = fix_leg->num_cpn - k;

    /*	Payer or receiver */
    call->pay_rec = pay_rec;

    /*	Fee */
    call->fee = fee[i];

    i++;
    j++;
  }

  err = am_check_calls(am);

FREE_RETURN:

  if (err) {
    am_free_calls(am);
  }

  return err;
}

/*	Check dates consistency */
Err am_check_calls(AM_STR am) {
  int i;

  /*	Check that ex and set dates are strictly increasing
                  Also check that funding and fix indices are strictly
     increasing  , i.e. there is no redundant calls */
  for (i = 1; i < am->num_calls; i++) {
    if (am->call[i].ex_date <= am->call[i - 1].ex_date) {
      return "Exercise dates should be increasing";
    }

    if (am->call[i].set_date <= am->call[i - 1].set_date) {
      return "Settlement dates should be increasing";
    }

    if (am->call[i].fund_idx < am->call[i - 1].fund_idx) {
      return "Number of funding coupons controlled by calls should be "
             "decreasing";
    }

    if (am->call[i].fix_idx < am->call[i - 1].fix_idx) {
      return "Number of fix coupons controlled by calls should be decreasing";
    }

    if (am->call[i].fund_idx <= am->call[i - 1].fund_idx &&
        am->call[i].fix_idx <= am->call[i - 1].fix_idx) {
      //			return serror ("Calls %d and %d -indexed after
      //today- are redundant"  , i-1  , i);
    }
  }

  /*	Check that set dates are after ex dates
                  Also check that the call date is before the start  , end and
     fixing dates of the coupons it controls */
  for (i = 0; i < am->num_calls; i++) {
    if (am->call[i].set_date < am->call[i].ex_date) {
      return "Settlement dates should be after exercise dates";
    }

    if (am->fund_leg->cpn[am->call[i].fund_idx].start_date <
            am->call[i].ex_date ||
        am->fund_leg->cpn[am->call[i].fund_idx].pay_date <
            am->call[i].ex_date) {
      return "A funding coupon starts before its exercise date";
    }

    if (am->fix_leg->cpn[am->call[i].fix_idx].start_date <
            am->call[i].ex_date ||
        am->fix_leg->cpn[am->call[i].fix_idx].pay_date < am->call[i].ex_date) {
      return "A fix coupon starts or fixes before its exercise date";
    }
  }

  /*	OK */
  return NULL;
}

/*	Free */
Err am_free_calls(AM_STR am) {
  if (am->call) {
    free(am->call);
    am->call = NULL;
  }

  return NULL;
}

/*	Functions for the underlying */

/*	Fill underlying structure from a predefined underlying */
Err am_fill_und(char *lgm2dund, char *vc, char *ref, char *swap_freq,
                char *swap_basis, AM_UND und) {
  double *sig_time = NULL, *sig = NULL;

  SrtUndPtr srtund;

  int i;
  Err err = NULL;
  double *lambda = NULL;
  double *lambda_time = NULL;
  int lambda_n;

  /*	Initialise */
  und->sigma_date = NULL;
  und->sigma_time = NULL;
  und->sigma = NULL;

  strcpy(und->name, lgm2dund);

  /*	Get term structures */
  err = Get_LGM2F_TermStructure2(
      lgm2dund, &(und->sigma), &(und->sigma_time), &(und->sigma_n), &(lambda),
      &(lambda_time), &(lambda_n), &(und->alpha), &(und->gamma), &(und->rho));
  if (err) {
    goto FREE_RETURN;
  }

  if (lambda_n > 1) {
    err = "Lambda term structure not allowed";
    goto FREE_RETURN;
  }

  und->lambda = lambda[0];

  /* write default values */
  und->lambda_n = lambda_n;
  und->pdlambda = 0;
  und->pdlambda_date = 0;
  und->pdlambda_time = 0;

  /*	Fill dates */
  srtund = lookup_und(lgm2dund);
  und->today = get_today_from_underlying(srtund);
  strcpy(und->yc, get_ycname_from_irund(srtund));
  strcpy(und->vc, vc);
  strcpy(und->ref, ref);
  strcpy(und->swap_freq, swap_freq);
  strcpy(und->swap_basis, swap_basis);

  und->sigma_date = (double *)calloc(und->sigma_n, sizeof(double));

  if (!und->sigma_date) {
    err = "Allocation error in am_fill_und";
    goto FREE_RETURN;
  }

  for (i = 0; i < und->sigma_n; i++) {
    und->sigma_date[i] =
        und->today + und->sigma_time[i] * DAYS_IN_YEAR + 1.0e-08;
  }

  und->has_inst_data = 0;
  cpd_init_calib_inst_data(&(und->inst_data));

FREE_RETURN:

  if (err) {
    am_free_und(und);
  }

  if (lambda)
    free(lambda);
  if (lambda_time)
    free(lambda_time);

  return err;
}

Err amortMidat_compute_diagonal_prices(
    char *yc, char *vc, char *ref, char *cFreq, char *cBasis, long today,
    int eod_flag, double mintime, double mininterval, int UseVol, int *nex,
    int *firstex, int **ex_bool, double **diag_prices, int model,
    Err (*get_correl)(char *vol_curve_name, double start_date, double end_date,
                      double strike, double *vol),
    char *CorrelName,
    Err (*get_cash_vol)(char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*
                            Err  (*GetVolForBadr)( Date  , Date  , double  ,
       SRT_Boolean  , double *)  , char *cVolType  ,
    */
    AM_STR am) {
  int i, j, k, l;
  Err err = NULL;

  long StartDate, EndDate, SpotDate, ExDate;

  int lNFixNot;
  double *dFixNotionals = NULL;
  double *dFixRates = NULL;
  double *dFee = NULL;
  long *lPayDates = NULL, *lStartDates = NULL, *lEndDates = NULL;

  long *lFixStartDates = NULL, *lFixEndDates = NULL;
  long *lFloatStartDates = NULL, *lFloatEndDates = NULL;

  double *dCoverages = NULL;

  int lNFloatNot;
  double *dFloatNotionals = NULL;
  double *dMargins = NULL;
  double *dSpreads = NULL;

  double **Correl = NULL;
  int NDimCorrel;

  double exer_fee;
  double dPrice;
  SrtCompounding srtFreq;
  SrtBasisCode srtBasis;

  // Output of Badr function
  double *dvPayDates = NULL;
  double *dvReplicatingStrikes = NULL;
  double *dvReplicatingNotionals = NULL;
  double *dvReplicatingSwaptions = NULL;
  double dFixedPV;
  double dFloatPV;
  double dSwapRate;
  SigKapTS *lgmSigKapTSPtrPtr = NULL;
  double dLGMVol;

  char *cRec = NULL;

  double *pdMargins = NULL;
  double *dFixNots = NULL;

  double lastcaltime;

  //	SrtDiffusionType srt_vol_type;

  cRec = "REC";

  err = interp_compounding(cFreq, &srtFreq);
  if (err) {
    smessage("Wrong Frequency");
    err = "Wrong Frequency";
    goto FREE_RETURN;
  }
  err = interp_basis(cBasis, &srtBasis);
  if (err) {
    smessage("Wrong Basis");
    err = "Wrong Basis";
    goto FREE_RETURN;
  }
  /*	err = interp_diffusion_type (cVolType  , &srt_vol_type);
          if(err)
          {
                  smessage("Wrong Diffusion Type");
                  err = "Wrong Diffusion Type";
                  goto FREE_RETURN;
          }
  */
  SpotDate = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

  lNFixNot = am->fix_leg->num_cpn;
  dFixNotionals = (double *)calloc(lNFixNot, sizeof(double));
  dFixRates = (double *)calloc(lNFixNot, sizeof(double));
  lFixStartDates = (long *)calloc(lNFixNot, sizeof(long));
  lFixEndDates = (long *)calloc(lNFixNot, sizeof(long));
  dFee = (double *)calloc(lNFixNot, sizeof(double));
  if ((!dFixNotionals) || (!dFixRates) || (!dFee)) {
    err = "Allocation failed in amortMidat_compute_diagonal_prices";
    smessage("Allocation failed in amortMidat_compute_diagonal_prices");
    goto FREE_RETURN;
  }

  // Inputs for Badr function
  lStartDates = (long *)calloc(lNFixNot + 1, sizeof(long));
  lEndDates = (long *)calloc(lNFixNot + 1, sizeof(long));
  lPayDates = (long *)calloc(lNFixNot + 1, sizeof(long));
  dCoverages = (double *)calloc(lNFixNot + 1, sizeof(double));
  dFixNots = (double *)calloc(lNFixNot + 1, sizeof(double));
  pdMargins = (double *)calloc(lNFixNot + 1, sizeof(double));

  // Outputs of Badr function
  dvPayDates = (double *)calloc(lNFixNot + 1, sizeof(double));
  dvReplicatingStrikes = (double *)calloc(lNFixNot + 1, sizeof(double));
  dvReplicatingNotionals = (double *)calloc(lNFixNot + 1, sizeof(double));
  dvReplicatingSwaptions = (double *)calloc(lNFixNot + 1, sizeof(double));
  if ((!dvPayDates) || (!dvReplicatingStrikes) || (!dvReplicatingNotionals) ||
      (!dvReplicatingSwaptions) || (!pdMargins) || (!lStartDates) ||
      (!lEndDates) || (!lPayDates) || (!dCoverages)) {
    err = "Allocation failed in amortMidat_compute_diagonal_prices";
    smessage("Allocation failed in amortMidat_compute_diagonal_prices");
    goto FREE_RETURN;
  }

  pdMargins[0] = 0;
  lEndDates[0] = 0;
  lStartDates[0] = 0;
  lStartDates[0] = 0;
  dFixNots[0] = 0;
  j = 0;
  for (i = 0; i < lNFixNot; ++i) {
    dFixNotionals[i] = am->fix_leg->notional[i];
    lFixStartDates[i] = am->fix_leg->cpn[i].start_date;
    lFixEndDates[i] = am->fix_leg->cpn[i].pay_date;

    dFixRates[i] = am->fix_leg->rate[i];
    dFee[i] = am->fix_leg->fee[i];

    // Badr function
    lStartDates[i + 1] = am->fix_leg->cpn[i].start_date;
    lEndDates[i + 1] = am->fix_leg->cpn[i].pay_date;
    lPayDates[i + 1] = am->fix_leg->cpn[i].pay_date;
    dCoverages[i + 1] = am->fix_leg->cpn[i].cvg;
    dFixNots[i + 1] = am->fix_leg->notional[i];
    while ((j < am->fund_leg->num_cpn - 1) &&
           (am->fund_leg->cpn[j].start_date <
            am->fix_leg->cpn[i].start_date - 10)) {
      ++j;
    }
    j = j - 1;
    pdMargins[i + 1] = am->fund_leg->margin[j];
  }

  lNFloatNot = am->fund_leg->num_cpn;

  lFloatStartDates = (long *)calloc(lNFloatNot, sizeof(long));
  lFloatEndDates = (long *)calloc(lNFloatNot, sizeof(long));
  dFloatNotionals = (double *)calloc(lNFloatNot, sizeof(double));

  dMargins = (double *)calloc(lNFloatNot, sizeof(double));
  dSpreads = (double *)calloc(lNFloatNot, sizeof(double));
  if ((!dFloatNotionals) || (!dMargins)) {
    err = "Allocation failed in amortMidat_compute_diagonal_prices";
    smessage("Allocation failed in amortMidat_compute_diagonal_prices");
    goto FREE_RETURN;
  }
  for (i = 0; i < lNFloatNot; ++i) {
    dFloatNotionals[i] = am->fund_leg->notional[i];

    lFloatStartDates[i] = am->fund_leg->cpn[i].start_date;
    lFloatEndDates[i] = am->fund_leg->cpn[i].pay_date;

    dMargins[i] = am->fund_leg->margin[i];
    dSpreads[i] = am->fund_leg->spread[i];
  }

  /*	Skip calls to be exercised before today */
  i = 0;
  j = 0;

  l = 0;

  *nex = am->fix_leg->num_cpn - j;
  *firstex = j;

  *ex_bool = (int *)calloc(*nex, sizeof(int));
  *diag_prices = (double *)calloc(*nex, sizeof(double));
  if ((!(*ex_bool)) || (!(*diag_prices))) {
    err = "allocation failed in amortMidat_compute_diagonal_prices";
    smessage("allocation failed in amortMidat_compute_diagonal_prices");
    goto FREE_RETURN;
  }

  EndDate = am->theoEndDate;

  lastcaltime = 0.0;
  if (am->fix_leg->cpn[j].start_time < DMAX(mintime, am->call[0].ex_time)) {
    (*ex_bool)[0] = 0;
    (*diag_prices)[0] = 0;
  } else {
    (*ex_bool)[0] = 1;
    exer_fee = dFee[j];
    StartDate = am->fix_leg->cpn[j].start_date;
    lastcaltime = am->fix_leg->cpn[j].start_time;

    while ((l < am->fund_leg->num_cpn) &&
           (am->fund_leg->cpn[l].start_time <
            am->fix_leg->cpn[j].start_time - 10.0 / 365.0)) {
      ++l;
    }

    if (model == 1) {
      err = Compute_CoInitalSwaps_Correl(SpotDate, CorrelName, get_correl,
                                         StartDate, EndDate, srtFreq, srtBasis,
                                         &NDimCorrel, &Correl);
      if (err) {
        goto FREE_RETURN;
      }

      err = AmortizedSwaptionShiftedLog(
          yc, vc, ref,

          StartDate, EndDate,

          srtFreq, srtBasis,

          exer_fee,

          lNFixNot - j, dFixNotionals + j, dFixRates + j,

          lNFloatNot - l, dFloatNotionals + l, dMargins + l,

          SRT_PUT,

          0, 10000, 1.0, 10,

          UseVol,

          Correl, &dPrice);
      if (err) {
        smessage("Pb in pricing European with Shifted Log SMM");
        goto FREE_RETURN;
      }

      if (Correl) {
        free_dmatrix(Correl, 0, NDimCorrel - 1, 0, NDimCorrel - 1);
        Correl = NULL;
      }
    } else {
      ExDate = add_unit(StartDate, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
      err = EuropeanAmortizingSwaption2(
          today, ExDate, EndDate, lNFixNot - j, lStartDates + j, lEndDates + j,
          lPayDates + j, dCoverages + j, dCoverages + j, dFixRates[0],
          dFixNots + j, yc, vc, cFreq, cBasis, cRec, ref, get_cash_vol,
          /*
                                                                          GetVolForBadr
             , srt_vol_type  ,
          */
          // Output
          dvPayDates, dvReplicatingStrikes, dvReplicatingNotionals, &dPrice,
          &dFixedPV, &dFloatPV, &dSwapRate, &lgmSigKapTSPtrPtr, &dLGMVol,
          dvReplicatingSwaptions,

          // Inputs
          pdMargins + j,
          0, // dFixRates+j  ,
          0, exer_fee);
      if (err) {
        smessage("Pb in pricing European with 1F");
        goto FREE_RETURN;
      }
    }

    (*diag_prices)[0] = dPrice;
  }

  for (k = 1; k < *nex; ++k) {
    if ((am->fix_leg->cpn[j + k].start_time - lastcaltime < mininterval) ||
        (am->fix_leg->cpn[j + k].start_time <
         DMAX(mintime, am->call[0].ex_time))) {
      (*ex_bool)[k] = 0;
      (*diag_prices)[k] = 0;
    } else {

      StartDate = am->fix_leg->cpn[j + k].start_date;
      lastcaltime = am->fix_leg->cpn[j + k].start_time;

      exer_fee = dFee[j + k];

      (*ex_bool)[k] = 1;

      while ((l < am->fund_leg->num_cpn) &&
             (am->fund_leg->cpn[l].start_time <
              am->fix_leg->cpn[j + k].start_time - 10.0 / 365.0)) {
        ++l;
      }

      if (model == 1) {
        err = Compute_CoInitalSwaps_Correl(SpotDate, CorrelName, get_correl,
                                           StartDate, EndDate, srtFreq,
                                           srtBasis, &NDimCorrel, &Correl);
        if (err) {
          goto FREE_RETURN;
        }

        err = AmortizedSwaptionShiftedLog(
            yc, vc, ref,

            StartDate, EndDate,

            srtFreq, srtBasis,

            exer_fee,

            lNFixNot - (k + j), dFixNotionals + k + j, dFixRates + k + j,

            lNFloatNot - l, dFloatNotionals + l, dMargins + l,

            SRT_PUT,

            0, 10000, 1.0, 10,

            UseVol,

            Correl, &dPrice);
        if (err) {
          smessage("Pb in pricing European with Shifted Log SMM");
          goto FREE_RETURN;
        }

        if (Correl) {
          free_dmatrix(Correl, 0, NDimCorrel - 1, 0, NDimCorrel - 1);
          Correl = NULL;
        }
      } else {
        ExDate = add_unit(StartDate, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
        err = EuropeanAmortizingSwaption2(
            today, ExDate, EndDate, lNFixNot - (k + j), lStartDates + (k + j),
            lEndDates + (k + j), lPayDates + (k + j), dCoverages + (k + j),
            dCoverages + (k + j), dFixRates[0], dFixNots + (k + j), yc, vc,
            cFreq, cBasis, cRec, ref, get_cash_vol,
            /*
                                                                            GetVolForBadr
               , srt_vol_type  ,
            */
            // Output
            dvPayDates, dvReplicatingStrikes, dvReplicatingNotionals, &dPrice,
            &dFixedPV, &dFloatPV, &dSwapRate, &lgmSigKapTSPtrPtr, &dLGMVol,
            dvReplicatingSwaptions,

            // Inputs
            pdMargins + k + j,
            0, // dFixRates+k+j  ,
            0, exer_fee);
        if (err) {
          smessage("Pb in pricing European with 1F");
          goto FREE_RETURN;
        }
      }

      (*diag_prices)[k] = dPrice;
    }
  }

FREE_RETURN:

  if (lStartDates)
    free(lStartDates);
  if (lEndDates)
    free(lEndDates);
  if (lPayDates)
    free(lPayDates);
  if (dCoverages)
    free(dCoverages);
  if (dFixNots)
    free(dFixNots);
  if (pdMargins)
    free(pdMargins);

  if (dvPayDates)
    free(dvPayDates);
  if (dvReplicatingStrikes)
    free(dvReplicatingStrikes);
  if (dvReplicatingNotionals)
    free(dvReplicatingNotionals);
  if (dvReplicatingSwaptions)
    free(dvReplicatingSwaptions);

  if (dFixNotionals)
    free(dFixNotionals);
  if (lFixStartDates)
    free(lFixStartDates);
  if (lFixEndDates)
    free(lFixEndDates);
  if (dFixRates)
    free(dFixRates);
  if (dFee)
    free(dFee);
  if (dFloatNotionals)
    free(dFloatNotionals);
  if (lFloatStartDates)
    free(lFloatStartDates);
  if (lFloatEndDates)
    free(lFloatEndDates);
  if (dMargins)
    free(dMargins);
  if (dSpreads)
    free(dSpreads);

  if (Correl) {
    free_dmatrix(Correl, 0, NDimCorrel - 1, 0, NDimCorrel - 1);
    Correl = NULL;
  }

  if (err) {
    if (*ex_bool) {
      free(*ex_bool);
      *ex_bool = NULL;
    }
    if (*diag_prices) {
      free(*diag_prices);
      *diag_prices = NULL;
    }
  }

  return err;
}

Err amortMidat_compute_diagonal_prices_new(
    char *yc, char *vc, char *ref, char *cFreq, char *cBasis, long today,
    int eod_flag, double mintime, double mininterval, int UseVol, int *nex,
    int *firstex, int **ex_bool, double **diag_prices, int model,
    Err (*get_correl)(char *vol_curve_name, double start_date, double end_date,
                      double strike, double *vol),
    char *CorrelName,
    Err (*get_cash_vol)(char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),

    AM_STR am) {
  int i, j, k, l;
  Err err = NULL;

  long StartDate, EndDate, SpotDate, ExDate;

  int lNFixNot;
  double *dFixNotionals = NULL;
  double *dFixRates = NULL;
  double *dFee = NULL;
  long *lPayDates = NULL, *lStartDates = NULL, *lEndDates = NULL;

  long *lFixStartDates = NULL, *lFixEndDates = NULL;
  long *lFloatStartDates = NULL, *lFloatEndDates = NULL;

  double *dCoverages = NULL;

  int lNFloatNot;
  double *dFloatNotionals = NULL;
  double *dMargins = NULL;
  double *dSpreads = NULL;

  double **Correl = NULL;
  int NDimCorrel;

  double exer_fee;
  double dPrice;
  SrtCompounding srtFreq;
  SrtBasisCode srtBasis;

  // Output of Badr function
  double *dvPayDates = NULL;
  double *dvReplicatingStrikes = NULL;
  double *dvReplicatingNotionals = NULL;
  double *dvReplicatingSwaptions = NULL;
  double dFixedPV;
  double dFloatPV;
  double dSwapRate;
  SigKapTS *lgmSigKapTSPtrPtr = NULL;
  double dLGMVol;

  char *cRec = NULL;

  double *pdMargins = NULL;
  double *dFixNots = NULL;

  double lastcaltime;

  /*	SrtDiffusionType srt_vol_type; */

  cRec = "REC";

  err = interp_compounding(cFreq, &srtFreq);
  if (err) {
    smessage("Wrong Frequency");
    err = "Wrong Frequency";
    goto FREE_RETURN;
  }
  err = interp_basis(cBasis, &srtBasis);
  if (err) {
    smessage("Wrong Basis");
    err = "Wrong Basis";
    goto FREE_RETURN;
  }

  SpotDate = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

  lNFixNot = am->fix_leg->num_cpn;
  dFixNotionals = (double *)calloc(lNFixNot, sizeof(double));
  dFixRates = (double *)calloc(lNFixNot, sizeof(double));
  lFixStartDates = (long *)calloc(lNFixNot, sizeof(long));
  lFixEndDates = (long *)calloc(lNFixNot, sizeof(long));
  dFee = (double *)calloc(lNFixNot, sizeof(double));
  if ((!dFixNotionals) || (!dFixRates) || (!dFee)) {
    err = "Allocation failed in amortMidat_compute_diagonal_prices";
    smessage("Allocation failed in amortMidat_compute_diagonal_prices");
    goto FREE_RETURN;
  }

  // Inputs for Badr function
  lStartDates = (long *)calloc(lNFixNot + 1, sizeof(long));
  lEndDates = (long *)calloc(lNFixNot + 1, sizeof(long));
  lPayDates = (long *)calloc(lNFixNot + 1, sizeof(long));
  dCoverages = (double *)calloc(lNFixNot + 1, sizeof(double));
  dFixNots = (double *)calloc(lNFixNot + 1, sizeof(double));
  pdMargins = (double *)calloc(lNFixNot + 1, sizeof(double));

  // Outputs of Badr function
  dvPayDates = (double *)calloc(lNFixNot + 1, sizeof(double));
  dvReplicatingStrikes = (double *)calloc(lNFixNot + 1, sizeof(double));
  dvReplicatingNotionals = (double *)calloc(lNFixNot + 1, sizeof(double));
  dvReplicatingSwaptions = (double *)calloc(lNFixNot + 1, sizeof(double));
  if ((!dvPayDates) || (!dvReplicatingStrikes) || (!dvReplicatingNotionals) ||
      (!dvReplicatingSwaptions) || (!pdMargins) || (!lStartDates) ||
      (!lEndDates) || (!lPayDates) || (!dCoverages)) {
    err = "Allocation failed in amortMidat_compute_diagonal_prices";
    smessage("Allocation failed in amortMidat_compute_diagonal_prices");
    goto FREE_RETURN;
  }

  pdMargins[0] = 0;
  lEndDates[0] = 0;
  lStartDates[0] = 0;
  lStartDates[0] = 0;
  dFixNots[0] = 0;
  j = 0;
  for (i = 0; i < lNFixNot; ++i) {
    dFixNotionals[i] = am->fix_leg->notional[i];
    lFixStartDates[i] = am->fix_leg->cpn[i].start_date;
    lFixEndDates[i] = am->fix_leg->cpn[i].end_date;

    dFixRates[i] = am->fix_leg->rate[i];
    dFee[i] = am->fix_leg->fee[i];

    // Badr function
    lStartDates[i + 1] = am->fix_leg->cpn[i].start_date;
    lEndDates[i + 1] = am->fix_leg->cpn[i].end_date;
    lPayDates[i + 1] = am->fix_leg->cpn[i].end_date;
    dCoverages[i + 1] = am->fix_leg->cpn[i].cvg;
    dFixNots[i + 1] = am->fix_leg->notional[i];
    while ((j < am->fund_leg->num_cpn - 1) &&
           (am->fund_leg->cpn[j].start_date <
            am->fix_leg->cpn[i].start_date - 10)) {
      ++j;
    }
    pdMargins[i + 1] = am->fund_leg->margin[j];
  }

  lNFloatNot = am->fund_leg->num_cpn;

  lFloatStartDates = (long *)calloc(lNFloatNot, sizeof(long));
  lFloatEndDates = (long *)calloc(lNFloatNot, sizeof(long));
  dFloatNotionals = (double *)calloc(lNFloatNot, sizeof(double));

  dMargins = (double *)calloc(lNFloatNot, sizeof(double));
  dSpreads = (double *)calloc(lNFloatNot, sizeof(double));
  if ((!dFloatNotionals) || (!dMargins)) {
    err = "Allocation failed in amortMidat_compute_diagonal_prices";
    smessage("Allocation failed in amortMidat_compute_diagonal_prices");
    goto FREE_RETURN;
  }
  for (i = 0; i < lNFloatNot; ++i) {
    dFloatNotionals[i] = am->fund_leg->notional[i];

    lFloatStartDates[i] = am->fund_leg->cpn[i].start_date;
    lFloatEndDates[i] = am->fund_leg->cpn[i].end_date;

    dMargins[i] = am->fund_leg->margin[i];
    dSpreads[i] = am->fund_leg->spread[i];
  }

  /*	Skip calls to be exercised before today */
  i = 0;
  j = 0;

  l = 0;

  *nex = am->fix_leg->num_cpn - j;
  *firstex = j;

  *ex_bool = (int *)calloc(*nex, sizeof(int));
  *diag_prices = (double *)calloc(*nex, sizeof(double));
  if ((!(*ex_bool)) || (!(*diag_prices))) {
    err = "allocation failed in amortMidat_compute_diagonal_prices";
    smessage("allocation failed in amortMidat_compute_diagonal_prices");
    goto FREE_RETURN;
  }

  EndDate = am->theoEndDate;

  lastcaltime = 0.0;
  if (am->fix_leg->cpn[j].start_time < DMAX(mintime, am->call[0].ex_time)) {
    (*ex_bool)[0] = 0;
    (*diag_prices)[0] = 0;
  } else {
    (*ex_bool)[0] = 1;
    exer_fee = dFee[j];
    StartDate = am->fix_leg->cpn[j].start_date;
    lastcaltime = am->fix_leg->cpn[j].start_time;

    while ((l < am->fund_leg->num_cpn) &&
           (am->fund_leg->cpn[l].start_time <
            am->fix_leg->cpn[j].start_time - 10.0 / 365.0)) {
      ++l;
    }

    if (model == 1) {
      NDimCorrel = lNFixNot - j;
      err = Compute_CoInitalSwaps_Correl2(SpotDate, CorrelName, get_correl,
                                          StartDate, EndDate, srtFreq, srtBasis,
                                          NDimCorrel, &Correl);
      if (err) {
        goto FREE_RETURN;
      }

      err = AmortizedSwaptionShiftedLogForMAD(
          yc, vc, ref,

          srtFreq, srtBasis,

          exer_fee,

          lFixStartDates + j, lFixEndDates + j, lNFixNot - j, dFixNotionals + j,
          dFixRates + j,

          lFloatStartDates + l, lFloatEndDates + l, lNFloatNot - l,
          dFloatNotionals + l, dMargins + l, dSpreads + l,

          SRT_PUT,

          10000,

          0, 0,

          UseVol,

          Correl, &dPrice);
      if (err) {
        smessage("Pb in pricing European with Shifted Log SMM");
        goto FREE_RETURN;
      }

      if (Correl) {
        free_dmatrix(Correl, 0, NDimCorrel - 1, 0, NDimCorrel - 1);
        Correl = NULL;
      }
    } else {
      ExDate = add_unit(StartDate, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
      err = EuropeanAmortizingSwaption2(
          today, ExDate, EndDate, lNFixNot - j, lStartDates + j, lEndDates + j,
          lPayDates + j, dCoverages + j, dCoverages + j, dFixRates[0],
          dFixNots + j, yc, vc, cFreq, cBasis, cRec, ref, get_cash_vol,

          // Output
          dvPayDates, dvReplicatingStrikes, dvReplicatingNotionals, &dPrice,
          &dFixedPV, &dFloatPV, &dSwapRate, &lgmSigKapTSPtrPtr, &dLGMVol,
          dvReplicatingSwaptions,

          // Inputs
          pdMargins + j,
          0, // dFixRates+j  ,
          0, exer_fee);
      if (err) {
        smessage("Pb in pricing European with 1F");
        goto FREE_RETURN;
      }
    }

    (*diag_prices)[0] = dPrice;
  }

  for (k = 1; k < *nex; ++k) {
    if ((am->fix_leg->cpn[j + k].start_time - lastcaltime < mininterval) ||
        (am->fix_leg->cpn[j + k].start_time <
         DMAX(mintime, am->call[0].ex_time))) {
      (*ex_bool)[k] = 0;
      (*diag_prices)[k] = 0;
    } else {

      StartDate = am->fix_leg->cpn[j + k].start_date;
      lastcaltime = am->fix_leg->cpn[j + k].start_time;

      exer_fee = dFee[j + k];

      (*ex_bool)[k] = 1;

      while ((l < am->fund_leg->num_cpn) &&
             (am->fund_leg->cpn[l].start_time <
              am->fix_leg->cpn[j + k].start_time - 10.0 / 365.0)) {
        ++l;
      }

      if (model == 1) {
        NDimCorrel = lNFixNot - (k + j);
        err = Compute_CoInitalSwaps_Correl2(SpotDate, CorrelName, get_correl,
                                            StartDate, EndDate, srtFreq,
                                            srtBasis, NDimCorrel, &Correl);
        if (err) {
          goto FREE_RETURN;
        }

        err = AmortizedSwaptionShiftedLogForMAD(
            yc, vc, ref,

            srtFreq, srtBasis,

            exer_fee,

            lFixStartDates + k + j, lFixEndDates + k + j, lNFixNot - (k + j),
            dFixNotionals + k + j, dFixRates + k + j,

            lFloatStartDates + l, lFloatEndDates + l, lNFloatNot - l,
            dFloatNotionals + l, dMargins + l, dSpreads + l,

            SRT_PUT,

            10000,

            0, 0,

            UseVol,

            Correl, &dPrice);
        if (err) {
          smessage("Pb in pricing European with Shifted Log SMM");
          goto FREE_RETURN;
        }

        if (Correl) {
          free_dmatrix(Correl, 0, NDimCorrel - 1, 0, NDimCorrel - 1);
          Correl = NULL;
        }
      } else {
        ExDate = add_unit(StartDate, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
        err = EuropeanAmortizingSwaption2(
            today, ExDate, EndDate, lNFixNot - (k + j), lStartDates + (k + j),
            lEndDates + (k + j), lPayDates + (k + j), dCoverages + (k + j),
            dCoverages + (k + j), dFixRates[0], dFixNots + (k + j), yc, vc,
            cFreq, cBasis, cRec, ref, get_cash_vol,

            // Output
            dvPayDates, dvReplicatingStrikes, dvReplicatingNotionals, &dPrice,
            &dFixedPV, &dFloatPV, &dSwapRate, &lgmSigKapTSPtrPtr, &dLGMVol,
            dvReplicatingSwaptions,

            // Inputs
            pdMargins + k + j,
            0, // dFixRates+k+j  ,
            0, exer_fee);
        if (err) {
          smessage("Pb in pricing European with 1F");
          goto FREE_RETURN;
        }
      }

      (*diag_prices)[k] = dPrice;
    }
  }

FREE_RETURN:

  if (lStartDates)
    free(lStartDates);
  if (lEndDates)
    free(lEndDates);
  if (lPayDates)
    free(lPayDates);
  if (dCoverages)
    free(dCoverages);
  if (dFixNots)
    free(dFixNots);
  if (pdMargins)
    free(pdMargins);

  if (dvPayDates)
    free(dvPayDates);
  if (dvReplicatingStrikes)
    free(dvReplicatingStrikes);
  if (dvReplicatingNotionals)
    free(dvReplicatingNotionals);
  if (dvReplicatingSwaptions)
    free(dvReplicatingSwaptions);

  if (dFixNotionals)
    free(dFixNotionals);
  if (lFixStartDates)
    free(lFixStartDates);
  if (lFixEndDates)
    free(lFixEndDates);
  if (dFixRates)
    free(dFixRates);
  if (dFee)
    free(dFee);
  if (dFloatNotionals)
    free(dFloatNotionals);
  if (lFloatStartDates)
    free(lFloatStartDates);
  if (lFloatEndDates)
    free(lFloatEndDates);
  if (dMargins)
    free(dMargins);
  if (dSpreads)
    free(dSpreads);

  if (Correl) {
    free_dmatrix(Correl, 0, NDimCorrel - 1, 0, NDimCorrel - 1);
    Correl = NULL;
  }

  if (err) {
    if (*ex_bool) {
      free(*ex_bool);
      *ex_bool = NULL;
    }
    if (*diag_prices) {
      free(*diag_prices);
      *diag_prices = NULL;
    }
  }

  return err;
}

/*	Fill underlying structure from calibration instruments */
Err am_calib_und(
    long today,
    /*	EOD Flag */
    int eod_flag,     /*	0: I  , 1: E */
    char *yc,         /*	yc */
    char *vc,         /*	vc (only if calib) */
    char *ref,        /*	ref rate (only if calib) */
    char *swap_freq,  /*	swap freq (only if calib) */
    char *swap_basis, /*	swap basis (only if calib) */
    double lambda,    /*	lambda if unique */
    double alpha,     /*	alpha */
    double gamma,     /*	gamma */
    double rho,       /*	rho */
    /*	Calib params */
    double mintime, double mininterval,

    int notperiod, double max_std_short, int one2F,
    int fix_lambda, /*	0: calib lambda to cap  , 1: fix lambda calib
                                                    to diagonal */
    int one_f_equi, /*	1F equivalent flag:
                                                    if set to 1  , then 2F
                       lambda will calibrate to the cap priced within calibrated
                       1F with the given lambda */
    int skip_last,  /*	If 1  , the last option is disregarded
                                                    and the forward volatility is
                       flat from option  n-1 */
    double max_var_jump,

    int strike_type, int european_model,

    Err (*get_correl)(char *correl_cube_name, double start_date,
                      double end_date, double strike, double *vol),
    char *CorrelName,

    Err (*GetVolForBadr)(Date, Date, double, SRT_Boolean, double *),
    char *cVolType,

    /*	End of calib params */
    AM_STR am,          /*	structure */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    AM_UND und) {
  int i, nex, firstex;
  long *lex = NULL;

  long end_struct; //  ,
  //				end_swap;
  double *long_strikes = NULL, *short_strikes = NULL;

  int Nfix;
  int *ex_bool = NULL;
  double *diag_prices = NULL;
  double *ex_fee = NULL;
  double *fixNotional = NULL;
  long *fix_start_dates = NULL;
  long *fix_end_dates = NULL;

  int Nfloat;
  double *floatNotional = NULL;
  long *float_start_dates = NULL;
  long *float_end_dates = NULL;
  double *float_margin = NULL;
  double *float_spread = NULL;

  long start_date;
  long end_date;

  int UseVol;

  AM_FUND_LEG fund_leg;
  AM_FIX_LEG fix_leg;

  Err err = NULL;

  /*	Initialise */
  fix_leg = am->fix_leg;
  fund_leg = am->fund_leg;

  und->sigma_date = NULL;
  und->sigma_time = NULL;
  und->sigma = NULL;

  und->has_inst_data = 0;
  cpd_init_calib_inst_data(&(und->inst_data));

  und->num_prices = 0;
  und->exercise_dates = NULL;
  und->mkt_prices = NULL;
  und->mdl_prices = NULL;

  und->lambda = lambda;

  /* write default values */
  und->lambda_n = 1;
  und->pdlambda = 0;
  und->pdlambda_date = 0;
  und->pdlambda_time = 0;

  und->alpha = alpha;
  und->gamma = gamma;
  und->rho = rho;
  und->today = today;

  strcpy(und->yc, yc);
  strcpy(und->vc, vc);
  strcpy(und->ref, ref);
  strcpy(und->swap_freq, swap_freq);
  strcpy(und->swap_basis, swap_basis);
  strcpy(und->name, "CALIB");

  /*	Exercise dates for calibration */
  end_struct = fund_leg->cpn[fund_leg->num_cpn - 1].pay_date;
  //	end_swap = end_struct;
  //	if (fix_leg->cpn[fix_leg->num_cpn-1].pay_date > end_swap)
  //	{
  //		end_swap = fix_leg->cpn[fix_leg->num_cpn-1].pay_date;
  //	}

  if (am->num_calls > 0 &&
      !(am->num_calls == 1 && am->call[0].ex_date <= today + eod_flag)) {
    /*	If call dates  , choose call dates as option expiries for calibration */
    nex = am->num_calls;
    lex = (long *)calloc(nex, sizeof(long));
    if (!lex) {
      err = "Allocation error (2) in am_calib_und";
      goto FREE_RETURN;
    }

    for (i = 0; i < nex; i++) {
      lex[i] = am->call[i].ex_date;
    }
  } else {
    /*	If no call dates  , exit */
    am->num_calls = 0;
    goto FREE_RETURN;
  }

  UseVol = 0;
  if (strike_type == 5) {
    UseVol = 1;
  }

  err = amortMidat_compute_diagonal_prices(
      yc, vc, ref, swap_freq, swap_basis, today, eod_flag, mintime, mininterval,
      UseVol, &nex, &firstex, &ex_bool, &diag_prices, european_model,
      get_correl, CorrelName, get_cash_vol,
      /*GetVolForBadr  , cVolType  ,*/ am);
  if (err) {
    goto FREE_RETURN;
  }

  und->num_prices = nex;
  und->exercise_dates = (double *)calloc(nex, sizeof(double));
  und->mkt_prices = (double *)calloc(nex, sizeof(double));
  und->mdl_prices = (double *)calloc(nex, sizeof(double));
  for (i = 0; i < nex; ++i) {
    und->exercise_dates[i] = am->fix_leg->cpn[firstex + i].start_date;
    und->mkt_prices[i] = diag_prices[i];
  }

  Nfix = am->fix_leg->num_cpn;
  short_strikes = NULL;
  fix_start_dates = (long *)calloc(am->fix_leg->num_cpn, sizeof(long));
  fix_end_dates = (long *)calloc(am->fix_leg->num_cpn, sizeof(long));
  fixNotional = (double *)calloc(am->fix_leg->num_cpn, sizeof(double));
  long_strikes = (double *)calloc(am->fix_leg->num_cpn, sizeof(double));
  ex_fee = (double *)calloc(am->fix_leg->num_cpn, sizeof(double));
  if ((!fixNotional) || (!long_strikes) || (!ex_fee)) {
    err = "memory allocation failed in am_calib_und";
    smessage("memory allocation failed in am_calib_und");
    goto FREE_RETURN;
  }
  for (i = 0; i < am->fix_leg->num_cpn; ++i) {
    fix_start_dates[i] = am->fix_leg->cpn[i].start_date;
    fix_end_dates[i] = am->fix_leg->cpn[i].end_date;
    fixNotional[i] = am->fix_leg->notional[i];
    long_strikes[i] = am->fix_leg->rate[i];
    ex_fee[i] = am->fix_leg->fee[i];
  }

  Nfloat = am->fund_leg->num_cpn;
  float_start_dates = (long *)calloc(am->fund_leg->num_cpn, sizeof(long));
  float_end_dates = (long *)calloc(am->fund_leg->num_cpn, sizeof(long));
  floatNotional = (double *)calloc(am->fund_leg->num_cpn, sizeof(double));
  float_margin = (double *)calloc(am->fund_leg->num_cpn, sizeof(double));
  float_spread = (double *)calloc(am->fund_leg->num_cpn, sizeof(double));
  if ((!floatNotional) || (!float_margin) || (!float_spread)) {
    err = "memory allocation failed in am_calib_und";
    smessage("memory allocation failed in am_calib_und");
    goto FREE_RETURN;
  }
  for (i = 0; i < am->fund_leg->num_cpn; ++i) {
    float_start_dates[i] = am->fund_leg->cpn[i].start_date;
    float_end_dates[i] = am->fund_leg->cpn[i].end_date;
    floatNotional[i] = am->fund_leg->notional[i];
    float_margin[i] = am->fund_leg->margin[i];
    float_spread[i] = am->fund_leg->spread[i];
  }

  start_date = am->fix_leg->cpn[0].start_date;
  end_date = am->theoEndDate;

  err = amortMidat_cpd_calib_diagonal(
      notperiod, yc, vc, ref,

      get_cash_vol,

      0, 0,

      ex_bool,

      start_date, end_date,

      long_strikes, short_strikes,

      strike_type,

      diag_prices, ex_fee, fixNotional, floatNotional, float_margin,

      max_std_short, ref, swap_freq, swap_basis,

      fix_lambda, one_f_equi,

      skip_last,

      max_var_jump,

      &lambda,

      one2F,

      alpha, gamma, rho, &(und->sigma_n), &(und->sigma_time), &(und->sigma),
      &(und->inst_data));

  if (err) {
    goto FREE_RETURN;
  }

  und->lambda = lambda;

  /* write default values */
  und->lambda_n = 1;
  und->pdlambda = 0;
  und->pdlambda_date = 0;
  und->pdlambda_time = 0;

  und->has_inst_data = 1;

  und->sigma_date = (double *)calloc(und->sigma_n, sizeof(double));
  if (!und->sigma_date) {
    err = "Allocation error (6) in am_calib_und";
    goto FREE_RETURN;
  }

  for (i = 0; i < und->sigma_n; i++) {
    und->sigma_date[i] = today + und->sigma_time[i] * DAYS_IN_YEAR + 1.0e-08;
  }

FREE_RETURN:

  if (err)
    am_free_und(und);

  if (fix_start_dates)
    free(fix_start_dates);
  if (fix_end_dates)
    free(fix_end_dates);
  if (fixNotional)
    free(fixNotional);
  if (long_strikes)
    free(long_strikes);
  if (ex_fee)
    free(ex_fee);

  if (float_start_dates)
    free(float_start_dates);
  if (float_end_dates)
    free(float_end_dates);
  if (floatNotional)
    free(floatNotional);
  if (float_margin)
    free(float_margin);
  if (float_spread)
    free(float_spread);

  if (lex)
    free(lex);
  if (short_strikes)
    free(short_strikes);

  if (ex_bool)
    free(ex_bool);
  if (diag_prices)
    free(diag_prices);

  return err;
}

/*	Fill underlying structure from calibration instruments */
Err am_calib_und_new(
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I  , 1: E */

    char *yc,                 /*	yc */
    char *vc,                 /*	vc */
    char *default_ref,        /*	ref rate */
    char *default_swap_freq,  /*	swap freq */
    char *default_swap_basis, /*	swap basis */

    char *ref,        /*	ref rate */
    char *swap_freq,  /*	swap freq */
    char *swap_basis, /*	swap basis */

    double lambda, /*	lambda if unique */
    double alpha,  /*	alpha */
    double gamma,  /*	gamma */
    double rho,    /*	rho */
    /*	Calib params */
    double mintime, double mininterval,

    int notperiod, double max_std_short, int one2F,
    int fix_lambda, /*	0: calib lambda to cap  , 1: fix lambda calib
                                                    to diagonal */
    int one_f_equi, /*	1F equivalent flag:
                                                    if set to 1  , then 2F
                       lambda will calibrate to the cap priced within calibrated
                       1F with the given lambda */
    int skip_last,  /*	If 1  , the last option is disregarded
                                                    and the forward volatility is
                       flat from option  n-1 */
    int use_jump, double max_var_jump,

    int strike_type, int european_model,

    Err (*get_correl)(char *correl_cube_name, double start_date,
                      double end_date, double strike, double *vol),
    char *CorrelName,
    /*
                    Err  (*GetVolForBadr)( Date  , Date  , double  , SRT_Boolean
       , double *)  , char *cVolType  ,
    */
    /*	End of calib params */
    AM_STR am,          /*	structure */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    AM_UND und) {
  int i, nex, firstex;
  long *lex = NULL;

  long end_struct; //  ,
                   // end_swap;
  double *long_strikes = NULL, *short_strikes = NULL;

  int Nfix;
  int *ex_bool = NULL;
  double *diag_prices = NULL;
  double *ex_fee = NULL;
  double *fixNotional = NULL;
  long *fix_start_dates = NULL;
  long *fix_end_dates = NULL;
  long *fix_pay_dates = NULL;

  int Nfloat;
  double *floatNotional = NULL;
  long *float_start_dates = NULL;
  long *float_end_dates = NULL;
  long *float_pay_dates = NULL;
  double *float_margin = NULL;
  double *float_spread = NULL;

  long start_date;
  long end_date;

  int UseVol;

  SrtCompounding floatfreq;
  SrtBasisCode floatbasis;

  char *float_freq = NULL;
  char *float_basis = NULL;

  AM_FUND_LEG fund_leg;
  AM_FIX_LEG fix_leg;

  Err err = NULL;

  /*	Initialise */
  fix_leg = am->fix_leg;
  fund_leg = am->fund_leg;

  und->sigma_date = NULL;
  und->sigma_time = NULL;
  und->sigma = NULL;

  und->has_inst_data = 0;
  cpd_init_calib_inst_data(&(und->inst_data));

  und->num_prices = 0;
  und->exercise_dates = NULL;
  und->mkt_prices = NULL;
  und->mdl_prices = NULL;

  und->lambda = lambda;

  /* write default values */
  und->lambda_n = 1;
  und->pdlambda = 0;
  und->pdlambda_date = 0;
  und->pdlambda_time = 0;

  und->alpha = alpha;
  und->gamma = gamma;
  und->rho = rho;
  und->today = today;

  strcpy(und->yc, yc);
  strcpy(und->vc, vc);
  strcpy(und->ref, ref);
  strcpy(und->swap_freq, swap_freq);
  strcpy(und->swap_basis, swap_basis);
  strcpy(und->name, "CALIB");

  /*	Exercise dates for calibration */
  end_struct = fund_leg->cpn[fund_leg->num_cpn - 1].pay_date;
  //	end_swap = end_struct;
  //	if (fix_leg->cpn[fix_leg->num_cpn-1].pay_date > end_swap)
  //	{
  //		end_swap = fix_leg->cpn[fix_leg->num_cpn-1].pay_date;
  //	}

  if (am->num_calls > 0 &&
      !(am->num_calls == 1 && am->call[0].ex_date <= today + eod_flag)) {
    /*	If call dates  , choose call dates as option expiries for calibration */
    nex = am->num_calls;
    lex = (long *)calloc(nex, sizeof(long));
    if (!lex) {
      err = "Allocation error (2) in am_calib_und";
      goto FREE_RETURN;
    }

    for (i = 0; i < nex; i++) {
      lex[i] = am->call[i].ex_date;
    }
  } else {
    /*	If no call dates  , exit */
    am->num_calls = 0;
    goto FREE_RETURN;
  }

  UseVol = 0;
  if (strike_type == 5) {
    UseVol = 1;
  }

  err = amortMidat_compute_diagonal_prices_new(
      yc, vc, ref, swap_freq, swap_basis, today, eod_flag, mintime, mininterval,
      UseVol, &nex, &firstex, &ex_bool, &diag_prices, european_model,
      get_correl, CorrelName, get_cash_vol,
      /*GetVolForBadr  , cVolType  ,*/ am);
  if (err) {
    goto FREE_RETURN;
  }

  und->num_prices = nex;
  und->exercise_dates = (double *)calloc(nex, sizeof(double));
  und->mkt_prices = (double *)calloc(nex, sizeof(double));
  und->mdl_prices = (double *)calloc(nex, sizeof(double));
  for (i = 0; i < nex; ++i) {
    und->exercise_dates[i] = am->fix_leg->cpn[firstex + i].start_date;
    und->mkt_prices[i] = diag_prices[i];
  }

  Nfix = am->fix_leg->num_cpn;

  //// Changed by Albert Wang on 12/19/03
  //// previously
  // short_strikes = NULL;
  short_strikes = (double *)calloc(am->fix_leg->num_cpn, sizeof(double));

  fix_start_dates = (long *)calloc(am->fix_leg->num_cpn, sizeof(long));
  fix_end_dates = (long *)calloc(am->fix_leg->num_cpn, sizeof(long));
  fix_pay_dates = (long *)calloc(am->fix_leg->num_cpn, sizeof(long));
  fixNotional = (double *)calloc(am->fix_leg->num_cpn, sizeof(double));
  long_strikes = (double *)calloc(am->fix_leg->num_cpn, sizeof(double));
  ex_fee = (double *)calloc(am->fix_leg->num_cpn, sizeof(double));
  if ((!fixNotional) || (!long_strikes) || (!ex_fee) || (!fix_start_dates) ||
      (!fix_end_dates) || (!fix_pay_dates)) {
    err = "memory allocation failed in am_calib_und";
    smessage("memory allocation failed in am_calib_und");
    goto FREE_RETURN;
  }
  for (i = 0; i < am->fix_leg->num_cpn; ++i) {
    fix_start_dates[i] = am->fix_leg->cpn[i].start_date;
    fix_end_dates[i] = am->fix_leg->cpn[i].end_date;
    fix_pay_dates[i] = am->fix_leg->cpn[i].end_date;
    fixNotional[i] = am->fix_leg->notional[i];
    long_strikes[i] = am->fix_leg->rate[i];
    short_strikes[i] = am->fix_leg->rate[i];
    ex_fee[i] = am->fix_leg->fee[i];
  }

  Nfloat = am->fund_leg->num_cpn;
  float_start_dates = (long *)calloc(am->fund_leg->num_cpn, sizeof(long));
  float_end_dates = (long *)calloc(am->fund_leg->num_cpn, sizeof(long));
  float_pay_dates = (long *)calloc(am->fund_leg->num_cpn, sizeof(long));
  floatNotional = (double *)calloc(am->fund_leg->num_cpn, sizeof(double));
  float_margin = (double *)calloc(am->fund_leg->num_cpn, sizeof(double));
  float_spread = (double *)calloc(am->fund_leg->num_cpn, sizeof(double));
  if ((!floatNotional) || (!float_margin) || (!float_spread) ||
      (!float_start_dates) || (!float_end_dates) || (!float_pay_dates)) {
    err = "memory allocation failed in am_calib_und";
    smessage("memory allocation failed in am_calib_und");
    goto FREE_RETURN;
  }
  for (i = 0; i < am->fund_leg->num_cpn; ++i) {
    float_start_dates[i] = am->fund_leg->cpn[i].start_date;
    float_end_dates[i] = am->fund_leg->cpn[i].end_date;
    float_pay_dates[i] = am->fund_leg->cpn[i].end_date;
    floatNotional[i] = am->fund_leg->notional[i];
    float_margin[i] = am->fund_leg->margin[i];
    float_spread[i] = am->fund_leg->spread[i];
  }

  start_date = am->fix_leg->cpn[0].start_date;
  end_date = am->theoEndDate;

  /*	err = amortMidat_cpd_calib_diagonal(
                                  notperiod  ,
                                  yc  ,
                                  vc  ,
                                  ref  ,

                                  get_cash_vol  ,

                                  0  ,
                                  0  ,

                                  ex_bool  ,

                                  start_date  ,
                                  end_date  ,

                                  long_strikes  ,
                                  short_strikes  ,

                                  strike_type  ,

                                  diag_prices  ,
                                  ex_fee  ,
                                  fixNotional  ,
                                  floatNotional  ,
                                  float_margin  ,

                                  max_std_short  ,
                                  ref  ,
                                  swap_freq  ,
                                  swap_basis  ,

                                  fix_lambda  ,
                                  one_f_equi  ,

                                  skip_last  ,

                                  max_var_jump  ,

                                  &lambda  ,

                                  one2F  ,

                                  alpha  ,
                                  gamma  ,
                                  rho  ,
                                  &(und->sigma_n)  ,
                                  &(und->sigma_time)  ,
                                  &(und->sigma)  ,
                                  &(und->inst_data));
  */

  err = swp_f_get_ref_rate_details(ref, &floatbasis, &floatfreq);
  err = translate_compounding(&float_freq, floatfreq);
  err = translate_basis(&float_basis, floatbasis);
  if (err) {
    goto FREE_RETURN;
  }

  err = amortMidat_cpd_calib_diagonal_new(
      notperiod, yc, vc, default_ref, default_swap_basis, default_swap_freq,

      get_cash_vol, // GetCpdAutocalCashVol

      0,
      1, /// shift type

      ex_bool, diag_prices, ex_fee,

      swap_freq, swap_basis, Nfix, fix_start_dates, fix_end_dates,
      fix_pay_dates, long_strikes, fixNotional,

      float_freq, float_basis, Nfloat, float_start_dates, float_end_dates,
      float_pay_dates, float_margin, float_spread, floatNotional,

      short_strikes,

      strike_type,

      max_std_short,

      fix_lambda, one_f_equi,

      skip_last,

      use_jump, max_var_jump,

      &lambda,

      one2F,

      alpha, gamma, rho, &(und->sigma_n), &(und->sigma_time), &(und->sigma));
  if (err) {
    goto FREE_RETURN;
  }

  und->lambda = lambda;

  /* write default values */
  und->lambda_n = 1;
  und->pdlambda = 0;
  und->pdlambda_date = 0;
  und->pdlambda_time = 0;

  und->has_inst_data = 1;

  und->sigma_date = (double *)calloc(und->sigma_n, sizeof(double));
  if (!und->sigma_date) {
    err = "Allocation error (6) in am_calib_und";
    goto FREE_RETURN;
  }

  for (i = 0; i < und->sigma_n; i++) {
    und->sigma_date[i] = today + und->sigma_time[i] * DAYS_IN_YEAR + 1.0e-08;
  }

FREE_RETURN:

  if (err)
    am_free_und(und);

  //	if (float_freq) free(float_freq);
  //	if (float_basis) free(float_basis);

  if (fix_start_dates)
    free(fix_start_dates);
  if (fix_end_dates)
    free(fix_end_dates);
  if (fix_pay_dates)
    free(fix_pay_dates);
  if (fixNotional)
    free(fixNotional);
  if (long_strikes)
    free(long_strikes);
  if (ex_fee)
    free(ex_fee);

  if (float_start_dates)
    free(float_start_dates);
  if (float_end_dates)
    free(float_end_dates);
  if (float_pay_dates)
    free(float_pay_dates);
  if (floatNotional)
    free(floatNotional);
  if (float_margin)
    free(float_margin);
  if (float_spread)
    free(float_spread);

  if (lex)
    free(lex);
  if (short_strikes)
    free(short_strikes);

  if (ex_bool)
    free(ex_bool);
  if (diag_prices)
    free(diag_prices);

  return err;
}

/*	Fill underlying structure from calibration instruments */
Err am_calib_und_new_ts(
    long today,
    /*	EOD Flag */
    int eod_flag, /*	0: I  , 1: E */

    char *yc,                 /*	yc */
    char *vc,                 /*	vc */
    char *default_ref,        /*	ref rate */
    char *default_swap_freq,  /*	swap freq */
    char *default_swap_basis, /*	swap basis */

    char *ref,        /*	ref rate */
    char *swap_freq,  /*	swap freq */
    char *swap_basis, /*	swap basis */

    int nlambda, double *pdlambda_time, double *pdlambda,

    double alpha, /*	alpha */
    double gamma, /*	gamma */
    double rho,   /*	rho */
    /*	Calib params */
    double mintime, double mininterval,

    int notperiod, double max_std_short, int one2F,
    int fix_lambda, /*	0: calib lambda to cap  , 1: fix lambda calib
                                                    to diagonal */
    int one_f_equi, /*	1F equivalent flag:
                                                    if set to 1  , then 2F
                       lambda will calibrate to the cap priced within calibrated
                       1F with the given lambda */
    int skip_last,  /*	If 1  , the last option is disregarded
                                                    and the forward volatility is
                       flat from option  n-1 */
    int use_jump, double max_var_jump,

    int strike_type, int european_model,

    Err (*get_correl)(char *correl_cube_name, double start_date,
                      double end_date, double strike, double *vol),
    char *CorrelName,
    /*
                    Err  (*GetVolForBadr)( Date  , Date  , double  , SRT_Boolean
       , double *)  , char *cVolType  ,
    */
    /*	End of calib params */
    AM_STR am,          /*	structure */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    AM_UND und) {
  int i, nex, firstex;
  long *lex = NULL;

  long end_struct; //  ,
                   // end_swap;
  double *long_strikes = NULL, *short_strikes = NULL;

  int Nfix;
  int *ex_bool = NULL;
  double *diag_prices = NULL;
  double *ex_fee = NULL;
  double *fixNotional = NULL;
  long *fix_start_dates = NULL;
  long *fix_end_dates = NULL;
  long *fix_pay_dates = NULL;

  int Nfloat;
  double *floatNotional = NULL;
  long *float_start_dates = NULL;
  long *float_end_dates = NULL;
  long *float_pay_dates = NULL;
  double *float_margin = NULL;
  double *float_spread = NULL;

  long start_date;
  long end_date;

  int UseVol;

  SrtCompounding floatfreq;
  SrtBasisCode floatbasis;

  char *float_freq = NULL;
  char *float_basis = NULL;

  AM_FUND_LEG fund_leg;
  AM_FIX_LEG fix_leg;

  Err err = NULL;

  /*	Initialise */
  fix_leg = am->fix_leg;
  fund_leg = am->fund_leg;

  und->sigma_date = NULL;
  und->sigma_time = NULL;
  und->sigma = NULL;

  und->has_inst_data = 0;
  cpd_init_calib_inst_data(&(und->inst_data));

  und->num_prices = 0;
  und->exercise_dates = NULL;
  und->mkt_prices = NULL;
  und->mdl_prices = NULL;

  // und->lambda = lambda;

  und->alpha = alpha;
  und->gamma = gamma;
  und->rho = rho;
  und->today = today;

  strcpy(und->yc, yc);
  strcpy(und->vc, vc);
  strcpy(und->ref, ref);
  strcpy(und->swap_freq, swap_freq);
  strcpy(und->swap_basis, swap_basis);
  strcpy(und->name, "CALIB");

  /*	Exercise dates for calibration */
  end_struct = fund_leg->cpn[fund_leg->num_cpn - 1].pay_date;
  //	end_swap = end_struct;
  //	if (fix_leg->cpn[fix_leg->num_cpn-1].pay_date > end_swap)
  //	{
  //		end_swap = fix_leg->cpn[fix_leg->num_cpn-1].pay_date;
  //	}

  if (am->num_calls > 0 &&
      !(am->num_calls == 1 && am->call[0].ex_date <= today + eod_flag)) {
    /*	If call dates  , choose call dates as option expiries for calibration */
    nex = am->num_calls;
    lex = (long *)calloc(nex, sizeof(long));
    if (!lex) {
      err = "Allocation error (2) in am_calib_und";
      goto FREE_RETURN;
    }

    for (i = 0; i < nex; i++) {
      lex[i] = am->call[i].ex_date;
    }
  } else {
    /*	If no call dates  , exit */
    am->num_calls = 0;
    goto FREE_RETURN;
  }

  UseVol = 0;
  if (strike_type == 5) {
    UseVol = 1;
  }

  err = amortMidat_compute_diagonal_prices_new(
      yc, vc, ref, swap_freq, swap_basis, today, eod_flag, mintime, mininterval,
      UseVol, &nex, &firstex, &ex_bool, &diag_prices, european_model,
      get_correl, CorrelName, get_cash_vol,
      /*GetVolForBadr  , cVolType  ,*/ am);
  if (err) {
    goto FREE_RETURN;
  }

  und->num_prices = nex;
  und->exercise_dates = (double *)calloc(nex, sizeof(double));
  und->mkt_prices = (double *)calloc(nex, sizeof(double));
  und->mdl_prices = (double *)calloc(nex, sizeof(double));
  for (i = 0; i < nex; ++i) {
    und->exercise_dates[i] = am->fix_leg->cpn[firstex + i].start_date;
    und->mkt_prices[i] = diag_prices[i];
  }

  Nfix = am->fix_leg->num_cpn;
  //	short_strikes = NULL;
  fix_start_dates = (long *)calloc(am->fix_leg->num_cpn, sizeof(long));
  fix_end_dates = (long *)calloc(am->fix_leg->num_cpn, sizeof(long));
  fix_pay_dates = (long *)calloc(am->fix_leg->num_cpn, sizeof(long));
  fixNotional = (double *)calloc(am->fix_leg->num_cpn, sizeof(double));
  long_strikes = (double *)calloc(am->fix_leg->num_cpn, sizeof(double));
  ex_fee = (double *)calloc(am->fix_leg->num_cpn, sizeof(double));
  if ((!fixNotional) || (!long_strikes) || (!ex_fee) || (!fix_start_dates) ||
      (!fix_end_dates) || (!fix_pay_dates)) {
    err = "memory allocation failed in am_calib_und";
    smessage("memory allocation failed in am_calib_und");
    goto FREE_RETURN;
  }
  for (i = 0; i < am->fix_leg->num_cpn; ++i) {
    fix_start_dates[i] = am->fix_leg->cpn[i].start_date;
    fix_end_dates[i] = am->fix_leg->cpn[i].end_date;
    fix_pay_dates[i] = am->fix_leg->cpn[i].end_date;
    fixNotional[i] = am->fix_leg->notional[i];
    long_strikes[i] = am->fix_leg->rate[i];
    ex_fee[i] = am->fix_leg->fee[i];
  }

  Nfloat = am->fund_leg->num_cpn;
  float_start_dates = (long *)calloc(am->fund_leg->num_cpn, sizeof(long));
  float_end_dates = (long *)calloc(am->fund_leg->num_cpn, sizeof(long));
  float_pay_dates = (long *)calloc(am->fund_leg->num_cpn, sizeof(long));
  floatNotional = (double *)calloc(am->fund_leg->num_cpn, sizeof(double));
  float_margin = (double *)calloc(am->fund_leg->num_cpn, sizeof(double));
  float_spread = (double *)calloc(am->fund_leg->num_cpn, sizeof(double));
  if ((!floatNotional) || (!float_margin) || (!float_spread) ||
      (!float_start_dates) || (!float_end_dates) || (!float_pay_dates)) {
    err = "memory allocation failed in am_calib_und";
    smessage("memory allocation failed in am_calib_und");
    goto FREE_RETURN;
  }
  for (i = 0; i < am->fund_leg->num_cpn; ++i) {
    float_start_dates[i] = am->fund_leg->cpn[i].start_date;
    float_end_dates[i] = am->fund_leg->cpn[i].end_date;
    float_pay_dates[i] = am->fund_leg->cpn[i].end_date;
    floatNotional[i] = am->fund_leg->notional[i];
    float_margin[i] = am->fund_leg->margin[i];
    float_spread[i] = am->fund_leg->spread[i];
  }

  start_date = am->fix_leg->cpn[0].start_date;
  end_date = am->theoEndDate;

  /*	err = amortMidat_cpd_calib_diagonal(
                                  notperiod  ,
                                  yc  ,
                                  vc  ,
                                  ref  ,

                                  get_cash_vol  ,

                                  0  ,
                                  0  ,

                                  ex_bool  ,

                                  start_date  ,
                                  end_date  ,

                                  long_strikes  ,
                                  short_strikes  ,

                                  strike_type  ,

                                  diag_prices  ,
                                  ex_fee  ,
                                  fixNotional  ,
                                  floatNotional  ,
                                  float_margin  ,

                                  max_std_short  ,
                                  ref  ,
                                  swap_freq  ,
                                  swap_basis  ,

                                  fix_lambda  ,
                                  one_f_equi  ,

                                  skip_last  ,

                                  max_var_jump  ,

                                  &lambda  ,

                                  one2F  ,

                                  alpha  ,
                                  gamma  ,
                                  rho  ,
                                  &(und->sigma_n)  ,
                                  &(und->sigma_time)  ,
                                  &(und->sigma)  ,
                                  &(und->inst_data));
  */

  err = swp_f_get_ref_rate_details(ref, &floatbasis, &floatfreq);
  err = translate_compounding(&float_freq, floatfreq);
  err = translate_basis(&float_basis, floatbasis);
  if (err) {
    goto FREE_RETURN;
  }

  err = amortMidat_cpd_calib_diagonal_new_ts(
      notperiod, yc, vc, default_ref, default_swap_basis, default_swap_freq,

      get_cash_vol, // GetCpdAutocalCashVol

      0, 1,

      ex_bool, diag_prices, ex_fee,

      swap_freq, swap_basis, Nfix, fix_start_dates, fix_end_dates,
      fix_pay_dates, long_strikes, fixNotional,

      float_freq, float_basis, Nfloat, float_start_dates, float_end_dates,
      float_pay_dates, float_margin, float_spread, floatNotional,

      short_strikes,

      strike_type,

      max_std_short,

      fix_lambda, one_f_equi,

      skip_last,

      use_jump, max_var_jump,

      nlambda, pdlambda, pdlambda_time,

      one2F,

      alpha, gamma, rho, &(und->sigma_n), &(und->sigma_time),
      &(und->sigma)); //// 121103
  if (err) {
    goto FREE_RETURN;
  }

  und->lambda_n = nlambda;
  und->pdlambda = (double *)calloc(und->lambda_n, sizeof(double));
  und->pdlambda_time = (double *)calloc(und->lambda_n, sizeof(double));
  und->pdlambda_date = (double *)calloc(und->lambda_n, sizeof(double));
  for (i = 0; i < und->lambda_n; ++i) {
    und->pdlambda[i] = pdlambda[i];
    und->pdlambda_time[i] = pdlambda_time[i];
    und->pdlambda_date[i] =
        today + und->pdlambda_time[i] * DAYS_IN_YEAR + 1.0e-08;
  }

  und->has_inst_data = 1;

  und->sigma_date = (double *)calloc(und->sigma_n, sizeof(double));
  if (!und->sigma_date) {
    err = "Allocation error (6) in am_calib_und";
    goto FREE_RETURN;
  }

  for (i = 0; i < und->sigma_n; i++) {
    und->sigma_date[i] = today + und->sigma_time[i] * DAYS_IN_YEAR + 1.0e-08;
  }

FREE_RETURN:

  if (err)
    am_free_und(und);

  //	if (float_freq) free(float_freq);
  //	if (float_basis) free(float_basis);

  if (fix_start_dates)
    free(fix_start_dates);
  if (fix_end_dates)
    free(fix_end_dates);
  if (fix_pay_dates)
    free(fix_pay_dates);
  if (fixNotional)
    free(fixNotional);
  if (long_strikes)
    free(long_strikes);
  if (ex_fee)
    free(ex_fee);

  if (float_start_dates)
    free(float_start_dates);
  if (float_end_dates)
    free(float_end_dates);
  if (float_pay_dates)
    free(float_pay_dates);
  if (floatNotional)
    free(floatNotional);
  if (float_margin)
    free(float_margin);
  if (float_spread)
    free(float_spread);

  if (lex)
    free(lex);
  if (short_strikes)
    free(short_strikes);

  if (ex_bool)
    free(ex_bool);
  if (diag_prices)
    free(diag_prices);

  return err;
}

void am_copy_und(AM_UND src, AM_UND dest) {
  int i;
  strcpy(dest->name, src->name);
  dest->today = src->today;
  strcpy(dest->yc, src->yc);
  strcpy(dest->vc, src->vc);
  strcpy(dest->ref, src->ref);
  strcpy(dest->swap_freq, src->swap_freq);
  strcpy(dest->swap_basis, src->swap_basis);

  dest->sigma_n = src->sigma_n;

  if (src->sigma_date) {
    dest->sigma_date = (double *)calloc(dest->sigma_n, sizeof(double));
    memcpy(dest->sigma_date, src->sigma_date, dest->sigma_n * sizeof(double));
  } else {
    dest->sigma_date = NULL;
  }

  if (src->sigma_time) {
    dest->sigma_time = (double *)calloc(dest->sigma_n, sizeof(double));
    memcpy(dest->sigma_time, src->sigma_time, dest->sigma_n * sizeof(double));
  } else {
    dest->sigma_time = NULL;
  }

  if (src->sigma) {
    dest->sigma = (double *)calloc(dest->sigma_n, sizeof(double));
    memcpy(dest->sigma, src->sigma, dest->sigma_n * sizeof(double));
  } else {
    dest->sigma = NULL;
  }

  if (src->has_inst_data) {
    cpd_copy_calib_inst_data(&(dest->inst_data), &(src->inst_data));
    dest->has_inst_data = 1;
  } else {
    dest->has_inst_data = 0;
  }

  dest->lambda_n = src->lambda_n;
  dest->lambda = src->lambda;

  if (dest->lambda_n > 1) {
    dest->pdlambda = (double *)calloc(dest->lambda_n, sizeof(double));
    dest->pdlambda_time = (double *)calloc(dest->lambda_n, sizeof(double));
    dest->pdlambda_date = (double *)calloc(dest->lambda_n, sizeof(double));

    for (i = 0; i < dest->lambda_n; ++i) {
      dest->pdlambda[i] = src->pdlambda[i];
      dest->pdlambda_time[i] = src->pdlambda_time[i];
      dest->pdlambda_date[i] = src->pdlambda_date[i];
    }
  }

  dest->alpha = src->alpha;
  dest->gamma = src->gamma;
  dest->rho = src->rho;
}

Err am_free_und(AM_UND und) {
  if (und->sigma_date)
    free(und->sigma_date);
  if (und->sigma_time)
    free(und->sigma_time);
  if (und->sigma)
    free(und->sigma);

  und->sigma_date = NULL;
  und->sigma_time = NULL;
  und->sigma = NULL;

  if (und->has_inst_data) {
    cpd_free_calib_inst_data(&und->inst_data);
    und->has_inst_data = 0;
  }

  if (und->num_prices > 0) {
    if (und->mkt_prices)
      free(und->mkt_prices);
    if (und->mdl_prices)
      free(und->mdl_prices);
    if (und->exercise_dates)
      free(und->exercise_dates);
    und->num_prices = 0;
  }

  return NULL;
}

/*	Function for constants used for reconstruction and evaluation  */

/*	Just round a double to the closest long */
static long strnd(double x) {
  static long y;
  y = (long)x;
  if (fabs(y + 1 - x) < fabs(y - x) && fabs(y + 1 - x) < fabs(y - 1 - x)) {
    return y + 1;
  } else if (fabs(y - x) < fabs(y - 1 - x)) {
    return y;
  } else {
    return y - 1;
  }
}

/*	THE function */
Err am_fill_eval_const(AM_UND und, AM_STR am,
                       /*	Index of the current call */
                       int call_idx, AM_EVAL_CONST eval_const) {
  AM_CALL call;
  AM_FUND_LEG fund_leg;
  AM_FUND_CPN fund_cpn;
  AM_FIX_LEG fix_leg;
  AM_FIX_CPN fix_cpn;
  long evt_date;
  double evt_time;
  double lam1, lam2, phi1, phi2, phi12;
  long today;
  long temp_date;
  double temp_time;
  int tmpidx;
  char *yc;
  int i, j;
  double t1, t2, dt;
  double gam1, gam2, gam12, fact1, fact1temp, fact2, fact2temp;
  Err err = NULL;
  double lambda_time;
  double *lambda_time_vect = NULL;
  double *lambda_vect = NULL;

  /*	Initialise
          ----------	*/

  /*	Extract info from call */
  call = am->call + call_idx;
  evt_date = call->ex_date;
  evt_time = call->ex_time;

  /*	Extract info for structure */
  fund_leg = am->fund_leg;
  fix_leg = am->fix_leg;

  /*	Extract info from und */
  today = und->today;
  yc = (char *)(und->yc);

  lambda_time = 0;

  lambda_time_vect = (double *)calloc(1, sizeof(double));
  lambda_vect = (double *)calloc(1, sizeof(double));
  if ((!lambda_vect) || (!lambda_time_vect)) {
    err = "Memory allocation failed in am_fill_eval_const";
    goto FREE_RETURN;
  }
  lambda_time_vect[0] = 1.0;
  lambda_vect[0] = und->lambda;

  err = LGM2FDtsPhi(evt_time, und->sigma_time, und->sigma_n, und->sigma,
                    lambda_time_vect, 1, lambda_vect, und->alpha, und->gamma,
                    und->rho, &phi1, &phi2, &phi12);

  /*	Get the maturities of the required discount factors
          --------------------------------------------------- */
  eval_const->num_df = 0;

  /*	Funding leg */
  if (call->num_fund_cpn > 0) {
    eval_const->do_fund = 1;
    eval_const->num_fund_cpn = call->num_fund_cpn;
    for (i = call->fund_idx; i < fund_leg->num_cpn; i++) {
      fund_cpn = fund_leg->cpn + i;
      eval_const->df_mat[eval_const->num_df] = fund_cpn->pay_time - evt_time;
      (eval_const->num_df)++;
    }
  } else {
    eval_const->do_fund = 0;
    eval_const->num_fund_cpn = 0;
    err = "No funding coupons left at call date";
    goto FREE_RETURN;
  }
  eval_const->df_mat[eval_const->num_df] =
      fund_leg->cpn[call->fund_idx].start_time - evt_time;
  (eval_const->num_df)++;

  /*	Fix Leg */
  if (call->num_fix_cpn > 0) {
    eval_const->do_fix_disc = 1;
    eval_const->do_fix_fwd = 1;
    eval_const->num_fix_cpn = call->num_fix_cpn;

    /*	Discounting */
    for (i = call->fix_idx; i < fix_leg->num_cpn; i++) {
      fix_cpn = fix_leg->cpn + i;
      eval_const->df_mat[eval_const->num_df] = fix_cpn->pay_time - evt_time;
      (eval_const->num_df)++;
    }

  } else {
    eval_const->do_fix_disc = 0;
    eval_const->do_fix_fwd = 0;
    eval_const->num_fix_cpn = 0;
    err = "No fix coupons left at call date";
    goto FREE_RETURN;
  }

  /*	Fee	*/
  eval_const->df_mat[eval_const->num_df] = call->set_time - evt_time;
  (eval_const->num_df)++;

  /*	Sort and unique */
  num_f_sort_vector(eval_const->num_df, eval_const->df_mat);
  num_f_unique_vector(&(eval_const->num_df), eval_const->df_mat);

  /*	Precalculate DFF  , gamma  , 0.5 * gamma * gamma
          --------------------------------------------	*/

  /* Find first used lambda */
  j = 0;

  gam1 = gam2 = gam12 = 0.0;
  lam1 = und->lambda;
  lam2 = lam1 + und->gamma;

  fact1 = 0.0;
  fact2 = 0.0;

  t1 = evt_time;

  for (i = 0; i < eval_const->num_df; i++) {
    temp_time = eval_const->df_mat[i];
    temp_date = evt_date + strnd(temp_time * 365.0);
    t2 = evt_time + temp_time;

    gam1 = 0.0;
    gam2 = 0.0;
    fact1temp = 0.0;
    fact2temp = 0.0;

    dt = (t2 - t1);

    gam1 += exp(-fact1temp) * (1.0 - exp(-lam1 * dt)) / lam1;
    gam2 += exp(-fact2temp) * (1.0 - exp(-lam2 * dt)) / lam2;

    fact1temp += lam1 * dt;
    fact2temp += lam2 * dt;

    if (i > 0) {
      eval_const->df_beta[i] = eval_const->df_beta[i - 1] + exp(-fact1) * gam1;
      eval_const->df_gamma[i] =
          eval_const->df_gamma[i - 1] + exp(-fact2) * gam2;
    } else {
      eval_const->df_beta[i] = gam1;
      eval_const->df_gamma[i] = gam2;
    }

    eval_const->df_alpha[i] =
        swp_f_zr(evt_date, temp_date, yc) * temp_time +
        0.5 * (eval_const->df_beta[i] * eval_const->df_beta[i] * phi1 +
               eval_const->df_gamma[i] * eval_const->df_gamma[i] * phi2) +
        und->rho * eval_const->df_beta[i] * eval_const->df_gamma[i] * phi12;

    fact1 += fact1temp;
    fact2 += fact2temp;

    t1 = t2;
  }

  /*	Find all dfs indices
          --------------------	*/

  /*	Funding Leg */
  tmpidx = 0;
  temp_time = fund_leg->cpn[call->fund_idx].start_time - evt_time;
  while (eval_const->df_mat[tmpidx] < temp_time - 1.0e-08)
    tmpidx++;
  eval_const->start_idx = tmpidx;
  if (eval_const->do_fund) {
    for (i = call->fund_idx; i < fund_leg->num_cpn; i++) {
      fund_cpn = fund_leg->cpn + i;
      temp_time = fund_cpn->pay_time - evt_time;
      while (eval_const->df_mat[tmpidx] < temp_time - 1.0e-08)
        tmpidx++;
      eval_const->fund_idx[i - call->fund_idx] = tmpidx;
    }
  }

  /*	Fix Leg */
  if (eval_const->do_fix_disc && eval_const->do_fix_fwd) {
    /*	Discounting */
    tmpidx = 0;
    for (i = call->fix_idx; i < fix_leg->num_cpn; i++) {
      fix_cpn = fix_leg->cpn + i;
      temp_time = fix_cpn->pay_time - evt_time;
      while (eval_const->df_mat[tmpidx] < temp_time - 1.0e-08)
        tmpidx++;
      eval_const->fix_disc_idx[i - call->fix_idx] = tmpidx;
    }

    /*	Fix Coupons	*/
    tmpidx = 0;
    for (i = 0; i < call->num_fix_cpn; i++) {
      fix_cpn = fix_leg->cpn + call->fix_idx + i;
    }
  }

  /*	Fee	*/
  tmpidx = 0;
  temp_time = call->set_time - evt_time;
  while (eval_const->df_mat[tmpidx] < temp_time - 1.0e-08)
    tmpidx++;
  eval_const->fee_idx = tmpidx;

  /*	End of precalculations
          ---------------------- */

FREE_RETURN:

  if (lambda_time_vect) {
    free(lambda_time_vect);
    lambda_time_vect = NULL;
  }
  if (lambda_vect) {
    free(lambda_vect);
    lambda_vect = NULL;
  }

  return err;
}

Err am_fill_eval_const_ts(int nLambda, double *pdLambdaValue,
                          double *pdLambdaTime, AM_UND und, AM_STR am,
                          /*	Index of the current call */
                          int call_idx, AM_EVAL_CONST eval_const) {
  AM_CALL call;
  AM_FUND_LEG fund_leg;
  AM_FUND_CPN fund_cpn;
  AM_FIX_LEG fix_leg;
  AM_FIX_CPN fix_cpn;
  long evt_date;
  double evt_time;
  double // lam1  ,
         // lam2  ,
      phi1,
      phi2, phi12;
  double dAvgL1_t1, dAvgL1_t2, dAvgL2_t1, dAvgL2_t2, dAvgL1, dAvgL2;

  long today;
  long temp_date;
  double temp_time;
  int tmpidx;
  char *yc;
  int i, j;
  double t1, t2, dt;
  double gam1, gam2, gam12, fact1, fact1temp, fact2, fact2temp;
  Err err = NULL;

  /*	Initialise
          ----------	*/

  /*	Extract info from call */
  call = am->call + call_idx;
  evt_date = call->ex_date;
  evt_time = call->ex_time;

  /*	Extract info for structure */
  fund_leg = am->fund_leg;
  fix_leg = am->fix_leg;

  /*	Extract info from und */
  today = und->today;
  yc = (char *)(und->yc);

#if 0
	lambda_time = 0;
	
	lambda_time_vect = (double *) calloc(1  , sizeof(double));
	lambda_vect = (double *) calloc(1  , sizeof(double));
	if((!lambda_vect) || (!lambda_time_vect) )
	{
		err = "Memory allocation failed in am_fill_eval_const";
		goto FREE_RETURN;
	}
	
	lambda_time_vect[0] = 1.0;
	lambda_vect[0] = und->lambda;
#endif

  err = LGM2FDtsPhi(evt_time, und->sigma_time, und->sigma_n, und->sigma,

                    pdLambdaTime, nLambda, pdLambdaValue,

                    und->alpha, und->gamma, und->rho, &phi1, &phi2, &phi12);

  /*	Get the maturities of the required discount factors
          --------------------------------------------------- */
  eval_const->num_df = 0;

  /*	Funding leg */
  if (call->num_fund_cpn > 0) {
    eval_const->do_fund = 1;
    eval_const->num_fund_cpn = call->num_fund_cpn;
    for (i = call->fund_idx; i < fund_leg->num_cpn; i++) {
      fund_cpn = fund_leg->cpn + i;
      eval_const->df_mat[eval_const->num_df] = fund_cpn->pay_time - evt_time;
      (eval_const->num_df)++;
    }
  } else {
    eval_const->do_fund = 0;
    eval_const->num_fund_cpn = 0;
    err = "No funding coupons left at call date";
    goto FREE_RETURN;
  }
  eval_const->df_mat[eval_const->num_df] =
      fund_leg->cpn[call->fund_idx].start_time - evt_time;
  (eval_const->num_df)++;

  /*	Fix Leg */
  if (call->num_fix_cpn > 0) {
    eval_const->do_fix_disc = 1;
    eval_const->do_fix_fwd = 1;
    eval_const->num_fix_cpn = call->num_fix_cpn;

    /*	Discounting */
    for (i = call->fix_idx; i < fix_leg->num_cpn; i++) {
      fix_cpn = fix_leg->cpn + i;
      eval_const->df_mat[eval_const->num_df] = fix_cpn->pay_time - evt_time;
      (eval_const->num_df)++;
    }

  } else {
    eval_const->do_fix_disc = 0;
    eval_const->do_fix_fwd = 0;
    eval_const->num_fix_cpn = 0;
    err = "No fix coupons left at call date";
    goto FREE_RETURN;
  }

  /*	Fee	*/
  eval_const->df_mat[eval_const->num_df] = call->set_time - evt_time;
  (eval_const->num_df)++;

  /*	Sort and unique */
  num_f_sort_vector(eval_const->num_df, eval_const->df_mat);
  num_f_unique_vector(&(eval_const->num_df), eval_const->df_mat);

  /*	Precalculate DFF  , gamma  , 0.5 * gamma * gamma
          --------------------------------------------	*/

  /* Find first used lambda */
  j = 0;

  gam1 = gam2 = gam12 = 0.0;

  /// previously
  // lam1 = und->lambda;
  // lam2 = lam1 + und->gamma;

  fact1 = 0.0;
  fact2 = 0.0;

  t1 = evt_time;

  for (i = 0; i < eval_const->num_df; i++) {
    temp_time = eval_const->df_mat[i];
    temp_date = evt_date + strnd(temp_time * 365.0);
    t2 = evt_time + temp_time;

    gam1 = 0.0;
    gam2 = 0.0;
    fact1temp = 0.0;
    fact2temp = 0.0;

    dt = (t2 - t1);

    dAvgL1_t1 = _average_lambda_(t1, nLambda, pdLambdaTime, pdLambdaValue);
    dAvgL1_t2 = _average_lambda_(t2, nLambda, pdLambdaTime, pdLambdaValue);

    dAvgL2_t1 = dAvgL1_t1 + und->gamma;
    dAvgL2_t2 = dAvgL1_t2 + und->gamma;

    dAvgL1 = (dAvgL1_t2 * t2 - dAvgL1_t1 * t1) / dt;
    dAvgL2 = (dAvgL2_t2 * t2 - dAvgL2_t1 * t1) / dt;

    gam1 += exp(-fact1temp) * (1.0 - exp(-dAvgL1 * dt)) / dAvgL1;
    gam2 += exp(-fact2temp) * (1.0 - exp(-dAvgL2 * dt)) / dAvgL2;

    fact1temp += dAvgL1 * dt; // lam1 * dt;
    fact2temp += dAvgL2 * dt; // lam2 * dt;

    if (i > 0) {
      eval_const->df_beta[i] = eval_const->df_beta[i - 1] + exp(-fact1) * gam1;
      eval_const->df_gamma[i] =
          eval_const->df_gamma[i - 1] + exp(-fact2) * gam2;
    } else {
      eval_const->df_beta[i] = gam1;
      eval_const->df_gamma[i] = gam2;
    }

    eval_const->df_alpha[i] =
        swp_f_zr(evt_date, temp_date, yc) * temp_time +
        0.5 * (eval_const->df_beta[i] * eval_const->df_beta[i] * phi1 +
               eval_const->df_gamma[i] * eval_const->df_gamma[i] * phi2) +
        und->rho * eval_const->df_beta[i] * eval_const->df_gamma[i] * phi12;

    fact1 += fact1temp;
    fact2 += fact2temp;

    t1 = t2;
  }

  /*	Find all dfs indices
          --------------------	*/

  /*	Funding Leg */
  tmpidx = 0;
  temp_time = fund_leg->cpn[call->fund_idx].start_time - evt_time;
  while (eval_const->df_mat[tmpidx] < temp_time - 1.0e-08)
    tmpidx++;
  eval_const->start_idx = tmpidx;
  if (eval_const->do_fund) {
    for (i = call->fund_idx; i < fund_leg->num_cpn; i++) {
      fund_cpn = fund_leg->cpn + i;
      temp_time = fund_cpn->pay_time - evt_time;
      while (eval_const->df_mat[tmpidx] < temp_time - 1.0e-08)
        tmpidx++;
      eval_const->fund_idx[i - call->fund_idx] = tmpidx;
    }
  }

  /*	Fix Leg */
  if (eval_const->do_fix_disc && eval_const->do_fix_fwd) {
    /*	Discounting */
    tmpidx = 0;
    for (i = call->fix_idx; i < fix_leg->num_cpn; i++) {
      fix_cpn = fix_leg->cpn + i;
      temp_time = fix_cpn->pay_time - evt_time;
      while (eval_const->df_mat[tmpidx] < temp_time - 1.0e-08)
        tmpidx++;
      eval_const->fix_disc_idx[i - call->fix_idx] = tmpidx;
    }

    /*	Fix Coupons	*/
    tmpidx = 0;
    for (i = 0; i < call->num_fix_cpn; i++) {
      fix_cpn = fix_leg->cpn + call->fix_idx + i;
    }
  }

  /*	Fee	*/
  tmpidx = 0;
  temp_time = call->set_time - evt_time;
  while (eval_const->df_mat[tmpidx] < temp_time - 1.0e-08)
    tmpidx++;
  eval_const->fee_idx = tmpidx;

  /*	End of precalculations
          ---------------------- */

FREE_RETURN:

  // if(lambda_time_vect)
  //{
  //	free(lambda_time_vect);
  //	lambda_time_vect = NULL;
  //}
  // if(lambda_vect)
  //{
  //	free(lambda_vect);
  //	lambda_vect = NULL;
  //}

  return err;
}

Err am_fill_adi_arg_ts(
    int nLambda, double *pdLambdaValue, double *pdLambdaTime, AM_UND und,
    AM_STR am,
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	Required number of steps */
    int req_stp, int req_stpx, AM_ADI_ARG adi_arg) {
  Err err = NULL;
  AM_PAY_ARG am_prm;
  int i, j;

  /*	Initialise */
  adi_arg->time = NULL;
  adi_arg->date = NULL;

  adi_arg->void_prm = NULL;
  adi_arg->is_event = NULL;
  adi_arg->ifr = NULL;

  /*	Compute time steps */

  /*	Copy event dates */

  adi_arg->nstp = am->num_calls;
  adi_arg->nstpx = req_stpx;

  adi_arg->time = (double *)calloc(adi_arg->nstp, sizeof(double));
  if (!adi_arg->time) {
    err = "Memory allocation error (1) in am_fill_adi_arg";
    goto FREE_RETURN;
  }
  for (i = 0; i < adi_arg->nstp; i++) {
    adi_arg->time[i] = am->call[i].ex_time;
  }

  /*	Fill time vector */

  /*	Add today if required */
  if (adi_arg->time[0] < -EPS) {
    err = "Past event date in am_fill_adi_arg";
    goto FREE_RETURN;
  }
  if (adi_arg->time[0] > EPS) {
    num_f_add_number(&(adi_arg->nstp), &(adi_arg->time), 0.0);
    num_f_sort_vector(adi_arg->nstp, adi_arg->time);
    num_f_unique_vector(&(adi_arg->nstp), adi_arg->time);
  }

  /*	If only one event today  , add empty event */
  if (adi_arg->nstp == 1) {
    num_f_add_number(&(adi_arg->nstp), (&adi_arg->time), 1.0);
  }

  /*	Fill the vector */
  num_f_fill_vector_newalgo(&(adi_arg->nstp), &(adi_arg->time), req_stp);

  /*	Make dates */
  adi_arg->date = (double *)calloc(adi_arg->nstp, sizeof(double));
  if (!adi_arg->date) {
    err = "Memory allocation error (2) in am_fill_adi_arg";
    goto FREE_RETURN;
  }

  for (i = 0; i < adi_arg->nstp; i++) {
    adi_arg->date[i] = und->today + DAYS_IN_YEAR * adi_arg->time[i];

    if (i > 0 && adi_arg->date[i] - adi_arg->date[i - 1] >= 1) {
      adi_arg->date[i] = (long)(adi_arg->date[i] + 1.0e-08);
      adi_arg->time[i] = YEARS_IN_DAY * (adi_arg->date[i] - und->today);
    }
  }

  /*	Sigma */

  adi_arg->sig_time = und->sigma_time;
  adi_arg->sig1 = und->sigma;
  adi_arg->nb_sig = und->sigma_n;

  /*	Lambdas and ratio */
  adi_arg->lam = und->lambda;
  adi_arg->alpha = und->alpha;
  adi_arg->gamma = und->gamma;
  adi_arg->rho = und->rho;

  /*	Spot fx and yield curves */

  strcpy(adi_arg->yc, und->yc);

  /*	Fill distributions */

  adi_arg->ifr = (double *)calloc(adi_arg->nstp, sizeof(double));

  if (!adi_arg->ifr) {
    err = "Memory allocation error (3) in am_fill_adi_arg";
    goto FREE_RETURN;
  }

  for (i = 0; i < adi_arg->nstp - 1; i++) {
    adi_arg->ifr[i] =
        swp_f_zr(adi_arg->date[i], adi_arg->date[i + 1], adi_arg->yc);
  }

  /*	Fill limit conditions (product) */

  adi_arg->is_event = (int *)calloc(adi_arg->nstp, sizeof(int));
  adi_arg->void_prm = (void **)calloc(adi_arg->nstp, sizeof(void *));

  if (!adi_arg->is_event || !adi_arg->void_prm) {
    err = "Memory allocation error (4) in am_fill_adi_arg";
    goto FREE_RETURN;
  }

  j = am->num_calls - 1;

  for (i = adi_arg->nstp - 1; i >= 0; i--) {
    if (j >= 0 && fabs(adi_arg->date[i] - am->call[j].ex_date) < 1.0e-08) {
      am_prm = malloc(sizeof(am_pay_arg));
      am_prm->und = und;
      am_prm->am = am;

      am_prm->call_idx = j;

      err = am_fill_eval_const_ts(nLambda, pdLambdaValue, pdLambdaTime, und, am,
                                  j, &(am_prm->eval_const));

      adi_arg->is_event[i] = 1;
      adi_arg->void_prm[i] = (void *)am_prm;

      j--;
    } else {
      adi_arg->is_event[i] = 0;
      adi_arg->void_prm[i] = NULL;
    }
  }

FREE_RETURN:

  if (err) {
    am_free_adi_arg(adi_arg);
  }

  return err;
}

Err am_fill_adi_arg(
    AM_UND und, AM_STR am,
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    /*	Required number of steps */
    int req_stp, int req_stpx, AM_ADI_ARG adi_arg) {
  Err err = NULL;
  AM_PAY_ARG am_prm;
  int i, j;

  /*	Initialise */
  adi_arg->time = NULL;
  adi_arg->date = NULL;

  adi_arg->void_prm = NULL;
  adi_arg->is_event = NULL;
  adi_arg->ifr = NULL;

  /*	Compute time steps */

  /*	Copy event dates */

  adi_arg->nstp = am->num_calls;
  adi_arg->nstpx = req_stpx;

  adi_arg->time = (double *)calloc(adi_arg->nstp, sizeof(double));
  if (!adi_arg->time) {
    err = "Memory allocation error (1) in am_fill_adi_arg";
    goto FREE_RETURN;
  }
  for (i = 0; i < adi_arg->nstp; i++) {
    adi_arg->time[i] = am->call[i].ex_time;
  }

  /*	Fill time vector */

  /*	Add today if required */
  if (adi_arg->time[0] < -EPS) {
    err = "Past event date in am_fill_adi_arg";
    goto FREE_RETURN;
  }
  if (adi_arg->time[0] > EPS) {
    num_f_add_number(&(adi_arg->nstp), &(adi_arg->time), 0.0);
    num_f_sort_vector(adi_arg->nstp, adi_arg->time);
    num_f_unique_vector(&(adi_arg->nstp), adi_arg->time);
  }

  /*	If only one event today  , add empty event */
  if (adi_arg->nstp == 1) {
    num_f_add_number(&(adi_arg->nstp), (&adi_arg->time), 1.0);
  }

  /*	Fill the vector */
  num_f_fill_vector_newalgo(&(adi_arg->nstp), &(adi_arg->time), req_stp);

  /*	Make dates */
  adi_arg->date = (double *)calloc(adi_arg->nstp, sizeof(double));
  if (!adi_arg->date) {
    err = "Memory allocation error (2) in am_fill_adi_arg";
    goto FREE_RETURN;
  }

  for (i = 0; i < adi_arg->nstp; i++) {
    adi_arg->date[i] = und->today + DAYS_IN_YEAR * adi_arg->time[i];

    if (i > 0 && adi_arg->date[i] - adi_arg->date[i - 1] >= 1) {
      adi_arg->date[i] = (long)(adi_arg->date[i] + 1.0e-08);
      adi_arg->time[i] = YEARS_IN_DAY * (adi_arg->date[i] - und->today);
    }
  }

  /*	Sigma */

  adi_arg->sig_time = und->sigma_time;
  adi_arg->sig1 = und->sigma;
  adi_arg->nb_sig = und->sigma_n;

  /*	Lambdas and ratio */
  adi_arg->lam = und->lambda;
  adi_arg->alpha = und->alpha;
  adi_arg->gamma = und->gamma;
  adi_arg->rho = und->rho;

  /*	Spot fx and yield curves */

  strcpy(adi_arg->yc, und->yc);

  /*	Fill distributions */

  adi_arg->ifr = (double *)calloc(adi_arg->nstp, sizeof(double));

  if (!adi_arg->ifr) {
    err = "Memory allocation error (3) in am_fill_adi_arg";
    goto FREE_RETURN;
  }

  for (i = 0; i < adi_arg->nstp - 1; i++) {
    adi_arg->ifr[i] =
        swp_f_zr(adi_arg->date[i], adi_arg->date[i + 1], adi_arg->yc);
  }

  /*	Fill limit conditions (product) */

  adi_arg->is_event = (int *)calloc(adi_arg->nstp, sizeof(int));
  adi_arg->void_prm = (void **)calloc(adi_arg->nstp, sizeof(void *));

  if (!adi_arg->is_event || !adi_arg->void_prm) {
    err = "Memory allocation error (4) in am_fill_adi_arg";
    goto FREE_RETURN;
  }

  j = am->num_calls - 1;

  for (i = adi_arg->nstp - 1; i >= 0; i--) {
    if (j >= 0 && fabs(adi_arg->date[i] - am->call[j].ex_date) < 1.0e-08) {
      am_prm = malloc(sizeof(am_pay_arg));
      am_prm->und = und;
      am_prm->am = am;

      am_prm->call_idx = j;

      err = am_fill_eval_const(und, am, j, &(am_prm->eval_const));

      adi_arg->is_event[i] = 1;
      adi_arg->void_prm[i] = (void *)am_prm;

      j--;
    } else {
      adi_arg->is_event[i] = 0;
      adi_arg->void_prm[i] = NULL;
    }
  }

FREE_RETURN:

  if (err) {
    am_free_adi_arg(adi_arg);
  }

  return err;
}

Err am_free_adi_arg(AM_ADI_ARG adi_arg) {
  int i;
  AM_PAY_ARG am_prm;

  if (adi_arg->time)
    free(adi_arg->time);
  if (adi_arg->date)
    free(adi_arg->date);
  if (adi_arg->ifr)
    free(adi_arg->ifr);
  if (adi_arg->is_event)
    free(adi_arg->is_event);

  if (adi_arg->void_prm) {
    for (i = 0; i < adi_arg->nstp; i++) {
      am_prm = (AM_PAY_ARG)(adi_arg->void_prm[i]);
      free(am_prm);
    }

    free(adi_arg->void_prm);
  }

  adi_arg->time = NULL;
  adi_arg->date = NULL;
  adi_arg->sig1 = NULL;
  adi_arg->sig_time = NULL;
  adi_arg->void_prm = NULL;
  adi_arg->is_event = NULL;

  return NULL;
}

/*	Main function to be called in order to fill and check all structures */
/*	==================================================================== */

Err am_fill_check_all_struct(
    /*	Today's date */
    long today, long theoEndDate,

    /*	The underlying */
    int use_calib, /*	0: use lgm2dund  , 1: calibrate */

    /*		if calib */
    char *yc,         /*	yc */
    char *vc,         /*	vc (only if calib) */
    char *ref,        /*	ref rate (only if calib) */
    char *swap_freq,  /*	swap freq (only if calib) */
    char *swap_basis, /*	swap basis (only if calib) */
    double lambda,    /*	lambda if unique */
    double alpha,     /*	alpha */
    double gamma,     /*	gamma */
    double rho,       /*	rho */
    /*	End of calib params */

    /*		if no calilb */
    char *lgm2dund,

    /*	The structure */

    /*		funding */
    char *fund_ref, double *fund_not, int fund_ncpn, long *fund_fix,
    long *fund_start, long *fund_end, long *fund_pay, char **fund_basis,
    double *fund_spr, double *fund_mrg,

    /*		cf */
    double *fix_not, int fix_ncpn, long *fix_start, long *fix_end,
    long *fix_pay, char **fix_basis, double *fix_rate, double *fix_fee,

    /*		calls */
    int ncall, int pay_rec, /*	0: rec pd  , 1: pay pd */
    long *ex_date, long *set_date, double *fee,

    /*	Numerical params */
    int req_stp, int req_stpx,

    /*	Calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),

    double mintime, double mininterval,

    int notperiod, int one2F, int use_jump, double max_var_jump,
    int strike_type, int european_model,

    Err (*get_correl)(char *vol_curve_name, double start_date, double end_date,
                      double strike, double *vol),
    char *CorrelName,

    double max_std_short, int fix_lambda, /*	0: calib lambda to cap  , 1: fix
                                             lambda calib to diagonal */
    int one_f_equi,                       /*	1F equivalent flag:
                                                                          if set to 1  , then 2F
                                             lambda will calibrate                       to the cap priced within calibrated
                                             1F                       with the given lambda */
    int skip_last, /*	If 1  , the last option is disregarded
                                                   and the forward volatility is
                      flat from option n-1 */
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I  , 1: E */
    int eod_ex_flag,  /*	0: I  , 1: E */

    /*	Results */
    AM_STR am, AM_UND und,
    int *call_feat, /*	0: No callable feature to be valued
                            1: Callable feature to be valued through adi */
    AM_ADI_ARG adi_arg) {
  Err err = NULL;
  char swap_freq_loc[256];
  char swap_basis_loc[256];
  int intfreq;

  /*	Initialisation */
  am->fund_leg = NULL;
  am->fix_leg = NULL;
  am->call = NULL;
  am->theoEndDate = theoEndDate;

  und->sigma_n = 0;
  und->sigma_date = NULL;
  und->sigma_time = NULL;
  und->sigma = NULL;

  strcpy(swap_basis_loc, fix_basis[0]);

  intfreq =
      (int)(12 * (fix_pay[fix_ncpn - 1] - fix_start[fix_ncpn - 1]) / 365.0 +
            0.5);
  if (intfreq == 1) {
    strcpy(swap_freq_loc, "M");
  } else if (intfreq == 3) {
    strcpy(swap_freq_loc, "Q");
  } else if (intfreq == 6) {
    strcpy(swap_freq_loc, "S");
  } else if (intfreq == 12) {
    strcpy(swap_freq_loc, "A");
  }

  cpd_init_calib_inst_data(&(und->inst_data));
  und->has_inst_data = 0;
  und->today = today;

  adi_arg->time = NULL;
  adi_arg->date = NULL;
  adi_arg->sig1 = NULL;
  adi_arg->sig_time = NULL;
  adi_arg->void_prm = NULL;
  adi_arg->is_event = NULL;
  adi_arg->ifr = NULL;

  /*	Funding leg */

  am->fund_leg = (AM_FUND_LEG)malloc(sizeof(am_fund_leg));
  if (!am->fund_leg) {
    err = "Memory allocation error (1) in am_fill_check_all_struct";
    goto FREE_RETURN;
  }

  err = am_fill_fund_leg(und->today, eod_fix_flag, fund_not, fund_ncpn,
                         fund_fix, fund_start, fund_end, fund_pay, fund_basis,
                         fund_spr, fund_mrg, am->fund_leg);
  if (err) {
    goto FREE_RETURN;
  }

  /*	fix leg */

  am->fix_leg = (AM_FIX_LEG)malloc(sizeof(am_fix_leg));
  if (!am->fix_leg) {
    err = "Memory allocation error (2) in am_fill_check_all_struct";
    goto FREE_RETURN;
  }

  err = am_fill_fix_leg(und->today, eod_fix_flag, fix_not, fix_ncpn,
                        //			fix_fix  ,
                        fix_start, fix_end, fix_pay, fix_basis, fix_rate,
                        fix_fee, am->fix_leg);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Calls */

  if (ncall > 0 && ex_date[ncall - 1] >= und->today + eod_ex_flag) {
    err = am_fill_calls(und->today, eod_ex_flag, ncall, pay_rec, ex_date,
                        set_date, fee, am);
    if (err) {
      goto FREE_RETURN;
    }
  } else {
    am->num_calls = 0;
    am->call = NULL;
  }

  /*	Underlying */
  if (am->num_calls > 0 && am->call) {
    if (use_calib) {
      err = am_calib_und_new(
          today, eod_ex_flag,

          yc, vc, ref, swap_freq, swap_basis,

          fund_ref, swap_freq_loc, swap_basis_loc, lambda, alpha, gamma, rho,

          mintime, mininterval,

          notperiod, max_std_short, one2F, fix_lambda,

          one_f_equi,

          skip_last,

          use_jump, max_var_jump, strike_type, european_model,

          get_correl, CorrelName,
          /*
                                                          GetVolForBadr  ,
                                                          cVolType  ,
          */
          am,

          get_cash_vol,

          und);

    } else {
      err = am_fill_und(lgm2dund, vc, ref, swap_freq, swap_basis, und);
    }
  }

  if (err) {
    goto FREE_RETURN;
  }

  /*	Tree */

  if (am->num_calls > 0 && am->call) {
    err = am_fill_adi_arg(und, am, get_cash_vol, req_stp, req_stpx, adi_arg);
    if (err) {
      goto FREE_RETURN;
    }

    *call_feat = 1;
  } else {
    *call_feat = 0;
  }

FREE_RETURN:

  if (err) {
    am_free_all_struct(am, und, *call_feat, adi_arg);
  }

  return err;
}

Err am_fill_check_all_struct_ts(
    /*	Today's date */
    long today, long theoEndDate,

    /*	The underlying */
    int use_calib, /*	0: use lgm2dund  , 1: calibrate */

    /*		if calib */
    char *yc,         /*	yc */
    char *vc,         /*	vc (only if calib) */
    char *ref,        /*	ref rate (only if calib) */
    char *swap_freq,  /*	swap freq (only if calib) */
    char *swap_basis, /*	swap basis (only if calib) */

    int nlambda, double *pdlambda_time, double *pdlambda,

    double alpha, /*	alpha */
    double gamma, /*	gamma */
    double rho,   /*	rho */
    /*	End of calib params */

    /*		if no calilb */
    char *lgm2dund,

    /*	The structure */

    /*		funding */
    char *fund_ref, double *fund_not, int fund_ncpn, long *fund_fix,
    long *fund_start, long *fund_end, long *fund_pay, char **fund_basis,
    double *fund_spr, double *fund_mrg,

    /*		cf */
    double *fix_not, int fix_ncpn, long *fix_start, long *fix_end,
    long *fix_pay, char **fix_basis, double *fix_rate, double *fix_fee,

    /*		calls */
    int ncall, int pay_rec, /*	0: rec pd  , 1: pay pd */
    long *ex_date, long *set_date, double *fee,

    /*	Numerical params */
    int req_stp, int req_stpx,

    /*	Calib params */
    Err (*get_cash_vol)(/*	function to get IR cash vol from the markets */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),

    double mintime, double mininterval,

    int notperiod, int one2F, int use_jump, double max_var_jump,
    int strike_type, int european_model,

    Err (*get_correl)(char *vol_curve_name, double start_date, double end_date,
                      double strike, double *vol),
    char *CorrelName,

    double max_std_short, int fix_lambda, /*	0: calib lambda to cap  , 1: fix
                                             lambda calib to diagonal */
    int one_f_equi,                       /*	1F equivalent flag:
                                                                          if set to 1  , then 2F
                                             lambda will calibrate                       to the cap priced within calibrated
                                             1F                       with the given lambda */
    int skip_last, /*	If 1  , the last option is disregarded
                                                   and the forward volatility is
                      flat from option n-1 */
    /*	EOD Flags */
    int eod_fix_flag, /*	0: I  , 1: E */
    int eod_ex_flag,  /*	0: I  , 1: E */

    /*	Results */
    AM_STR am, AM_UND und,
    int *call_feat, /*	0: No callable feature to be valued
                            1: Callable feature to be valued through adi */
    AM_ADI_ARG adi_arg) {
  Err err = NULL;
  char swap_freq_loc[256];
  char swap_basis_loc[256];
  int intfreq;

  /*	Initialisation */
  am->fund_leg = NULL;
  am->fix_leg = NULL;
  am->call = NULL;
  am->theoEndDate = theoEndDate;

  und->sigma_n = 0;
  und->sigma_date = NULL;
  und->sigma_time = NULL;
  und->sigma = NULL;

  strcpy(swap_basis_loc, fix_basis[0]);

  intfreq =
      (int)(12 * (fix_pay[fix_ncpn - 1] - fix_start[fix_ncpn - 1]) / 365.0 +
            0.5);
  if (intfreq == 1) {
    strcpy(swap_freq_loc, "M");
  } else if (intfreq == 3) {
    strcpy(swap_freq_loc, "Q");
  } else if (intfreq == 6) {
    strcpy(swap_freq_loc, "S");
  } else if (intfreq == 12) {
    strcpy(swap_freq_loc, "A");
  }

  cpd_init_calib_inst_data(&(und->inst_data));
  und->has_inst_data = 0;
  und->today = today;

  adi_arg->time = NULL;
  adi_arg->date = NULL;
  adi_arg->sig1 = NULL;
  adi_arg->sig_time = NULL;
  adi_arg->void_prm = NULL;
  adi_arg->is_event = NULL;
  adi_arg->ifr = NULL;

  /*	Funding leg */

  am->fund_leg = (AM_FUND_LEG)malloc(sizeof(am_fund_leg));
  if (!am->fund_leg) {
    err = "Memory allocation error (1) in am_fill_check_all_struct";
    goto FREE_RETURN;
  }

  err = am_fill_fund_leg(und->today, eod_fix_flag, fund_not, fund_ncpn,
                         fund_fix, fund_start, fund_end, fund_pay, fund_basis,
                         fund_spr, fund_mrg, am->fund_leg);
  if (err) {
    goto FREE_RETURN;
  }

  /*	fix leg */

  am->fix_leg = (AM_FIX_LEG)malloc(sizeof(am_fix_leg));
  if (!am->fix_leg) {
    err = "Memory allocation error (2) in am_fill_check_all_struct";
    goto FREE_RETURN;
  }

  err = am_fill_fix_leg(und->today, eod_fix_flag, fix_not, fix_ncpn,
                        //			fix_fix  ,
                        fix_start, fix_end, fix_pay, fix_basis, fix_rate,
                        fix_fee, am->fix_leg);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Calls */

  if (ncall > 0 && ex_date[ncall - 1] >= und->today + eod_ex_flag) {
    err = am_fill_calls(und->today, eod_ex_flag, ncall, pay_rec, ex_date,
                        set_date, fee, am);
    if (err) {
      goto FREE_RETURN;
    }
  } else {
    am->num_calls = 0;
    am->call = NULL;
  }

  /*	Underlying */
  if (am->num_calls > 0 && am->call) {
    if (use_calib) {
      err = am_calib_und_new_ts(
          today, eod_ex_flag,

          yc, vc, ref, swap_freq, swap_basis,

          fund_ref, swap_freq_loc, swap_basis_loc, nlambda, pdlambda_time,
          pdlambda, alpha, gamma, rho,

          mintime, mininterval,

          notperiod, max_std_short, one2F, fix_lambda,

          one_f_equi,

          skip_last,

          use_jump, max_var_jump, strike_type, european_model,

          get_correl, CorrelName,
          /*
                                                          GetVolForBadr  ,
                                                          cVolType  ,
          */
          am,

          get_cash_vol,

          und);

    } else {
      err = am_fill_und(lgm2dund, vc, ref, swap_freq, swap_basis, und);
    }
  }

  if (err) {
    goto FREE_RETURN;
  }

  /*	Tree */

  if (am->num_calls > 0 && am->call) {
    err = am_fill_adi_arg_ts(nlambda, pdlambda, pdlambda_time, und, am,
                             get_cash_vol, req_stp, req_stpx, adi_arg);
    if (err) {
      goto FREE_RETURN;
    }

    *call_feat = 1;
  } else {
    *call_feat = 0;
  }

FREE_RETURN:

  if (err) {
    am_free_all_struct(am, und, *call_feat, adi_arg);
  }

  return err;
}
/*	Free all structures */
Err am_free_all_struct(AM_STR am, AM_UND und, int call_feat,
                       AM_ADI_ARG adi_arg) {
  am_free_und(und);
  am_free_calls(am);

  if (am->fund_leg) {
    am_free_fund_leg(am->fund_leg);
    free(am->fund_leg);
  }

  if (am->fix_leg) {
    am_free_fix_leg(am->fix_leg);
    free(am->fix_leg);
  }

  am_free_adi_arg(adi_arg);

  return NULL;
}

static Err LGM2FCalcPhi(int nstept, double *time, double *sigma,
                        double *sigma_time, int nb_sigma, double lambda,
                        double alpha, double gamma, double rho, double *phi1,
                        double *phi2, double *phi12) {
  double t1, t2, ta, tb, dt;
  double alpha2, gamma2;
  double lamtemp, sig1, lam1;

  LGMSVSolFunc Phi1_, Phi2_, Phi12_;

  LGMSVSolFunc *Phi1 = &Phi1_, *Phi2 = &Phi2_, *Phi12 = &Phi12_;

  double vector[20];
  int i, j, nb_sig_minus1;

  double *ts_time = NULL, *sig = NULL, *lam = NULL;

  int nb_ts;

  Err err = NULL;

  memset(vector, 0, 20);

  ts_time = calloc(nb_sigma, sizeof(double));

  if (!ts_time) {
    err = "Mermory allocation failure in LGM2FExpectations2";
    goto FREE_RETURN;
  }

  memcpy(ts_time, sigma_time, nb_sigma * sizeof(double));
  nb_ts = nb_sigma;

  sig = calloc(nb_ts, sizeof(double));
  lam = calloc(nb_ts, sizeof(double));

  if (!sig || !lam) {
    err = "Mermory allocation failure in LGM2FExpectations2";
    goto FREE_RETURN;
  }

  for (i = 0; i < nb_ts; i++) {
    sig[i] = sigma[i];
    lam[i] = lambda;
  }

  nb_sig_minus1 = nb_ts - 1;
  alpha2 = alpha * alpha;
  gamma2 = 2.0 * gamma;

  /* Initialisation */
  phi1[0] = 0.0;
  phi2[0] = 0.0;
  phi12[0] = 0.0;

  Phi1->bIsft1 = 1;
  Phi1->b = 0.0;
  Phi1->bIsgt1 = 1;
  Phi1->c = 0.0;
  Phi1->bIsht1 = 1;

  Phi2->bIsft1 = 1;
  Phi2->b = 0.0;
  Phi2->bIsgt1 = 1;
  Phi2->c = 0.0;
  Phi2->bIsht1 = 1;

  Phi12->bIsft1 = 1;
  Phi12->b = 0.0;
  Phi12->bIsgt1 = 1;
  Phi12->c = 0.0;
  Phi12->bIsht1 = 1;

  lam1 = lam[0];
  sig1 = sig[0];

  Phi1->a = sig1 * sig1;
  Phi1->dLambda = 2.0 * lam1;
  Phi2->a = alpha2 * Phi1->a;
  Phi2->dLambda = Phi1->dLambda + gamma2;
  Phi12->a = alpha * Phi1->a;
  Phi12->dLambda = Phi1->dLambda + gamma;

  t1 = 0.0;
  j = 0;

  for (i = 1; i < nstept; i++) {
    t2 = time[i];
    dt = t2 - t1;

    ta = t1;
    tb = ts_time[j];

    lamtemp = 0.0;

    Phi1->dXt1 = phi1[i - 1];
    Phi2->dXt1 = phi2[i - 1];
    Phi12->dXt1 = phi12[i - 1];

    while (tb < t2 && j < nb_sig_minus1) {
      dt = (tb - ta);

      /* Calculation at intermedary time tb */
      LGMSVFuncValue2(Phi1, dt, vector, 0, &Phi1->dXt1);
      LGMSVFuncValue2(Phi2, dt, vector, 0, &Phi2->dXt1);
      LGMSVFuncValue2(Phi12, dt, vector, 0, &Phi12->dXt1);

      lamtemp += lam1 * dt;

      j++;
      ta = tb;
      tb = ts_time[j];
      lam1 = lam[j];
      sig1 = sig[j];

      /* Update the model parameters */
      Phi1->a = sig1 * sig1;
      Phi1->dLambda = 2.0 * lam1;
      Phi2->a = alpha2 * Phi1->a;
      Phi2->dLambda = Phi1->dLambda + gamma2;
      Phi12->a = alpha * Phi1->a;
      Phi12->dLambda = Phi1->dLambda + gamma;
    }

    dt = (t2 - ta);

    LGMSVFuncValue2(Phi1, dt, vector, 0, &phi1[i]);
    LGMSVFuncValue2(Phi2, dt, vector, 0, &phi2[i]);
    LGMSVFuncValue2(Phi12, dt, vector, 0, &phi12[i]);

    lamtemp += lam1 * dt;

    t1 = t2;
  }

FREE_RETURN:

  if (sig)
    free(sig);
  if (lam)
    free(lam);
  if (ts_time)
    free(ts_time);

  return err;
}

Err am_calc_mdl_iv_fwd(AM_STR am, AM_UND und, AM_ADI_ARG adi_arg,

                       int num_hermite,
                       /*	Result */
                       double *premium) {

  int i, j, k, index;
  double *time = NULL, *x = NULL, *w = NULL, *phi1 = NULL, *phi2 = NULL,
         *phi12 = NULL, *r1 = NULL, **r2 = NULL, *r3 = NULL, ***pv = NULL;

  double std1, std2, std3;
  double coef, sum, sum_part, df, fee;

  AM_CALL call;

  Err err = NULL;

  /* memory allocation */

  time = (double *)calloc(am->num_calls + 1, sizeof(double));

  x = dvector(1, num_hermite);
  w = dvector(1, num_hermite);

  phi1 = (double *)calloc(am->num_calls + 1, sizeof(double));
  phi2 = (double *)calloc(am->num_calls + 1, sizeof(double));
  phi12 = (double *)calloc(am->num_calls + 1, sizeof(double));

  r1 = dvector(0, num_hermite - 1);
  r2 = dmatrix(0, num_hermite - 1, 0, num_hermite - 1);
  r3 = dvector(0, num_hermite - 1);

  pv = f3tensor(0, num_hermite - 1, 0, num_hermite - 1, 0, 0);

  if (!time || !x || !w || !phi1 || !phi2 || !phi12 || !r1 || !r2 || !r3 ||
      !pv) {
    err = "Memory allocation failure in am_calc_mdl_iv_fwd";
    goto FREE_RETURN;
  }

  /* Hermite calculation */
  err = HermiteStandard(x, w, num_hermite);

  if (err) {
    goto FREE_RETURN;
  }

  time[0] = 0.0;
  for (i = 0; i < am->num_calls; i++) {
    call = am->call + i;
    time[i + 1] = call->ex_time;
  }

  /* Phi Calculation */
  err = LGM2FCalcPhi(am->num_calls + 1, time, und->sigma, und->sigma_time,
                     und->sigma_n, und->lambda, und->alpha, und->gamma,
                     und->rho, phi1, phi2, phi12);

  if (err) {
    goto FREE_RETURN;
  }

  index = 0;

  for (i = 0; i < am->num_calls; i++) {
    /* first remove the fee */
    call = am->call + i;
    fee = call->fee;
    call->fee = 0.0;

    std1 = sqrt(phi1[i + 1]);
    std2 = sqrt(phi2[i + 1]);

    coef = -und->rho * phi12[i + 1] / phi1[i + 1];

    /* r3 = r2 + coef * r1 */
    std3 = sqrt(phi2[i + 1] + coef * coef * phi1[i + 1] +
                2.0 * und->rho * coef * std1 * std2);

    /* construct the grid */
    for (j = 0; j < num_hermite; j++) {
      r1[j] = std1 * x[j + 1];
      r3[j] = std3 * x[j + 1];
    }

    /* reconstruction */
    for (j = 0; j < num_hermite; j++) {
      for (k = 0; k < num_hermite; k++) {
        r2[j][k] = (r3[k] - coef * r1[j]);
        pv[j][k][0] = -BIG_BIG_NUMBER;
      }
    }

    /* launch evaluation */

    while (!((adi_arg->void_prm)[index]) && index < 10000) {
      index++;
    }

    err = am_payoff_4_lgm2f_adi(
        call->ex_date, call->ex_time, (adi_arg->void_prm)[index], und->yc,
        &und->lambda, NULL, 1, und->gamma, und->rho, phi1[i + 1], phi2[i + 1],
        phi12[i + 1], 0, num_hermite - 1, 0, num_hermite - 1, r1, r2, 1, pv);

    index++;

    /* put the fee back */
    call->fee = fee;

    /* convolution */

    sum = 0.0;
    sum_part = 0.0;

    for (j = 0; j < num_hermite; j++) {
      sum_part = 0.0;
      for (k = 0; k < num_hermite; k++) {
        sum_part += w[k + 1] * pv[j][k][0];
      }

      sum = sum + w[j + 1] * sum_part;
    }

    df = swp_f_df(und->today, call->ex_date, und->yc);

    sum *= df;

    premium[i] = sum;
  }

FREE_RETURN:

  if (time)
    free(time);

  if (x)
    free_dvector(x, 1, num_hermite);
  if (w)
    free_dvector(w, 1, num_hermite);

  if (phi1)
    free(phi1);
  if (phi2)
    free(phi2);
  if (phi12)
    free(phi12);

  if (r1)
    free_dvector(r1, 0, num_hermite - 1);
  if (r2)
    free_dmatrix(r2, 0, num_hermite - 1, 0, num_hermite - 1);
  if (r3)
    free_dvector(r3, 0, num_hermite - 1);

  if (pv)
    free_f3tensor(pv, 0, num_hermite - 1, 0, num_hermite - 1, 0, 0);

  return err;
}

Err am_payoff_4_lgm2f_adi(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    void *yc, double *lam, double *ts_time, int nb_ts, double gamma, double rho,
    double phi1, double phi2, double phi12,
    /* Nodes data */
    int l1, int u1, int l2, int u2, double *r1, double **r2, int nprod,
    /* Vector of results to be updated */
    double ***prod_val) {
  AM_PAY_ARG am_arg;
  AM_STR am;
  AM_CALL call;
  AM_UND und;
  AM_FUND_LEG fund_leg;
  AM_FUND_CPN fund_cpn;
  AM_FIX_LEG fix_leg;
  AM_FIX_CPN fix_cpn;
  AM_EVAL_CONST eval_const;
  int call_idx;

  int i, j, l;

  double R1, R2, *r2i;

  double fund_leg_pv, df, coupon, fix_leg_pv, iv;

  int num_fund_cpn, num_fix_cpn;

  double fee;

  Err err = NULL;

  /*	Get the event */

  am_arg = (AM_PAY_ARG)func_parm;
  eval_const = (AM_EVAL_CONST)(&(am_arg->eval_const));
  am = am_arg->am;
  call_idx = am_arg->call_idx;
  call = am->call + call_idx;
  und = (AM_UND)(am_arg->und);
  fund_leg = am->fund_leg;
  fix_leg = am->fix_leg;

  num_fund_cpn = call->num_fund_cpn;
  num_fix_cpn = call->num_fix_cpn;

  /*	Eval payoff
          ----------- */

  for (i = l1; i <= u1; i++) {
    R1 = r1[i];
    r2i = r2[i];
    for (j = l2; j <= u2; j++) {
      R2 = r2i[j];

      /*	Calculate discount factors	*/
      for (l = 0; l < eval_const->num_df; l++) {
        eval_const->df[l] =
            exp(-eval_const->df_alpha[l] - eval_const->df_beta[l] * R1 -
                eval_const->df_gamma[l] * R2);
      }

      /*	PV of funding leg */
      fund_leg_pv = 0.0;

      /*	DF initial */

      fund_leg_pv += eval_const->df[eval_const->start_idx] *
                     fund_leg->notional[call->fund_idx];

      /*	Fund Spread Coupons */
      for (l = 0; l < num_fund_cpn; l++) {
        fund_cpn = am->fund_leg->cpn + call->fund_idx + l;
        fund_leg_pv += eval_const->df[eval_const->fund_idx[l]] * fund_cpn->cpn;
      }

      /*	PV of FIX leg */
      fix_leg_pv = 0.0;

      /*	Coupons */
      for (l = 0; l < num_fix_cpn; l++) {
        /*	Coupon access */
        fix_cpn = fix_leg->cpn + call->fix_idx + l;

        /*	Discount */
        df = eval_const->df[eval_const->fix_disc_idx[l]];

        /*	Midat */
        coupon = fix_cpn->cpn;

        /*	Coupon pv */
        fix_leg_pv += df * coupon;
      }

      /*	Intrinsic value */
      if (call->pay_rec == 0) {
        iv =
            fix_leg_pv - fund_leg_pv -
            eval_const->df[eval_const->start_idx] * fix_leg->fee[call->fix_idx];
      } else {
        iv =
            fund_leg_pv - fix_leg_pv +
            eval_const->df[eval_const->start_idx] * fix_leg->fee[call->fix_idx];
      }

      /*	Fee */
      if (fabs(call->fee) > 1.0e-08) {
        fee = eval_const->df[eval_const->fee_idx] * call->fee;
      } else {
        fee = 0.0;
      }

      /*	Process max */

      if (prod_val[i][j][0] < iv - fee) {
        prod_val[i][j][0] = iv - fee;
      }
    }
  }

  /*	End of payoff valuation
          ----------------------- */

  return err;
}

/*	Main pricing function */

/*	Launch the pde */
Err am_launch_adi(AM_STR am, AM_UND und, AM_ADI_ARG adi_arg,
                  /*	Result */
                  double *prem) {
  double temp_val[2];
  Err err;
  int one_lam = 1;
  double fixed_lam;

  /* launch the corresponding ADI */

  fixed_lam = adi_arg->lam;

  err = lgm2f_adi(adi_arg->nstp, adi_arg->time, adi_arg->date, adi_arg->nstpx,
                  fixed_lam, adi_arg->sig_time, adi_arg->sig1, adi_arg->nb_sig,
                  adi_arg->alpha, adi_arg->gamma, adi_arg->rho,
                  adi_arg->void_prm, adi_arg->is_event, adi_arg->ifr,
                  adi_arg->yc, am_payoff_4_lgm2f_adi, 1, (double *)temp_val);

  *prem = temp_val[0];

  return err;
}

/*	Launch the pde */
Err am_launch_adi_ts(double *pdLambdaValue, double *pdLambdaTime,
                     int nLambdaSize, AM_STR am, AM_UND und, AM_ADI_ARG adi_arg,
                     /*	Result */
                     double *prem) {
  double temp_val[2];
  Err err;
  int one_lam = 1;
  // double fixed_lam;
  int ndisc_method = 1; /// 0 linear discretization  , 1. normal discretization

  /* launch the corresponding ADI */

  // fixed_lam = adi_arg->lam;

  err = lgm2fTau_adi(
      adi_arg->nstp, adi_arg->time, adi_arg->date,

      adi_arg->nstpx, ndisc_method,

      adi_arg->sig1, adi_arg->sig_time, adi_arg->nb_sig,

      pdLambdaValue,
      pdLambdaTime, //// NB: calendar time  , not year fraction !!!!!!!
      nLambdaSize,

      adi_arg->alpha, adi_arg->gamma, adi_arg->rho, adi_arg->void_prm,
      adi_arg->is_event, adi_arg->ifr, adi_arg->yc,
      am_payoff_4_lgm2f_adi, // payoff_lgm2fTau_pde  ,//am_payoff_4_lgm2f_adi  ,
      1, (double *)temp_val);

  *prem = temp_val[0];

  return err;
}
