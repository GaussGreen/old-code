/* ==========================================================
        FILENAME :			OTCcaller.cxx

        PURPOSE:			Calculates One Time Callables using a
   smile on the Fx by a copula based method

        AUTHOR:				J. Dinh

        DATE:				25-FEB-2004
   ========================================================== */

#include "OTCcaller.h"
#include "CPDVol.h"
#include "CPDcaller.h"
#include "Fx3FBetaDLMCalculations.h"
#include "Fx3FBetaDLMCalibration.h"
#include "Fx3FBetaDLMUtil.h"
#include "GenericMidatAutocal.h"
#include "MCEBOptimisation.h"
#include "OTCutils.h"

#include "math.h"
#include "num_h_interp.h"
#include "opfnctns.h"
#include "opsabrgenericinterp.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srtaccess.h"

/*	Caller for callable power duals */
/*	------------------------------- */

#define POS_VAL(X) ((X) > 0 ? (X) : 0)

#define CALL_VAL(FWD, STRIKE, D, S) ((FWD)*norm((D) + (S)) - (STRIKE)*norm((D)))

#define PUT_VAL(FWD, STRIKE, D, S)                                             \
  (-(FWD)*norm(-(D) - (S)) + (STRIKE)*norm(-(D)))

Err otc_caller(
    /*	Today's date */
    long today,
    /*	The underlying */
    int use_calib, /*	0: use fx3dund        , 1: calibrate */
    /*		if calib */
    double fx_spot,                   /*	2bd fwd */
    long fx_spot_date, int dom_calib, /*	Calibrate domestic underlying */
    char *dom_und,        /*	If no        , domestic underlying to be used */
    char *dom_yc,         /*	Domestic yc */
    char *dom_vc,         /*	Domestic vc (only if calib) */
    char *dom_ref,        /*	Domestic ref rate (only if calib) */
    char *dom_swap_freq,  /*	Domestic swap freq (only if calib) */
    char *dom_swap_basis, /*	Domestic swap basis (only if calib) */
    double dom_lam,       /*	Domestic lambda */
    int for_calib,        /*	Same for foreign */
    char *for_und, char *for_yc, char *for_vc, char *for_ref,
    char *for_swap_freq, char *for_swap_basis, double for_lam,

    double min_fact, /*	Maximum down jump on variance */
    double max_fact, /*	Maximum up jump on variance */
    int use_jumps,   /*	Allow vol term structure to jump */

    double *corr_times, double *correl_dom_for, /*	Correlations */
    double *correl_dom_fx, double *correl_for_fx,
    double **cumulative_corr_matrix, // Dom-For || Dom-Fx || For-Fx
    long corr_n_times, CPDBETADLMPARAMS cpd_dlm_params,
    Err (*get_ir_cash_vol)(/*	Function to get IR cash vol from the markets */
                           char *vol_curve_name, double start_date,
                           double end_date, double cash_strike, int zero,
                           char *ref_rate_name, double *vol, double *power),
    /*	Fx vol from the market */
    long *fx_mkt_vol_date, double *fx_mkt_vol, int num_fx_mkt_vol,
    /*	Fx SABR parameters from the market */
    double *fx_mkt_smile_alpha, double *fx_mkt_smile_beta,
    double *fx_mkt_smile_rho, double *fx_mkt_smile_pi,
    int use_sabr, /*	0: no smile adj-t        , 1: only for underlying , 2:
               fees adj-t for the call/KO */
    int smile_spec_type, //	0: lognormal vol + SABR params        , 1:
                         //sigma-beta
                         //+ SABR params        , 2: BMM (not yet supported)
    /*		if no calilb */
    char *fx3dund,
    /*	The structure */
    long start_date, /*	Date at which initial notional exchange occurs */
    /*		funding */
    double fund_not, /*	Notional */
    int fund_ccy,    /*	0: domestic        , 1: foreign 2: other */
    char *
        fund_ccy_yc, /*	If different from domestic or foreign (fund_ccy = 2) */
    double fx_fund_dom, /*	If different from domestic or foreign (fund_ccy
                     = 2) 2 bd fwd */
    long fx_fund_dom_spot_date, int fund_ncpn, /*	Number of coupons */
    long *fund_fix,                            /*	Fixing dates */
    long *fund_start,                          /*	Start dates */
    long *fund_pay,                            /*	Pay dates */
    char **fund_basis,                         /*	Basis */
    double *fund_spr,                          /*	Forward spreads */
    double *fund_mrg,                          /*	Margins */
    double *fund_fix_cpn, /*	Past coupon fixing if relevant        ,
                                                    includes spr        , but
                       not mrg        , cvg and notional */
    /*		pd */
    double pd_not,   /*	Notional */
    int pd_ncpn,     /*	Number of coupons */
    long *pd_fix,    /*	Fx fixing dates */
    long *pd_start,  /*	Start dates */
    long *pd_pay,    /*	Pay dates */
    char **pd_basis, /*	Basis */
    double
        *pd_alpha, /*	Coupon = alpha + beta * fx [capped        , floored] */
    double *pd_beta, int *pd_floored, double *pd_floor, int *pd_capped,
    double *pd_cap,
    /*		pd interp coupon specification */
    int *pd_nfxpts,         /*	Number of coupon interpolation points */
    double **pd_fxpts,      /*	Coupon interpolation points */
    double **pd_cpn_at_pts, /*	Coupon values at interpolation points
                         (interpolation is linear) */
    int *pd_lin_xtrpl_l,    /*	0 = flat extrapolation to the left        , 1 =
                         linear    w first 2 pts slope */
    int *pd_lin_xtrpl_r,    /*	0 = flat extrapolation to the right        , 1 =
                         linear w last 2 pts slope */
    double *pd_fix_fx,      /*	Past Fx fixing if relevant */
    /*		pd not refund */
    long *pd_not_ref_fix,    //	fx fixing dates FOR EACH CALL DATE + in the end
                             // if relevant
    double pd_not_ref_alpha, /*	Final notional on PD leg */
    double pd_not_ref_beta, int pd_not_ref_floored, double pd_not_ref_floor,
    int pd_not_ref_capped, double pd_not_ref_cap,
    /*		pd not interp specification */
    int *pd_not_ref_nfxpts, /*	Number of notional interpolation points FOR EACH
                         CALL DATE + in the end */
    double **pd_not_ref_fxpts, double **pd_not_ref_cpn_at_pts,
    int *pd_not_ref_lin_xtrpl_l, /*	0 = flat extrapolation to the left , 1
                              = linear w first 2 pts slope */
    int *pd_not_ref_lin_xtrpl_r, /*	0 = flat extrapolation to the right , 1
                              = linear w last 2 pts slope */
    double *pd_not_ref_fix_fx,   /*	fx fixings FOR EACH CALL DATE + in the
                              end if relevant */
    /*		calls */
    int *call_type,  /*	0: call        , 1: KO */
    int ncall,       /*	Number of calls */
    int pay_rec,     /*	0: rec pd        , 1: pay pd */
    long *ex_date,   /*	Call dates */
    long *set_date,  /*	Settlement dates */
    double *barrier, /*	in case of a pure KO or a Callable KO */
    int *bar_type,   /*	0: up and in        , 1: down and in */
    double *fees,    /*  fees if deal is called in domestic currency */
    int TARN_Do,     //	1: Prices a Target note Powerdual
    int use_GMA,     //	0: Does not output the multicallable price | 1: GMA1 |
                     // 1: GMA2 | 3: GMA1 & GMA2
    /*	Numerical params */
    long req_stp,      /*	Number of time steps in the tree */
    long req_pth,      /*	Number of paths in the MC */
    double bar_smooth, /*	Smoothing factor for barriers */
    int do_pecs,       /*	Do PECS in the MC */
    int forcetree,     /*	If equal to 1 then the valuation is done in a tree
                        */
    int do_optim,      /*	If equal to 1 then the call are replaced by optimal KO
                        */
    int force_optim, /*	If equal to 1 then all call will be replaced by optimal
                  KO	*/
    int fx_bound,    /*	If equal to 1 then optimisation on the Fx        , on
                  the IV    otherwise	*/
    int use_bound,   /*	If equal to 1 then prices the call as UO on the Fx using
                  a provided boundary	*/
    int do_infos,    /*	infos on callable right */
    /*	EOD Fixing Flag */
    int eod_fix_flag, /*	0: I        , 1: E */
    /*	EOD Payment Flag */
    int eod_pay_flag, /*	0: I        , 1: E */
    /*	EOD Exercise Flag */
    int eod_ex_flag, /*	0: I        , 1: E */
    /*	Exercised flag */
    int exercised,    /*	Flag */
    long ex_date_ex,  /*	Date when exercised */
    long ex_date_set, /*	Corresponding settlement date */
    /*  Parameters */
    int nSimul, int Do_pecs, int nPoints, int nStd, double CummulPrecision,
    int CummulLinear, double std, int nIter, int PayoffFunction,
    double fwdSmileVisuNstd,
    /* otc */
    int fwdVolMethod, /*	0=Sliding 3F vol; 1=Converging 3F vol; 2=Sliding
                   Cvg Sbeta; 3=Cvg Cvg Sbeta */
    int smileOtc,     /*	use smile parameters for the OTC */
    int smileFee,     /*	use smile parameters for the OTC Fees */
    int otc,          /*	which call to keep */
    int smileModel,   /*	0=SSL; 1=Log Mix; 2=Beta Mix */
    double BMpi,      /*	if smileModel = 1 or 2 then BMpi is the probability of
                   the state 1 */
    /* Correl */
    int firstIndex, int secondIndex, int firstLong, int secondLong,
    long correlTstar, char *CPDsigma, char *CPDalpha, char *CPDbeta,
    char *CPDrho,
    /* Change the Funding for speed */
    int FundingSpeedUp,
    /* Change the calibration strategy */
    int LongShort,
    /* Do not calculate all OTC*/
    int nStart, int oneOutOfN,
    /*Fast MC*/
    int FMC_do, double FMC_precision, int FMC_min_paths,
    /*	Results */
    double **Results) {

  /* Declaration */
  double time1, time2;
  int i, otcStart = otc, otcEnd = otc;
  cpd_str cpd_, *cpd = &cpd_;
  cpd_und und_, *und = &und_;
  cpd_tree_arg tree_arg_, *tree_arg = &tree_arg_;
  cpd_mc_arg mc_arg_, *mc_arg = &mc_arg_;

  int call_feat;

  int free_struct = 0;
  int skip_fill = 0;

  int for_fund;
  long fund_start_date;
  double eq_final_ex, eq_init_ex;

  double thirdCcyfundNot;
  char Copyfund_ccy_yc[255];

  /* For the Copula */
  double *mergeTimes = NULL, *sigDom = NULL, *sigFor = NULL, *sigFx = NULL,
         *corDF = NULL, *corDFx = NULL, *corFFx = NULL;
  int mergeNtimes;

  /* for coupon interpolation */
  int max_interp_pts = 2, alloc_pd_ncpn = pd_ncpn, alloc_ncall = ncall;
  int ifst_pt, ilst_pt;
  int *pd_num_strikes = NULL, *pd_not_num_strikes = NULL;
  double *pd_wcst = NULL, *pd_wspot = NULL, **pd_strikes = NULL,
         **pd_weights = NULL;
  double *pd_not_wcst = NULL, *pd_not_wspot = NULL, **pd_not_strikes = NULL,
         **pd_not_weights = NULL;
  int use_cpn_opt_str = (pd_nfxpts != NULL),
      use_not_opt_str = (pd_not_ref_nfxpts != NULL);
  double fwd, temp, smile_std;

  otcpd_params *OTC_params = NULL;
  otcpd_precalc *OTC_precalc = NULL;

  // VOL MKT STRUCTURE
  SMILE_VOL_MARKET smile_mkt = NULL;
  SMILE_PARAMETERS smile_params = NULL;
  double *fwd_at_vol_dates = NULL;

  Err err = NULL;

  if ((use_cpn_opt_str || use_not_opt_str) && exercised) {
    err = serror("Exercised flag is temporarily not supported in case of "
                 "option string coupon or redemption");
    goto FREE_RETURN;
  }

  // Calculate max_interp_pts:

  if (use_cpn_opt_str)
    for (i = 0; i < pd_ncpn; i++)
      if (pd_nfxpts[i] > max_interp_pts)
        max_interp_pts = pd_nfxpts[i];
  if (use_not_opt_str)
    for (i = 0; i <= ncall; i++)
      if (pd_not_ref_nfxpts[i] > max_interp_pts)
        max_interp_pts = pd_not_ref_nfxpts[i];

  // Transform smile specification if not default:
  if (smile_spec_type == 1 ||
      smile_spec_type == 3) // SABR Beta -> SABR Log || BMM Beta -> BMM Log
  {
    for (i = 0; i < num_fx_mkt_vol; i++) {
      fwd = fx_spot * swp_f_df(fx_spot_date, fx_mkt_vol_date[i], for_yc) /
            swp_f_df(fx_spot_date, fx_mkt_vol_date[i], dom_yc);

      temp = (add_unit(fx_mkt_vol_date[i], -2, SRT_BDAY, MODIFIED_SUCCEEDING) -
              today) *
             YEARS_IN_DAY;

      if (smile_spec_type == 1) {
        err = srt_f_optsarbvol(fwd, fwd, temp, fx_mkt_vol[i],
                               fx_mkt_smile_alpha[i], fx_mkt_smile_beta[i],
                               fx_mkt_smile_rho[i], SRT_BETAVOL, SRT_LOGNORMAL,
                               &smile_std);
        if (err)
          goto FREE_RETURN;
      } else {
        err = srt_f_optbmmvol(fwd, fwd, temp, fx_mkt_vol[i],
                              fx_mkt_smile_alpha[i], fx_mkt_smile_beta[i],
                              fx_mkt_smile_rho[i], fx_mkt_smile_pi[i],
                              SRT_BETAVOL, SRT_LOGNORMAL, &smile_std);
      }

      fx_mkt_vol[i] = smile_std;
    }
    smile_spec_type--; // 1->0 || 3->2
  }

  /* ==================================
  If tierce funding
  ================================== */
  if (fund_ccy == 2) {
    /* save the initial fund not in the third ccy
       and the third Ccy yield curve name
       used only when exercised flag is TRUE	   */
    thirdCcyfundNot = fund_not;
    strcpy(Copyfund_ccy_yc, fund_ccy_yc);

    /* Convert into domestic */
    fund_ccy = 0;
    for_fund = 1;
    err = convert_funding_to_domestic(
        today, start_date, eod_fix_flag, eod_pay_flag, fx_fund_dom,
        fx_fund_dom_spot_date, pd_not, dom_yc, fund_ncpn, fund_fix, fund_start,
        fund_pay, fund_basis, fund_ccy_yc, &fund_not, fund_spr, fund_mrg,
        fund_fix_cpn, &fund_start_date, &eq_final_ex, &eq_init_ex);
    if (err) {
      return err;
    }
  } else {
    for_fund = 0;
  }

  // Allocate memory for strings of calls:

  pd_num_strikes = (int *)calloc(alloc_pd_ncpn, sizeof(int));
  pd_wcst = (double *)calloc(alloc_pd_ncpn, sizeof(double));
  pd_wspot = (double *)calloc(alloc_pd_ncpn, sizeof(double));
  pd_strikes = dmatrix(0, alloc_pd_ncpn - 1, 0, max_interp_pts - 1);
  pd_weights = dmatrix(0, alloc_pd_ncpn - 1, 0, max_interp_pts - 1);
  pd_not_num_strikes = (int *)calloc(alloc_ncall + 1, sizeof(int));
  pd_not_wcst = (double *)calloc(alloc_ncall + 1, sizeof(double));
  pd_not_wspot = (double *)calloc(alloc_ncall + 1, sizeof(double));
  pd_not_strikes = dmatrix(0, alloc_ncall, 0, max_interp_pts - 1);
  pd_not_weights = dmatrix(0, alloc_ncall, 0, max_interp_pts - 1);

  if (!pd_num_strikes || !pd_wcst || !pd_wspot || !pd_strikes || !pd_weights ||
      !pd_not_num_strikes || !pd_not_wcst || !pd_not_wspot || !pd_not_strikes ||
      !pd_not_weights) {
    err = serror("Memory failure in cpd_autocal");
    goto FREE_RETURN;
  }

  // Now convert all types of coupon specifications into const + fwd + sum of
  // calls:

  // Start with coupons:
  if (use_cpn_opt_str) {
    for (i = 0; i < pd_ncpn; i++) {
      err = transform_interp_coupon(pd_nfxpts[i], pd_fxpts[i], pd_cpn_at_pts[i],
                                    pd_lin_xtrpl_l[i], pd_lin_xtrpl_r[i],
                                    &pd_wcst[i], &pd_wspot[i], &ifst_pt,
                                    &ilst_pt, pd_weights[i]);
      if (err)
        goto FREE_RETURN;

      pd_num_strikes[i] = ilst_pt - ifst_pt + 1;
      memcpy(pd_strikes[i], pd_fxpts[i] + ifst_pt,
             pd_num_strikes[i] * sizeof(double));
      if (ifst_pt > 0)
        memmove(pd_weights[i], pd_weights[i] + ifst_pt,
                pd_num_strikes[i] * sizeof(double));
    }
  } else // in this case transform the old alpha-beta specification into a new
         // one:
  {
    for (i = 0; i < pd_ncpn; i++) {
      transform_oldspec_into_newspec(
          pd_alpha[i], pd_beta[i], pd_floored[i], pd_floor[i], pd_capped[i],
          pd_cap[i], &pd_num_strikes[i], &pd_wcst[i], &pd_wspot[i],
          pd_strikes[i], pd_weights[i]);
    }
  }

  // Now convert the redemption:
  if (use_not_opt_str) {
    for (i = 0; i <= ncall; i++) {
      err = transform_interp_coupon(
          pd_not_ref_nfxpts[i], pd_not_ref_fxpts[i], pd_not_ref_cpn_at_pts[i],
          pd_not_ref_lin_xtrpl_l[i], pd_not_ref_lin_xtrpl_r[i], &pd_not_wcst[i],
          &pd_not_wspot[i], &ifst_pt, &ilst_pt, pd_not_weights[i]);
      if (err)
        goto FREE_RETURN;

      pd_not_num_strikes[i] = ilst_pt - ifst_pt + 1;
      memcpy(pd_not_strikes[i], pd_not_ref_fxpts[i] + ifst_pt,
             pd_not_num_strikes[i] * sizeof(double));
      if (ifst_pt > 0)
        memmove(pd_not_weights[i], pd_not_weights[i] + ifst_pt,
                pd_not_num_strikes[i] * sizeof(double));
    }
  } else {
    for (i = 0; i < ncall;
         i++) // Only domestic redemption for exercised deal in the old case
    {
      pd_not_num_strikes[i] = 0;
      pd_not_wcst[i] = 1.0;
      pd_not_wspot[i] = 0.0;
    }
    transform_oldspec_into_newspec(
        pd_not_ref_alpha, pd_not_ref_beta, pd_not_ref_floored, pd_not_ref_floor,
        pd_not_ref_capped, pd_not_ref_cap, &pd_not_num_strikes[ncall],
        &pd_not_wcst[ncall], &pd_not_wspot[ncall], pd_not_strikes[ncall],
        pd_not_weights[ncall]);
  }

  /* ==================================
  Initialise structures
  ================================== */

  //=========================================
  //
  //		CREATE A VOL MARKET STRUCTURE
  //
  //=========================================
  smile_mkt = calloc(1, sizeof(smile_vol_market));
  smile_params = calloc(1, sizeof(smile_parameters));
  fwd_at_vol_dates = calloc(num_fx_mkt_vol, sizeof(double));

  if (!smile_mkt || !smile_params || !fwd_at_vol_dates) {
    err = "CPDcaller: Memory allocation failure";
    if (err)
      goto FREE_RETURN;
  }

  // memory allocation
  err = cpd_alloc_smile_vol_market(num_fx_mkt_vol, smile_mkt);
  if (err)
    goto FREE_RETURN;

  // create the vector of forwards
  for (i = 0; i < num_fx_mkt_vol; i++)
    fwd_at_vol_dates[i] = fx_spot *
                          swp_f_df(fx_spot_date, fx_mkt_vol_date[i], for_yc) /
                          swp_f_df(fx_spot_date, fx_mkt_vol_date[i], dom_yc);

  // fill the mkt structure
  err = cpd_fill_smile_vol_market(
      today, smile_spec_type, num_fx_mkt_vol, fwd_at_vol_dates, fx_mkt_vol,
      fx_mkt_smile_alpha, fx_mkt_smile_beta, fx_mkt_smile_rho, fx_mkt_smile_pi,
      NULL, fx_mkt_vol_date, smile_mkt);

  if (err)
    goto FREE_RETURN;

  // free the vect of forwards
  if (fwd_at_vol_dates)
    free(fwd_at_vol_dates);

  //=========================================
  //
  //		INITIALISE CPD STRUCTURES
  //
  //=========================================
  skip_fill = 0;

  err = cpd_fill_check_all_struct(
      today, use_calib, fx_spot, fx_spot_date, dom_calib, dom_und, dom_yc,
      dom_vc, dom_ref, dom_swap_freq, dom_swap_basis, dom_lam, for_calib,
      for_und, for_yc, for_vc, for_ref, for_swap_freq, for_swap_basis, for_lam,
      min_fact, max_fact, use_jumps, corr_times, correl_dom_for, correl_dom_fx,
      correl_for_fx, corr_n_times, cpd_dlm_params, get_ir_cash_vol,
      fx_mkt_vol_date, fx_mkt_vol, num_fx_mkt_vol, fx3dund, fund_not, fund_ccy,
      fund_ncpn, fund_fix, fund_start, fund_pay, fund_basis, fund_spr, fund_mrg,
      pd_not, pd_ncpn, pd_fix, pd_start, pd_pay, pd_basis, pd_alpha, pd_beta,
      pd_floored, pd_floor, pd_capped, pd_cap, use_cpn_opt_str, pd_num_strikes,
      pd_wcst, pd_wspot, pd_strikes, pd_weights, pd_not_ref_fix,
      pd_not_ref_alpha, pd_not_ref_beta, pd_not_ref_floored, pd_not_ref_floor,
      pd_not_ref_capped, pd_not_ref_cap, use_not_opt_str, pd_not_num_strikes,
      pd_not_wcst, pd_not_wspot, pd_not_strikes, pd_not_weights, call_type,
      ncall, pay_rec, ex_date, set_date, barrier, bar_type, fees,
      TARN_Do, // Is TARN?
      req_stp, req_pth, do_pecs, forcetree, do_optim, force_optim, fx_bound,
      use_bound, bar_smooth, eod_fix_flag, eod_ex_flag, cpd, und, &call_feat,
      tree_arg, mc_arg, 0.0, 0.0, 0.0, skip_fill);

  if (err)
    goto FREE_RETURN;

  /* ==================================
  Make changes in the call notional exchange for past call dates
  ================================== */

  /*	Skip calls to be exercised before today */
  i = 0;
  while (i < ncall && ex_date[i] < today + eod_fix_flag)
    i++;

  free_struct = 1;

  /* ==================================
  Allocation
  ================================== */

  // Test if the funding change is admissible
  if (cpd->fund_leg->dom_for) {
    if (cpd->call[0].ex_date >
        cpd->pd_leg->cpn[cpd->call[0].pd_idx].fx_val_date) {
      FundingSpeedUp = 0;
    }
  }

  OTC_params = calloc(1, sizeof(otcpd_params));

  cpd_init_otc_params(OTC_params, cpd, ncall, nSimul, Do_pecs, nPoints, nStd,
                      CummulPrecision, CummulLinear, std, nIter, PayoffFunction,
                      fwdSmileVisuNstd,
                      /* otc */
                      fwdVolMethod, /*	0=Sliding 3F vol; 1=Converging 3F vol;
                                 2=Sliding Cvg Sbeta; 3=Cvg Cvg Sbeta */
                      smileOtc,     /*	use smile parameters for the OTC */
                      smileFee,     /*	use smile parameters for the OTC Fees */
                      smileModel,   /*	0=SSL; 1=BM; 2=HestonD */
                      otc,          /*	0=SSL; 1=Log Mix; 2=Beta Mix */
                      BMpi,
                      /* Correl */
                      firstIndex, secondIndex, firstLong, secondLong,
                      correlTstar, CPDsigma, CPDalpha, CPDbeta, CPDrho,
                      /* Time Dimension */
                      FundingSpeedUp, nStart, oneOutOfN,
                      /*Fast MC*/
                      FMC_do, FMC_precision, FMC_min_paths,
                      /*Which GMA*/
                      use_GMA, LongShort);

  /* ==================================
  Merge the TS
  ================================== */
  err = merge_rates_fx_corr_ts(
      und->sigma_time_rates, und->sigma_dom, und->sigma_n_rates,
      und->sigma_time_rates, und->sigma_for, und->sigma_n_rates,
      und->sigma_time_fx, und->sigma_fx, und->sigma_n_fx, und->corr_times,
      und->correl_dom_for, und->correl_dom_fx, und->correl_for_fx,
      und->corr_n_times, &mergeTimes, &sigDom, &sigFor, &sigFx, &corDF, &corDFx,
      &corFFx, &mergeNtimes);

  if (err)
    goto FREE_RETURN;

  /* ==================================
  Precalculations
  ================================== */
  time1 = clock();
  OTC_precalc = calloc(1, sizeof(otcpd_precalc));
  cpd_init_otc_precalc(OTC_precalc);

  err = cpd_otc_precalc(OTC_precalc, OTC_params, cpd, und, smile_mkt,
                        mergeTimes, mergeNtimes, sigDom, sigFor, sigFx, corDF,
                        corDFx, corFFx);

  if (err)
    goto FREE_RETURN;
  time2 = clock();
  smessage("Phase 0 -Precalculations        , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  /* ==================================
  Launch the pricers
  ================================== */
  if (!cpd->type) {
    err = otc_pricer(
        cpd, und, OTC_params, OTC_precalc, smile_mkt, smile_params, pd_not,
        mergeNtimes, mergeTimes, sigDom, sigFor, sigFx, corDF, corDFx, corFFx,
        start_date, /*	Date at which initial notional exchange occurs */
        cumulative_corr_matrix, 0, NULL,
        /*	Results */
        Results);
  } else {
    err = ko_pricer(
        cpd, und, OTC_params, OTC_precalc, smile_mkt, smile_params, pd_not,
        mergeNtimes, mergeTimes, sigDom, sigFor, sigFx, corDF, corDFx, corFFx,
        start_date, /*	Date at which initial notional exchange occurs */
        /*	Results */
        Results);
  }

FREE_RETURN:

  if (smile_mkt) {
    cpd_free_smile_vol_market(smile_mkt);
    free(smile_mkt);
  }

  if (smile_params)
    free(smile_params);

  if (pd_num_strikes)
    free(pd_num_strikes);
  if (pd_wcst)
    free(pd_wcst);
  if (pd_wspot)
    free(pd_wspot);
  if (pd_strikes)
    free_dmatrix(pd_strikes, 0, alloc_pd_ncpn - 1, 0, max_interp_pts - 1);
  if (pd_weights)
    free_dmatrix(pd_weights, 0, alloc_pd_ncpn - 1, 0, max_interp_pts - 1);
  if (pd_not_num_strikes)
    free(pd_not_num_strikes);
  if (pd_not_wcst)
    free(pd_not_wcst);
  if (pd_not_wspot)
    free(pd_not_wspot);
  if (pd_not_strikes)
    free_dmatrix(pd_not_strikes, 0, alloc_ncall, 0, max_interp_pts - 1);
  if (pd_not_weights)
    free_dmatrix(pd_not_weights, 0, alloc_ncall, 0, max_interp_pts - 1);

  if (OTC_params) {
    cpd_free_otc_params(OTC_params);
    free(OTC_params);
  }

  if (OTC_precalc) {
    cpd_free_otc_precalc(OTC_precalc, cpd->num_calls, cpd->pd_leg->num_cpn,
                         cpd->fund_leg->num_cpn);
    free(OTC_precalc);
  }

  if (free_struct)
    cpd_free_all_struct(cpd, und, tree_arg, mc_arg);

  if (mergeTimes)
    free(mergeTimes);
  if (sigDom)
    free(sigDom);
  if (sigFor)
    free(sigFor);
  if (sigFx)
    free(sigFx);
  if (corDF)
    free(corDF);
  if (corDFx)
    free(corDFx);
  if (corFFx)
    free(corFFx);

  return err;
}

Err FxFwdSmilePrice(
    long today, long spot_date, double spot_fx, long forwardFixDate,
    double forwardFixTime, long forwardValDate, double forwardValTime,
    long fixDate, double fixTime, long valDate, double valTime, char *dom_yc,
    char *for_yc, double *sigma_time_dom, int sigma_n_dom, double *sigma_dom,
    double lda_dom, double *sigma_time_for, int sigma_n_for, double *sigma_for,
    double lda_for, double *sigma_time_fx, double *sigma_fx, int sigma_n_fx,
    double *corr_times, double *correl_dom_for, double *correl_dom_fx,
    double *correl_for_fx, int corr_n_times, double smile_vol,
    double smile_alpha, double smile_beta, double smile_rho, double SigmaFwd,
    double AlphaFwd, double BetaFwd, double RhoFwd, int nSimul, int nPoints,
    int nIter, int nStd, double std, int smileModel, int nStrikes,
    double *Strikes, int isCall, int OutputVol, int StochVolFudge,
    double SwitchLevel, double Proba, double AlphaFudge, double *Results) {
  double forward;
  double *mergeTimes = NULL, *sigDom = NULL, *sigFor = NULL, *sigFx = NULL,
         *corDF = NULL, *corDFx = NULL, *corFFx = NULL;
  int mergeNtimes;

  double **matrix = NULL;

  /* For the measure */
  double tstarTime;
  int tstarDate;

  Err err = NULL;

  /* ==================================
  get the end time and date
  ================================== */
  tstarDate = forwardFixDate;
  tstarTime = forwardFixTime;

  forward = spot_fx * swp_f_df(spot_date, tstarDate, for_yc) /
            swp_f_df(spot_date, tstarDate, dom_yc);

  matrix = dmatrix(0, nSimul - 1, 0, 3 - 1);

  /* ==================================
  Merge the TS
  ================================== */
  err = merge_rates_fx_corr_ts(
      sigma_time_dom, sigma_dom, sigma_n_dom, sigma_time_for, sigma_for,
      sigma_n_for, sigma_time_fx, sigma_fx, sigma_n_fx, corr_times,
      correl_dom_for, correl_dom_fx, correl_for_fx, corr_n_times, &mergeTimes,
      &sigDom, &sigFor, &sigFx, &corDF, &corDFx, &corFFx, &mergeNtimes);

  if (err)
    goto FREE_RETURN;

  err =
      OTCgetCopula(forward, forwardFixTime, tstarTime, mergeTimes, mergeNtimes,
                   sigDom, lda_dom, sigFor, lda_for, sigFx, corDF, corDFx,
                   corFFx, smile_vol, smile_alpha, smile_beta, smile_rho,
                   nSimul, nPoints, nIter, nStd, std, smileModel, matrix);

  if (err)
    goto FREE_RETURN;

  /* ==================================
  Payoff
  ================================== */
  forward = spot_fx * swp_f_df(spot_date, valDate, for_yc) /
            swp_f_df(spot_date, valDate, dom_yc);

  if (!StochVolFudge) {
    err = fwdSABRpayoff(today, forwardFixDate, forwardFixTime, valDate, fixTime,
                        valTime, forward, SigmaFwd, AlphaFwd, BetaFwd, RhoFwd,
                        matrix, nSimul, tstarTime, tstarDate, dom_yc, for_yc,
                        mergeTimes, mergeNtimes, sigDom, lda_dom, sigFor,
                        lda_for, mergeTimes, sigFx, mergeNtimes, mergeTimes,
                        corDF, corDFx, corFFx, mergeNtimes, nStrikes, Strikes,
                        isCall, OutputVol, Results);

    if (err)
      goto FREE_RETURN;
  } else {
    err = fwdSABRpayoffAlphaFudge(
        today, forwardFixDate, forwardFixTime, valDate, fixTime, valTime,
        forward, SigmaFwd, AlphaFwd, BetaFwd, RhoFwd, matrix, nSimul, tstarTime,
        tstarDate, dom_yc, for_yc, mergeTimes, mergeNtimes, sigDom, lda_dom,
        sigFor, lda_for, mergeTimes, sigFx, mergeNtimes, mergeTimes, corDF,
        corDFx, corFFx, mergeNtimes, nStrikes, Strikes, isCall, OutputVol,
        SwitchLevel, Proba, AlphaFudge, Results);

    if (err)
      goto FREE_RETURN;
  }

FREE_RETURN:
  if (matrix)
    free_dmatrix(matrix, 0, nSimul - 1, 0, 3 - 1);
  if (mergeTimes)
    free(mergeTimes);
  if (sigDom)
    free(sigDom);
  if (sigFor)
    free(sigFor);
  if (sigFx)
    free(sigFx);
  if (corDF)
    free(corDF);
  if (corDFx)
    free(corDFx);
  if (corFFx)
    free(corFFx);

  return err;
}

Err otc_pricer(
    CPD_STR cpd, CPD_UND und, otcpd_params *OTC_params,
    otcpd_precalc *OTC_precalc, SMILE_VOL_MARKET smile_mkt,
    SMILE_PARAMETERS smile_params, double pd_not, int mergeNtimes,
    double *mergeTimes, double *sigDom, double *sigFor, double *sigFx,
    double *corDF, double *corDFx, double *corFFx,
    long start_date, /*	Date at which initial notional exchange occurs */
    double **cumulative_corr_matrix,
    int erasing_call_done, //	Don't pass to GMA all the OTC
    int *erased_call_list, //	List of erasedz
    double **Results) {
  /* Declaration */
  int i, j, jj, k, otc, otcStart = OTC_params->OTC, otcEnd = OTC_params->OTC;

  /* For the Copula */
  double time1, time2;
  double stdDom, cumulDomForCor, stdFor, cumulForFxCor, meanFor, cumulDomFxCor;

  double **partialResults = NULL;

  /* For the SSL */
  double error, forward;

  double tstarTime;
  int tstarDate;

  /* =====================
          GMA
  ==================== */
  double *dExeTimes = NULL, *dForward = NULL, *dLongOption = NULL,
         *dShortOption = NULL;

  genmidat_calibparams *calib_params = NULL;
  genmidat_pdeparams *pde_params = NULL;
  genmidat_autocalparams *autocal_params = NULL;
  genmidat_autocalinfos *sInfos = NULL;
  genmidat_model *sModel = NULL;

  Err err = NULL;

  /* ==================================
  gaussian simulation
  ================================== */
  /*	If the payoff == 6        , then we price both smile and flat
          simultaneously and hence we need an additional dimension to the matrix
   */
  OTC_params->COPULAmatrix = dmatrix(0, OTC_params->COPULAnSimul - 1, 0,
                                     3 - 1 + OTC_params->OTCsmileAndFlat);
  OTC_params->COPULAgauss = dmatrix(0, OTC_params->COPULAnSimul - 1, 0, 3 - 1);
  OTC_params->CUMUL = dmatrix(0, 2 - 1, 0, OTC_params->CUMULnPoints - 1);

  if (!OTC_params->COPULAmatrix || !OTC_params->COPULAgauss ||
      !OTC_params->CUMUL) {
    err = "otc_caller: memory allocation error";
    goto FREE_RETURN;
  }

  time1 = clock();

  err = balsam_generation(OTC_params->COPULAnSimul, 3, OTC_params->COPULAgauss);

  if (err)
    goto FREE_RETURN;

  time2 = clock();
  smessage("Phase 1 -BalSam generation        , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);
  time1 = clock();

  if (OTC_params->COPULAdo_pecs) {
    adjust_emp_covar(OTC_params->COPULAgauss, OTC_params->COPULAnSimul, 3);
  }

  time2 = clock();
  smessage("Phase 2 -PECS adjustment        , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  /* ==================================
  LOOP ON ALL THE REQUESTED OTC
  ================================== */
  if (OTC_params->OTCpayoff == 4 || OTC_params->OTCpayoff == 6 ||
      OTC_params->OTCpayoff == 123) {
    otcStart = 0;
    otcEnd = cpd->num_calls - 1;
  }

  for (otc = otcStart; otc <= otcEnd;
       otc++) // TO do Check if the call has been erased
  {
    if (!erasing_call_done || (erasing_call_done && !erased_call_list[otc])) {
      OTC_params->OTC = otc;

      if (otc > cpd->num_calls - 1) {
        err = "OTCcaller: OTC index out of range";
        goto FREE_RETURN;
      }

      /* ==================================
      get the measure time and date
      ================================== */
      tstarDate = cpd->call[otc].ex_date;
      tstarTime = cpd->call[otc].ex_time;

      /* ==================================
      get the 3F Cumulative correlations at the call date
      ================================== */
      // Use the inputed correlation or the 3F correlation
      err = OTCgetMoments(cpd->call[otc].ex_time, und->lda_dom, und->lda_for,
                          tstarTime, mergeTimes, mergeNtimes, sigDom, sigFor,
                          sigFx, corDF, corDFx, corFFx, &stdDom, &meanFor,
                          &stdFor, &cumulDomForCor, &cumulDomFxCor,
                          &cumulForFxCor);

      if (cumulative_corr_matrix) {
        cumulDomForCor = cumulative_corr_matrix[otc][0];
        cumulDomFxCor = cumulative_corr_matrix[otc][1];
        cumulForFxCor = cumulative_corr_matrix[otc][2];
      }

      if (err)
        goto FREE_RETURN;

      err = OTCcorrelateVariables(stdDom, stdFor, cumulDomForCor, cumulDomFxCor,
                                  cumulForFxCor, OTC_params->COPULAnSimul,
                                  OTC_params->COPULAgauss,
                                  OTC_params->COPULAmatrix);

      if (err)
        goto FREE_RETURN;

      if (OTC_params->OTCsmileAndFlat) {
        for (i = 0; i < OTC_params->COPULAnSimul; i++)
          OTC_params->COPULAmatrix[i][3] = OTC_params->COPULAmatrix[i][2];
      }

      /* ==================================
      get the Sabr parameters at the call date
      ================================== */
      forward = und->spot_fx * swp_f_df(und->today, tstarDate, und->for_yc) /
                swp_f_df(und->today, tstarDate, und->dom_yc);

      err = OTCfwdSABR(forward, 0.0, cpd->call[otc].ex_time, tstarTime,
                       smile_mkt, smile_params, mergeTimes, mergeNtimes, sigDom,
                       und->lda_dom, sigFor, und->lda_for, sigFx, corDF, corDFx,
                       corFFx, 1234);

      if (err)
        goto FREE_RETURN;

      /* ==================================
      Calibrate the Fx smile and generate the smiled realisations
      ================================== */
      if (fabs(cpd->call[otc].ex_time) > 1.0e-10) {
        if (OTC_params->OTCsmile || OTC_params->OTCsmileAndFlat) {
          time1 = clock();

          err =
              OTCcalibSmileModel(cpd->call[otc].ex_time, /* Maturity in years */
                                 forward,                /* Forward */
                                 OTC_params, smile_params, &error);

          if (err)
            goto FREE_RETURN;

          /* ==================================
          Get the Cumulative Fx distribution under Qtstar
          ==================================*/
          err = OTCgetCumulative(forward, cpd->call[otc].ex_time, OTC_params);

          if (err)
            goto FREE_RETURN;

          err = GenerateRealisations(
              0.0, meanFor, OTC_params->CUMUL, OTC_params->COPULAnSimul,
              OTC_params->CUMULnPoints, 1, 2, OTC_params->CUMULlinear,
              OTC_params->COPULAmatrix);

          if (err)
            goto FREE_RETURN;

          time2 = clock();
          smessage("Phase 3 -Simulation Finished        , time in sec: %.2f",
                   (double)(time2 - time1) / CLOCKS_PER_SEC);
          time1 = clock();
        }

        /* ==================================
        Generate the lognormal realisations
        ================================== */
        if (!OTC_params->OTCsmile || OTC_params->OTCsmileAndFlat) {
          time1 = clock();

          k = 2 + OTC_params->OTCsmileAndFlat;

          if (!OTC_params->OTCsmileAndFlat) {
            err = GenerateRealisations(
                0.0, meanFor, OTC_params->CUMUL, OTC_params->COPULAnSimul,
                OTC_params->CUMULnPoints, 0, k, OTC_params->CUMULlinear,
                OTC_params->COPULAmatrix);

            if (err)
              goto FREE_RETURN;
          }

          for (i = 0; i < OTC_params->COPULAnSimul; i++) {
            OTC_params->COPULAmatrix[i][k] =
                forward *
                exp(-0.5 * smile_params->sigma * smile_params->sigma *
                        cpd->call[otc].ex_time +
                    OTC_params->COPULAmatrix[i][k] * smile_params->sigma *
                        sqrt(cpd->call[otc].ex_time));
          }

          time2 = clock();
          smessage("Phase 3 -Simulation Finished        , time in sec: %.2f",
                   (double)(time2 - time1) / CLOCKS_PER_SEC);
          time1 = clock();
        }
      } else {
        err = GenerateRealisations(
            0.0, meanFor, OTC_params->CUMUL, OTC_params->COPULAnSimul,
            OTC_params->CUMULnPoints, 0, 2, OTC_params->CUMULlinear,
            OTC_params->COPULAmatrix);

        for (i = 0; i < OTC_params->COPULAnSimul; i++)
          OTC_params->COPULAmatrix[i][2] = forward;

        if (OTC_params->OTCsmileAndFlat) {
          for (i = 0; i < OTC_params->COPULAnSimul; i++)
            OTC_params->COPULAmatrix[i][3] = forward;
        }
      }

      /* ==================================
      =====================================

                              PAYOFF

      =====================================
      ================================== */
      if (OTC_params->OTCpayoff == 0) {
        err = OTCfwdSABRCalib(
            cpd, und, OTC_precalc->pd_fwdsmile_vol[OTC_params->OTC],
            OTC_precalc->pd_fwdsmile_alpha[OTC_params->OTC],
            OTC_precalc->pd_fwdsmile_beta[OTC_params->OTC],
            OTC_precalc->pd_fwdsmile_rho[OTC_params->OTC], OTC_params->OTC,
            OTC_params->COPULAnSimul, OTC_params->COPULAmatrix, tstarTime,
            tstarDate, OTC_params->OTCfwdSmileVisuStd,
            OTC_precalc->FxAdj[OTC_params->OTC], Results);

        if (err)
          goto FREE_RETURN;
      } else if (OTC_params->OTCpayoff == 1) {
        err = "Wrong Payoff Function ";
        goto FREE_RETURN;
      } else if (OTC_params->OTCpayoff == 2) {
        err = "Wrong Payoff Function ";
        goto FREE_RETURN;
      } else if (OTC_params->OTCpayoff == 3 || OTC_params->OTCpayoff == 4) {
        /*	OUTPUT

                lg fund_leg		pl fund_leg		st fund_leg
           lg fund_fee		st fund_fee
                lg pd_leg		pl pd_leg		st pd_leg
           lg pd_fee		st pd_fee lg otc			pl otc
           st otc lg std			pl std			st std
        */
        partialResults = dmatrix(0, 3, 0, 4);
        if (!partialResults) {
          err = "otc_caller: memory allocation error";
          goto FREE_RETURN;
        }

        err = OTCPDpayoffFeeCombinedSmileAndFlat(
            cpd, und, OTC_params, OTC_precalc, smile_mkt, smile_params, pd_not,
            OTC_params->COPULAmatrix, tstarTime, tstarDate, partialResults);

        if (err)
          goto FREE_RETURN;

        if (OTC_params->OTCpayoff == 3) {
          Results[0][0] = partialResults[0][0];
          Results[0][1] = partialResults[0][1];
          Results[0][2] = partialResults[0][2];
          Results[0][3] = partialResults[0][3];
          Results[0][4] = partialResults[0][4];

          Results[1][0] = partialResults[1][0];
          Results[1][1] = partialResults[1][1];
          Results[1][2] = partialResults[1][2];
          Results[1][3] = partialResults[1][3];
          Results[1][4] = partialResults[1][4];

          Results[2][0] = partialResults[2][0];
          Results[2][1] = partialResults[2][1];
          Results[2][2] = partialResults[2][2];

          Results[3][0] = partialResults[3][0];
          Results[3][1] = partialResults[3][1];
          Results[3][2] = partialResults[3][2];
        } else {

          if (!cpd->call[OTC_params->OTC].pay_rec) {
            Results[otc + 1][0] = partialResults[1][0] - partialResults[0][0];
            Results[otc + 1][4] = partialResults[1][3] - partialResults[0][3];
            Results[otc + 1][5] = partialResults[1][4] - partialResults[0][4];
          } else {
            Results[otc + 1][0] = -partialResults[1][0] + partialResults[0][0];
            Results[otc + 1][4] = -partialResults[1][3] + partialResults[0][3];
            Results[otc + 1][5] = -partialResults[1][4] + partialResults[0][4];
          }

          Results[otc + 1][1] = partialResults[2][0];
          Results[otc + 1][2] = partialResults[2][1];
          Results[otc + 1][3] = partialResults[2][2];
          Results[otc + 1][6] = tstarTime;
        }

        if (partialResults) {
          free_dmatrix(partialResults, 0, 3, 0, 4);
          partialResults = NULL;
        }
      } else if (OTC_params->OTCpayoff == 5) {
        err = OTCPDcorrel(
            cpd, und, pd_not, OTC_precalc->pd_fwdsmile_vol[OTC_params->OTC],
            OTC_precalc->pd_fwdsmile_alpha[OTC_params->OTC],
            OTC_precalc->pd_fwdsmile_beta[OTC_params->OTC],
            OTC_precalc->pd_fwdsmile_rho[OTC_params->OTC],
            OTC_params->OTCsmileFee, OTC_params->OTC, OTC_params->COPULAnSimul,
            OTC_params->COPULAmatrix, tstarTime, tstarDate,
            OTC_precalc->FxAdj[OTC_params->OTC], OTC_precalc->FeeFxAdj,
            OTC_params->CORRELfirstIndex, OTC_params->CORRELsecondIndex,
            OTC_params->CORRELfirstLong, OTC_params->CORRELsecondLong,
            OTC_params->CORRELtStar, Results);

        if (err)
          goto FREE_RETURN;
      } else if (OTC_params->OTCpayoff == 6) {
        err = OTCPDpayoffFeeCombinedSmileAndFlat(
            cpd, und, OTC_params, OTC_precalc, smile_mkt, smile_params, pd_not,
            OTC_params->COPULAmatrix, tstarTime, tstarDate, Results);

        if (err)
          goto FREE_RETURN;
      } else if (OTC_params->OTCpayoff == 123) {
        err = "Wrong Payoff Function ";

        if (err)
          goto FREE_RETURN;

      } else if (OTC_params->OTCpayoff == 1234) {
        forward = und->spot_fx *
                  swp_f_df(und->today, cpd->call[otc].set_date, und->for_yc) /
                  swp_f_df(und->today, cpd->call[otc].set_date, und->dom_yc);

        Results[0][0] = 0.0;

        for (j = 1; j < 9; j++) {
          Results[j][0] =
              forward *
              exp(OTC_params->CUMULnStd * ((double)j - 4) / 9.0 *
                  OTC_params->SSLvolUp * sqrt(cpd->call[otc].set_time));
        }

        for (j = 0; j < 9; j++) {
          err = SSLprice(forward, cpd->call[otc].set_time, Results[j][0],
                         OTC_params->SSLvolUp, OTC_params->SSLvolDown,
                         OTC_params->SSLshift, 1, &(Results[j][2]));

          if (err)
            goto FREE_RETURN;

          err =
              srt_f_optsarbvol(forward, Results[j][0], cpd->call[otc].set_time,
                               smile_params->sigma, smile_params->alpha,
                               smile_params->beta, smile_params->rho,
                               SRT_LOGNORMAL, SRT_LOGNORMAL, &(Results[j][3]));

          if (err)
            goto FREE_RETURN;

          Results[j][3] =
              srt_f_optblksch(forward, Results[j][0], Results[j][3],
                              cpd->call[otc].set_time, 1.0, SRT_CALL, PREMIUM);
        }

        err = OTCtestPayoff(und, OTC_params->COPULAmatrix,
                            OTC_params->COPULAnSimul, tstarTime, tstarDate, 9,
                            Results);

        if (err)
          goto FREE_RETURN;
      } else {
        err = "otc_caller: Invalid PayoffFunction";
        goto FREE_RETURN;
      }

      time2 = clock();
      smessage("Phase 4 -Payoff evaluation Finished        , time in sec: %.2f",
               (double)(time2 - time1) / CLOCKS_PER_SEC);
    }
  }

  /*===================================
          IF THE PAYOFF FUNCTION IS #4 (ALL EUROPEANS)
          CALL GMA TO GET THE MIDAT PRICE FROM THE EUROPEANS
  ===================================*/
  if (OTC_params->OTCpayoff == 4 || OTC_params->OTCpayoff == 6) {
    calib_params =
        (genmidat_calibparams *)calloc(1, sizeof(genmidat_calibparams));
    pde_params = (genmidat_pdeparams *)calloc(1, sizeof(genmidat_pdeparams));
    autocal_params =
        (genmidat_autocalparams *)calloc(1, sizeof(genmidat_autocalparams));
    sInfos = (genmidat_autocalinfos *)calloc(1, sizeof(genmidat_autocalinfos));
    sModel = (genmidat_model *)calloc(1, sizeof(genmidat_model));

    dExeTimes = (double *)calloc(cpd->num_calls, sizeof(double));
    dForward = (double *)calloc(cpd->num_calls, sizeof(double));
    dLongOption = (double *)calloc(cpd->num_calls, sizeof(double));
    dShortOption = (double *)calloc(cpd->num_calls, sizeof(double));

    if (!calib_params || !pde_params || !autocal_params || !sModel || !sInfos ||
        !dExeTimes || !dForward || !dLongOption || !dShortOption) {
      err = "Memory allocation faillure in otc_caller";
      goto FREE_RETURN;
    }
  }

  k = 0;
  if (erasing_call_done) {
    for (i = 0; i < cpd->num_calls; i++) {
      if (erased_call_list[i])
        k++;
    }

    if (k == cpd->num_calls) // no calls left
    {
      Results[0][0] = 0.0;
      Results[0][1] = 0.0;
      Results[0][2] = 0.0;
      Results[0][3] = 0.0;

      goto FREE_RETURN;
    }
  }

  if (OTC_params->OTCpayoff == 4) {
    j = 0;
    for (i = 0; i < cpd->num_calls; i++) {
      if (!erasing_call_done || (erasing_call_done && !erased_call_list[i])) {
        dExeTimes[j] = Results[i + 1][6];
        dForward[j] = Results[i + 1][0];
        dLongOption[j] = Results[i + 1][1];
        dShortOption[j] = Results[i + 1][2];
        j++;
      }
    }

    genmidat_init_autocalinfos(sInfos);
    /* Create the model */
    err = genmidat_alloc_model(cpd->num_calls - k, 1, sModel);

    // Initialisation
    genmidat_set_autocalparams_default(autocal_params);
    genmidat_set_calibparams_default(calib_params);
    genmidat_set_pdeparams_default(pde_params);

    calib_params->iShortIsCap = 1;
    pde_params->lNbTime = 300;
    pde_params->lNbX = 300;

    for (i = 0; i < cpd->num_calls - k; i++) {
      dForward[i] /= sModel->dNumeraire;
    }
    if (err)
      goto FREE_RETURN;

    err = genmidat_init_model(0.0, 0.0, 0.0, dExeTimes, dForward, NULL, NULL,
                              NULL, NULL, sModel->dStartBeta2,
                              sModel->dNumeraire, sModel);

    if (err)
      goto FREE_RETURN;

    err = genmidat_allocate_autocalinfos(cpd->num_calls - k, autocal_params,
                                         sInfos);

    if (err)
      goto FREE_RETURN;

    /* GMA 1 */
    if (OTC_params->use_GMA == 1 || OTC_params->use_GMA == 3) {
      err = GenericMidatAutocal(cpd->num_calls - k, dLongOption, dShortOption,
                                NULL, sModel, autocal_params, calib_params,
                                pde_params, &(Results[0][0]), sInfos);

      if (err)
        goto FREE_RETURN;
    }

    /* GMA 2 */
    if (OTC_params->use_GMA == 2 || OTC_params->use_GMA == 3) {
      calib_params->iShortIsCap = 0;
      j = 0;
      for (i = 0; i < cpd->num_calls; i++) {
        if (!erasing_call_done || (erasing_call_done && !erased_call_list[i])) {
          dShortOption[j] = Results[i + 1][3];
          j++;
        }
      }

      err = GenericMidatAutocal(cpd->num_calls - k, dLongOption, dShortOption,
                                NULL, sModel, autocal_params, calib_params,
                                pde_params, &(Results[0][1]), sInfos);

      if (err)
        goto FREE_RETURN;
    }
  } else if (OTC_params->OTCpayoff == 6) {
    j = 0;
    for (i = 0; i < cpd->num_calls; i++) {
      if (!erasing_call_done || (erasing_call_done && !erased_call_list[i])) {
        dExeTimes[j] = Results[i + 1][0];
        dForward[j] = Results[i + 1][1];
        dLongOption[j] = Results[i + 1][2];
        dShortOption[j] = Results[i + 1][3];
        j++;
      }
    }

    genmidat_init_autocalinfos(sInfos);
    /* Create the model */
    err = genmidat_alloc_model(cpd->num_calls - k, 1, sModel);

    // Initialisation
    genmidat_set_autocalparams_default(autocal_params);
    genmidat_set_calibparams_default(calib_params);
    genmidat_set_pdeparams_default(pde_params);

    calib_params->iShortIsCap = 1;
    pde_params->lNbTime = 300;
    pde_params->lNbX = 300;

    for (i = 0; i < cpd->num_calls - k; i++) {
      dForward[i] /= sModel->dNumeraire;
    }
    if (err)
      goto FREE_RETURN;

    err = genmidat_init_model(0.0, 0.0, 0.0, dExeTimes, dForward, NULL, NULL,
                              NULL, NULL, sModel->dStartBeta2,
                              sModel->dNumeraire, sModel);

    if (err)
      goto FREE_RETURN;

    err = genmidat_allocate_autocalinfos(cpd->num_calls - k, autocal_params,
                                         sInfos);

    if (err)
      goto FREE_RETURN;

    /* if all OTC at the end are worth 0        , then we remove them */
    jj = 0;
    while (fabs(dLongOption[cpd->num_calls - 1 - jj - k]) < 1.0e-20 &&
           jj + k < cpd->num_calls &&
           (!erasing_call_done ||
            (erasing_call_done && !erased_call_list[cpd->num_calls - 1 - jj])))
      jj++;

    if (jj == cpd->num_calls) {
      Results[0][0] = 0.0; // GMA1 with smile
      Results[0][1] = 0.0; // GMA2 with smile
    } else {
      /* GMA 1 */
      if (OTC_params->use_GMA == 1 || OTC_params->use_GMA == 3) {
        err = GenericMidatAutocal(cpd->num_calls - jj - k, dLongOption,
                                  dShortOption, NULL, sModel, autocal_params,
                                  calib_params, pde_params,
                                  &(Results[0][0]), // GMA1 with smile
                                  sInfos);

        if (err)
          goto FREE_RETURN;
      }
      /* GMA 2 */
      if (OTC_params->use_GMA == 2 || OTC_params->use_GMA == 3) {
        calib_params->iShortIsCap = 0;
        j = 0;
        for (i = 0; i < cpd->num_calls; i++) {
          if (!erasing_call_done ||
              (erasing_call_done && !erased_call_list[i])) {
            dShortOption[j] = Results[i + 1][4];
            j++;
          }
        }

        err = GenericMidatAutocal(cpd->num_calls - jj - k, dLongOption,
                                  dShortOption, NULL, sModel, autocal_params,
                                  calib_params, pde_params,
                                  &(Results[0][1]), // GMA2 with smile
                                  sInfos);

        if (err)
          goto FREE_RETURN;
      }
    }

    j = 0;
    for (i = 0; i < cpd->num_calls; i++) {
      if (!erasing_call_done || (erasing_call_done && !erased_call_list[i])) {
        dForward[j] = Results[i + 1][5];
        dLongOption[j] = Results[i + 1][6];
        dShortOption[j] = Results[i + 1][7];
        j++;
      }
    }

    err = genmidat_init_model(0.0, 0.0, 0.0, dExeTimes, dForward, NULL, NULL,
                              NULL, NULL, sModel->dStartBeta2,
                              sModel->dNumeraire, sModel);

    if (err)
      goto FREE_RETURN;

    /* if all OTC at the end are worth 0        , then we remove them */
    jj = 0;
    while (fabs(dLongOption[cpd->num_calls - 1 - jj - k]) < 1.0e-20 &&
           jj + k < cpd->num_calls &&
           (!erasing_call_done ||
            (erasing_call_done && !erased_call_list[cpd->num_calls - 1 - jj])))
      jj++;

    if (jj == cpd->num_calls) {
      Results[0][2] = 0.0; // GMA1 with 3F
      Results[0][3] = 0.0; // GMA2 with 3F
    } else {
      /* GMA 1 */
      if (OTC_params->use_GMA == 1 || OTC_params->use_GMA == 3) {
        calib_params->iShortIsCap = 1;
        err = GenericMidatAutocal(cpd->num_calls - jj - k, dLongOption,
                                  dShortOption, NULL, sModel, autocal_params,
                                  calib_params, pde_params,
                                  &(Results[0][2]), // GMA1 with 3F
                                  sInfos);

        if (err)
          goto FREE_RETURN;
      }
      /* GMA 2 */
      if (OTC_params->use_GMA == 2 || OTC_params->use_GMA == 3) {
        calib_params->iShortIsCap = 0;
        j = 0;
        for (i = 0; i < cpd->num_calls; i++) {
          if (!erasing_call_done ||
              (erasing_call_done && !erased_call_list[i])) {
            dShortOption[j] = Results[i + 1][8];
            j++;
          }
        }

        err = GenericMidatAutocal(cpd->num_calls - jj - k, dLongOption,
                                  dShortOption, NULL, sModel, autocal_params,
                                  calib_params, pde_params,
                                  &(Results[0][3]), // GMA2 with 3F
                                  sInfos);

        if (err)
          goto FREE_RETURN;
      }
    }
  }

FREE_RETURN:

  if (dExeTimes)
    free(dExeTimes);
  if (dForward)
    free(dForward);
  if (dLongOption)
    free(dLongOption);
  if (dShortOption)
    free(dShortOption);

  if (calib_params)
    free(calib_params);
  if (pde_params)
    free(pde_params);
  if (autocal_params)
    free(autocal_params);

  if (sModel) {
    genmidat_free_model(sModel);
    free(sModel);
  }

  if (autocal_params) {
    genmidat_free_autocalinfos(sInfos);
    free(sInfos);
  }

  if (partialResults)
    free_dmatrix(partialResults, 0, 3, 0, 4);

  return err;
}

Err ko_pricer(
    CPD_STR cpd, CPD_UND und, otcpd_params *OTC_params,
    otcpd_precalc *OTC_precalc, SMILE_VOL_MARKET smile_mkt,
    SMILE_PARAMETERS smile_params, double pd_not, int mergeNtimes,
    double *mergeTimes, double *sigDom, double *sigFor, double *sigFx,
    double *corDF, double *corDFx, double *corFFx,
    long start_date, /*	Date at which initial notional exchange occurs */
    /*	Results */
    double **Results) {
  /* Declaration */
  int i, j, k, l, otc;

  /* For the Copula */
  double time1, time2;
  double stdDom, cumulDomForCor, stdFor, cumulForFxCor, meanFor, cumulDomFxCor;

  /* For the SSL */
  double error, forward;
  /* For the cvxt adjustment */

  double tstarTime;
  int tstarDate;

  double **RealisationMatrix = NULL;
  double ***correlations = NULL, **CORR = NULL, **cholMat = NULL;
  double var1, var2, var3, covar;

  double ***save_values = NULL;

  int *iOptimise = NULL;
  MCEBPARAMS sMCEBParams = NULL;

  Err err = NULL;

  if (OTC_params->KO_do_optim) {
    save_values =
        f3tensor(0, cpd->num_calls, 0, 1, 0, OTC_params->COPULAnSimul - 1);

    if (!save_values) {
      err = "Memory allocation failure in ko_pricer";
      goto FREE_RETURN;
    }
  }

  /* ==================================
  gaussian simulation
  ================================== */
  OTC_params->KO_FxMatrix =
      dmatrix(0, OTC_params->COPULAnSimul - 1, 0, cpd->num_calls - 1);
  OTC_params->KO_FxGauss =
      dmatrix(0, OTC_params->COPULAnSimul - 1, 0, cpd->num_calls - 1);
  OTC_params->KO_RatesMatrix =
      dmatrix(0, OTC_params->COPULAnSimul - 1, 0, 2 - 1);
  OTC_params->KO_RatesGauss =
      dmatrix(0, OTC_params->COPULAnSimul - 1, 0, 2 - 1);
  OTC_params->CUMUL = dmatrix(0, 2 - 1, 0, OTC_params->CUMULnPoints - 1);
  RealisationMatrix =
      dmatrix(0, OTC_params->COPULAnSimul - 1, 0, 2 + cpd->num_calls - 1);

  if (!RealisationMatrix || !OTC_params->KO_FxMatrix ||
      !OTC_params->KO_FxGauss || !OTC_params->KO_RatesMatrix ||
      !OTC_params->KO_RatesGauss) {
    err = "ko_pricer: memory allocation error (1)";
    goto FREE_RETURN;
  }

  time1 = clock();

  err = balsam_generation(OTC_params->COPULAnSimul, cpd->num_calls + 2,
                          RealisationMatrix);

  if (err)
    goto FREE_RETURN;

  time2 = clock();
  smessage("Phase 1 -BalSam generation        , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);
  time1 = clock();

  if (OTC_params->COPULAdo_pecs) {
    adjust_emp_covar(RealisationMatrix, OTC_params->COPULAnSimul,
                     cpd->num_calls + 2);
  }

  time2 = clock();
  smessage("Phase 2 -PECS adjustment        , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  /* Copy the Realisation Matrix */
  for (i = 0; i < OTC_params->COPULAnSimul; i++) {
    for (j = 0; j < cpd->num_calls; j++) {
      OTC_params->KO_FxGauss[i][j] = RealisationMatrix[i][j];
    }

    OTC_params->KO_RatesGauss[i][0] = RealisationMatrix[i][cpd->num_calls];
    OTC_params->KO_RatesGauss[i][1] = RealisationMatrix[i][cpd->num_calls + 1];
  }

  if (RealisationMatrix)
    free_dmatrix(RealisationMatrix, 0, OTC_params->COPULAnSimul - 1, 0,
                 2 + cpd->num_calls - 1);

  /*
  if(1 == 2)
  {
          OTCutilsSaveMatrix(	OTC_params->Options_number        ,
                                                  OTC_params->COPULAnSimul ,
                                                  OTC_params->KO_FxGauss , 2 ,
                                                  OTC_params->COPULAnSimul ,
                                                  OTC_params->KO_RatesGauss , 0
  ,0        ,NULL);
  }
  */

  /* ==================================
  GET ALL THE FX CORRELATIONS
  ================================== */
  correlations = f3tensor(0, 9, 0, 9, 0, mergeNtimes - 1);

  if (!correlations) {
    err = "Memory allocation failure in ko_pricer (2)";
    goto FREE_RETURN;
  }

  for (i = 0; i < mergeNtimes; i++) {
    correlations[0][1][i] = corDF[i];  // Dom1 For1
    correlations[0][2][i] = corDFx[i]; // Dom1 Fx1
    correlations[0][3][i] = 1.0;       // Dom1 Dom2
    correlations[0][4][i] = corDF[i];  // Dom1 For2
    correlations[0][5][i] = corDFx[i]; // Dom1 Fx2

    correlations[1][2][i] = corFFx[i]; // For1 Fx1
    correlations[1][3][i] = corDF[i];  // For1 Dom2
    correlations[1][4][i] = 1.0;       // For1 For2
    correlations[1][5][i] = corFFx[i]; // For1 Fx2

    correlations[2][3][i] = corDFx[i]; // Fx1 Dom2
    correlations[2][4][i] = corFFx[i]; // Fx1 For2
    correlations[2][5][i] = 1.0;       // Fx1 Fx2

    correlations[3][4][i] = corDF[i];  // Dom2 For2
    correlations[3][5][i] = corDFx[i]; // Dom2 Fx2

    correlations[4][5][i] = corFFx[i]; // For2 Fx2
  }

  cholMat = dmatrix(0, cpd->num_calls - 1, 0, cpd->num_calls - 1);
  CORR = dmatrix(0, cpd->num_calls - 1, 0, cpd->num_calls - 1);

  if (!cholMat || !CORR) {
    err = "ko_pricer: memory allocation error (3)";
    goto FREE_RETURN;
  }

  for (i = 0; i < cpd->num_calls; i++) {
    CORR[i][i] = 1.0;

    for (j = i + 1; j < cpd->num_calls; j++) {
      err =
          Fx3DtsFxFwdCov_corr(cpd->call[i].ex_time, cpd->call[j].ex_time, 0.0,
                              cpd->call[i].ex_time, mergeTimes, mergeNtimes,
                              /* FX1 or LGM1 */
                              sigDom, und->lda_dom, sigFor, und->lda_for, sigFx,
                              /* FX2 or LGM2 */
                              sigDom, und->lda_dom, sigFor, und->lda_for, sigFx,
                              correlations, &var1, &var2, &covar);

      if (err)
        goto FREE_RETURN;

      err = Fx3DtsImpliedVol_corr(cpd->call[j].ex_time, cpd->call[i].ex_time,
                                  cpd->call[j].ex_time, mergeTimes, mergeNtimes,
                                  sigDom, und->lda_dom, sigFor, und->lda_for,
                                  mergeTimes, sigFx, mergeNtimes, mergeTimes,
                                  corDF, corDFx, corFFx, mergeNtimes, &var3);

      if (err)
        goto FREE_RETURN;

      var2 += var3 * var3 * (cpd->call[j].ex_time - cpd->call[i].ex_time);

      if (fabs(var1 * var2) < 1.0e-15) {
        CORR[i][j] = 0.0;
      } else {
        CORR[i][j] = covar / sqrt(var1 * var2);
      }

      CORR[j][i] = CORR[i][j];
    }
  }

  /* ==================================
  CORRELATE ALL THE FFX GAUSSIAN MATRIX
  ================================== */
  nr_choldc(cpd->num_calls, CORR, cholMat);

  for (i = 0; i < OTC_params->COPULAnSimul; i++) {
    for (k = 0; k < cpd->num_calls; k++) {
      for (l = 0; l <= k; l++) {
        OTC_params->KO_FxMatrix[i][k] +=
            cholMat[k][l] * OTC_params->KO_FxGauss[i][l];
      }
    }
  }

  /* Copy the correlated brownian in the gaussian matrix */
  for (i = 0; i < OTC_params->COPULAnSimul; i++) {
    for (k = 0; k < cpd->num_calls; k++) {
      OTC_params->KO_FxGauss[i][k] = OTC_params->KO_FxMatrix[i][k];
    }
  }

  /* ==================================
  GENERATE ALL FX UNDER THE CORRESPONDING QTval FOR EACH FX
  ================================== */
  for (otc = 0; otc < cpd->num_calls; otc++) {
    /* GET SABR PARAMS AT KO DATE */
    /* CALIBRATE A SMILE MODEL */
    /* GET THE CUMUL */
    /* INVERT THE CUMUL */

    OTC_params->OTC = otc;

    /* ==================================
    get the measure time and date
    ================================== */
    tstarDate = cpd->call[otc].ex_date;
    tstarTime = cpd->call[otc].ex_time;

    /* ==================================
    get the Sabr parameters at the KO date
    ================================== */
    forward = und->spot_fx * swp_f_df(und->today, tstarDate, und->for_yc) /
              swp_f_df(und->today, tstarDate, und->dom_yc);

    err =
        OTCfwdSABR(forward, 0.0, cpd->call[otc].ex_time, tstarTime, smile_mkt,
                   smile_params, mergeTimes, mergeNtimes, sigDom, und->lda_dom,
                   sigFor, und->lda_for, sigFx, corDF, corDFx, corFFx, 1234);

    if (err)
      goto FREE_RETURN;

    /* ==================================
    Calibrate the Fx smile
    ================================== */
    if (cpd->call[otc].ex_time > 1.0e-10) {
      if (OTC_params->OTCsmile) {
        err = OTCcalibSmileModel(cpd->call[otc].ex_time, /* Maturity in years */
                                 forward,                /* Forward */
                                 OTC_params, smile_params, &error);

        if (err)
          goto FREE_RETURN;

        time1 = clock();

        /* ==================================
        Get the Cumulative Fx distribution under Qtstar
        ==================================*/
        err = OTCgetCumulative(forward, cpd->call[otc].ex_time, OTC_params);

        if (err)
          goto FREE_RETURN;

        /* ==================================
        Get the Fx under QTval
        ==================================*/
        err = OTCutilsGetFxRealisations(OTC_params);

        if (err)
          goto FREE_RETURN;
      } else {
        for (i = 0; i < OTC_params->COPULAnSimul; i++) {
          OTC_params->KO_FxMatrix[i][OTC_params->OTC] =
              forward *
              exp(-0.5 * smile_params->sigma * smile_params->sigma *
                      cpd->call[otc].ex_time +
                  OTC_params->KO_FxGauss[i][OTC_params->OTC] *
                      smile_params->sigma * sqrt(cpd->call[otc].ex_time));
        }
      }

      time2 = clock();
      smessage("Phase 3 -Simulation Finished        , time in sec: %.2f",
               (double)(time2 - time1) / CLOCKS_PER_SEC);
      time1 = clock();
    } else {
      for (i = 0; i < OTC_params->COPULAnSimul; i++)
        OTC_params->KO_FxMatrix[i][OTC_params->OTC] = forward;
    }
  }

  /*
  if(1 == 2)
  {
          OTCutilsSaveMatrix(	OTC_params->Options_number        ,
                                                  OTC_params->COPULAnSimul ,
                                                  OTC_params->KO_FxGauss ,
                                                  OTC_params->Options_number ,
                                                  OTC_params->COPULAnSimul ,
                                                  OTC_params->KO_FxMatrix ,
                                                  OTC_params->Options_number ,
                                                  OTC_params->Options_number ,
                                                  CORR);
  }
  */

  /* ==================================
  LOOP ON ALL COUPONS
  ================================== */
  for (otc = 0; otc < cpd->num_calls; otc++) {
    OTC_params->OTC = otc;

    /* ==================================
    get the measure time and date
    ================================== */
    tstarDate = cpd->call[otc].ex_date;
    tstarTime = cpd->call[otc].ex_time;

    /* ==================================
    get the 3F Cumulative correlations at the call date
    ================================== */
    err =
        OTCgetMoments(cpd->call[otc].ex_time, und->lda_dom, und->lda_for,
                      tstarTime, mergeTimes, mergeNtimes, sigDom, sigFor, sigFx,
                      corDF, corDFx, corFFx, &stdDom, &meanFor, &stdFor,
                      &cumulDomForCor, &cumulDomFxCor, &cumulForFxCor);

    if (err)
      goto FREE_RETURN;

    /* ==================================
    Correlate the rates and Fx
    ==================================*/
    err = OTCutilsCorrelateFxAndRates(OTC_params, 0.0, stdDom, meanFor, stdFor,
                                      cumulDomForCor, cumulDomFxCor,
                                      cumulForFxCor);

    if (err)
      goto FREE_RETURN;

    /*
    if(1 == 2)
    {
            OTCutilsSaveMatrix(	OTC_params->Options_number        ,
                                                    OTC_params->COPULAnSimul ,
                                                    OTC_params->KO_FxGauss , 2 ,
                                                    OTC_params->COPULAnSimul ,
                                                    OTC_params->KO_RatesGauss ,
                                                    2        ,
                                                    OTC_params->COPULAnSimul ,
                                                    OTC_params->KO_RatesMatrix);
    }
    */

    /* ==================================
    Payoff
    ================================== */

    if (OTC_params->KO_do_optim) {
      err = KOPDpayoffFee(cpd, und, OTC_params, OTC_precalc, smile_mkt,
                          smile_params, start_date, pd_not, tstarTime,
                          tstarDate, Results, save_values);

      if (err)
        goto FREE_RETURN;

      time2 = clock();
      smessage(
          "Phase 4 -Precalculation for optimisation        , time in sec: %.2f",
          (double)(time2 - time1) / CLOCKS_PER_SEC);
      time1 = clock();
    } else {
      err = KOPDpayoffFee(cpd, und, OTC_params, OTC_precalc, smile_mkt,
                          smile_params, start_date, pd_not, tstarTime,
                          tstarDate, Results, NULL);

      if (err)
        goto FREE_RETURN;
      time2 = clock();
      smessage("Phase 4 -Payoff evaluation Finished        , time in sec: %.2f",
               (double)(time2 - time1) / CLOCKS_PER_SEC);
    }
  }

  if (OTC_params->KO_do_optim) {
    time2 = clock();
    smessage(
        "Phase 4 -Precalculation for optimisation        , time in sec: %.2f",
        (double)(time2 - time1) / CLOCKS_PER_SEC);
    time1 = clock();
  } else {
    time2 = clock();
    smessage("Phase 4 -Payoff evaluation Finished        , time in sec: %.2f",
             (double)(time2 - time1) / CLOCKS_PER_SEC);
  }

  if (OTC_params->KO_do_optim) {
    OTC_params->KO_do_optim = 0;
    iOptimise = calloc(cpd->num_calls + 1, sizeof(int));
    sMCEBParams = calloc(1, sizeof(MCEBParams));

    if (!iOptimise || !sMCEBParams) {
      err = "Memory allocation faillure in Generic MidatAutocal";
      goto FREE_RETURN;
    }

    iOptimise[0] = 0;

    for (i = 0; i < cpd->num_calls; i++) {
      if (cpd->call[i].cxxall_type == 0)
        iOptimise[i + 1] = 1;
    }

    mceb_set_default_params(sMCEBParams);

    /* Set MCEB Params */
    sMCEBParams->iCallCurrent = 1;
    sMCEBParams->iIsKO = 1;

    sMCEBParams->iColPay = 0;

    sMCEBParams->iColBound = 0;
    sMCEBParams->iMultiIndex = 0;
    sMCEBParams->iNbIndex = 1;
    sMCEBParams->iRemoveLastOnLast = 1;
    sMCEBParams->iDoInfos = 0;
    sMCEBParams->iFindBestOptim = 0;

    err = mceb_allocate_params(sMCEBParams, cpd->num_calls + 1);
    if (err)
      goto FREE_RETURN;

    err = find_and_optimise_boundary(
        save_values, cpd->num_calls + 1, OTC_params->COPULAnSimul, iOptimise,
        sMCEBParams, &(Results[0][0]), &(Results[0][1]));

    if (err)
      goto FREE_RETURN;

    for (i = 0; i < cpd->num_calls; i++) {
      if (iOptimise[i + 1]) {
        cpd->call[i].barrier = log(sMCEBParams->dBarrier[i + 1] / und->spot_fx);
        cpd->call[i].orig_barrier = sMCEBParams->dBarrier[i + 1];
        Results[2 + i][7] = sMCEBParams->dBarrier[i + 1];
      } else
        Results[2 + i][7] = cpd->call[i].orig_barrier;
    }

    Results[0][0] = Results[0][1] = 0.0;

    time2 = clock();
    smessage("Phase 5 -Optimal Barrier found        , time in sec: %.2f",
             (double)(time2 - time1) / CLOCKS_PER_SEC);
    time1 = clock();

    for (otc = 0; otc < cpd->num_calls; otc++) {
      OTC_params->OTC = otc;

      /* ==================================
      get the measure time and date
      ================================== */
      tstarDate = cpd->call[otc].ex_date;
      tstarTime = cpd->call[otc].ex_time;

      /* ==================================
      get the 3F Cumulative correlations at the call date
      ================================== */
      err = OTCgetMoments(cpd->call[otc].ex_time, und->lda_dom, und->lda_for,
                          tstarTime, mergeTimes, mergeNtimes, sigDom, sigFor,
                          sigFx, corDF, corDFx, corFFx, &stdDom, &meanFor,
                          &stdFor, &cumulDomForCor, &cumulDomFxCor,
                          &cumulForFxCor);

      if (err)
        goto FREE_RETURN;

      /* ==================================
      Correlate the rates and Fx
      ==================================*/
      err = OTCutilsCorrelateFxAndRates(OTC_params, 0.0, stdDom, meanFor,
                                        stdFor, cumulDomForCor, cumulDomFxCor,
                                        cumulForFxCor);

      if (err)
        goto FREE_RETURN;

      /* ==================================
      Payoff
      ================================== */
      err = KOPDpayoffFee(cpd, und, OTC_params, OTC_precalc, smile_mkt,
                          smile_params, start_date, pd_not, tstarTime,
                          tstarDate, Results, NULL);

      if (err)
        goto FREE_RETURN;
    }
    time2 = clock();
    smessage("Phase 6 -Payoff evaluation Finished        , time in sec: %.2f",
             (double)(time2 - time1) / CLOCKS_PER_SEC);
  }

FREE_RETURN:

  if (correlations)
    free_f3tensor(correlations, 0, 9, 0, 9, 0, mergeNtimes - 1);
  if (cholMat)
    free_dmatrix(cholMat, 0, cpd->num_calls - 1, 0, cpd->num_calls - 1);
  if (CORR)
    free_dmatrix(CORR, 0, cpd->num_calls - 1, 0, cpd->num_calls - 1);

  if (iOptimise)
    free(iOptimise);

  if (sMCEBParams) {
    mceb_free_params(sMCEBParams);
    free(sMCEBParams);
  }

  if (save_values)
    free_f3tensor(save_values, 0, cpd->num_calls, 0, 1, 0,
                  OTC_params->COPULAnSimul - 1);

  return err;
}