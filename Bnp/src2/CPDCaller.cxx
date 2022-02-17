
#include "Fx3FBetaDLMCalculations.h"
#include "Fx3FBetaDLMCalibration.h"
#include "Fx3FBetaDLMUtil.h"

#include "CPDVol.h"
#include "OTCcaller.h"
#include "OTCutils.h"
#include "math.h"
#include "num_h_interp.h"
#include "opfnctns.h"
#include "opsabrgenericinterp.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srtaccess.h"

#define CPD_HUGE_FEES 9.0E+20
#define CPD_HUGE_BARRIER 999999

/*
Convert the funding leg into a domestic one
                -	spreads and margins
                -	past fixings
                -	final notional exchange
                -	initial notional exchange
*/
Err convert_funding_to_domestic(
    /*	Inputs */
    long today,                 /*	Today */
    long not_ex_date,           /*	Date at which the
                                                                  initial notional
                                                                  exchange takes
                             place           (or has taken place) */
    int eod_fix_flag,           /*	0: I        , 1: E */
    int eod_pay_flag,           /*	0: I        , 1: E */
    double fx_fund_dom,         /*	Fx fund/dom        , 2bd fwd */
    long fx_fund_dom_spot_date, /*	Spot date for Fx */
    double dom_not,             /*	Domestic notional */
    char *dom_yc,               /*	Domestic discount curve */
    int fund_ncpn,              /*	Number of coupon */
    long *fund_fix,             /*	Fixing dates */
    long *fund_start,           /*	Start dates */
    long *fund_pay,             /*	Pay dates */
    char **fund_basis,          /*	Basis */
    /*	The following are modified */
    char *fund_yc,        /*	Funding discount curve        ,
                                            changed to domestic */
    double *fund_not,     /*	Funding notional
                                                    in funding ccy        ,
                                            converted to domestic */
    double *fund_spr,     /*	Spread in funding ccy        ,
                                            put to 0 */
    double *fund_mrg,     /*	Margin in funding ccy        ,
                                            converted to margin over
                                            cash libor in domestic currency */
    double *fund_fix_cpn, /*	Fixing: contains spread
                                                    but not margin        ,
                                            converted to equivalent
                                            domestic cash-flows */
    /*	The following are returned */
    long *fund_start_date, /*	Start date of the funding */
    double *eq_final_ex,   /*	Domestic cash-flow equivalent
                                                             to final
                        exchange   (to be delivered   at funding start date) */
    double *
        eq_init_ex) /*	Domestic cash-flow equivalent
                                                            to initial exchange
                                                    (to be delivered
                                                    at initial exchange date) */
{
  int i, i0;
  double ffx;
  double not_ratio;

  /*	Find the start date */
  i = 0;
  while (i < fund_ncpn && fund_fix[i] < today + eod_fix_flag) {
    i++;
  }
  i0 = i;

  if (i0 == fund_ncpn) {
    *fund_start_date = fund_pay[i0 - 1];
  } else {
    *fund_start_date = fund_start[i0];
  }

  /*	Precalc notional ratio */
  not_ratio = *fund_not / dom_not;

  /*	Adjust fixings */
  for (i = 0; i < i0; i++) {
    if (fund_pay[i] >= today + eod_pay_flag) {

      ffx = fx_fund_dom *
            swp_f_df(fx_fund_dom_spot_date, fund_pay[i], fund_yc) /
            swp_f_df(fx_fund_dom_spot_date, fund_pay[i], dom_yc);
      fund_fix_cpn[i] = not_ratio * ffx * (fund_fix_cpn[i] + fund_mrg[i]);
    } else {
      fund_fix_cpn[i] = 0.0;
    }
    fund_mrg[i] = 0.0;
  }

  /*	Adjust spreads and margins */
  for (i = i0; i < fund_ncpn; i++) {
    ffx = fx_fund_dom * swp_f_df(fx_fund_dom_spot_date, fund_pay[i], fund_yc) /
          swp_f_df(fx_fund_dom_spot_date, fund_pay[i], dom_yc);

    fund_mrg[i] = not_ratio * ffx * (fund_spr[i] + fund_mrg[i]);
    fund_spr[i] = 0.0;
  }

  /*	Compute the domestic cash flow at start date equivalent to final
   * exchange */
  if (*fund_start_date >= today + eod_pay_flag) {
    ffx = fx_fund_dom *
          swp_f_df(fx_fund_dom_spot_date, *fund_start_date, fund_yc) /
          swp_f_df(fx_fund_dom_spot_date, *fund_start_date, dom_yc);
    *eq_final_ex = *fund_not * ffx;
  } else {
    *eq_final_ex = 0.0;
  }

  /*	Compute the domestic cash flow at start date equivalent to initial
   * exchange */
  if (not_ex_date >= today + eod_pay_flag) {
    ffx = fx_fund_dom * swp_f_df(fx_fund_dom_spot_date, not_ex_date, fund_yc) /
          swp_f_df(fx_fund_dom_spot_date, not_ex_date, dom_yc);
    *eq_init_ex = *fund_not * ffx;
  } else {
    *eq_init_ex = 0.0;
  }

  /*	Change funding yc and notional to domestic */
  strcpy(fund_yc, dom_yc);
  *fund_not = dom_not;

  return NULL;
}

/*	Caller for callable power duals */
/*	------------------------------- */

#define POS_VAL(X) ((X) > 0 ? (X) : 0)

#define CALL_VAL(FWD, STRIKE, D, S) ((FWD)*norm((D) + (S)) - (STRIKE)*norm((D)))

#define PUT_VAL(FWD, STRIKE, D, S)                                             \
  (-(FWD)*norm(-(D) - (S)) + (STRIKE)*norm(-(D)))

#define OPT_VAL_MACRO(TYPE, FWD, STRIKE, STD, HALF_STD)                        \
  ((TYPE) == 0                                                                 \
       ? 0.0                                                                   \
       : ((TYPE) == 1                                                          \
              ? POS_VAL((FWD) - (STRIKE))                                      \
              : ((TYPE) == 2                                                   \
                     ? POS_VAL((STRIKE) - (FWD))                               \
                     : ((TYPE) == 3 ? CALL_VAL((FWD), (STRIKE),                \
                                               log((FWD) / (STRIKE)) / (STD) - \
                                                   (HALF_STD),                 \
                                               (STD))                          \
                                    : PUT_VAL((FWD), (STRIKE),                 \
                                              log((FWD) / (STRIKE)) / (STD) -  \
                                                  (HALF_STD),                  \
                                              (STD))))))

Err transform_interp_coupon(int npts, double *pts, double *cpn, int linxtr_l,
                            int linxtr_r, double *cst, double *wspot,
                            int *first_pt, int *last_pt, double *weights) {
  int i;
  double tangent = 0.0;

  if ((linxtr_l || linxtr_r) && npts < 2)
    return serror("At least 2 points needed for extrapolation");

  *first_pt = 0;
  *wspot = 0.0;
  if (linxtr_l) {
    *first_pt = 1;
    *wspot = tangent = (cpn[1] - cpn[0]) / (pts[1] - pts[0]);
  }

  *cst = cpn[*first_pt] - *wspot * pts[*first_pt];
  for (i = *first_pt; i < npts - 1; i++) {
    weights[i] = -tangent;
    weights[i] += (tangent = (cpn[i + 1] - cpn[i]) / (pts[i + 1] - pts[i]));
  }
  if (linxtr_r)
    *last_pt = npts - 2;
  else {
    *last_pt = npts - 1;
    weights[npts - 1] = -tangent;
  }
  return NULL;
}

void transform_oldspec_into_newspec(double alpha, double beta, int floored,
                                    double floor, int capped, double cap,
                                    int *nstrikes, double *wcst, double *wspot,
                                    double *strikes, double *weights)

{
  *nstrikes = floored + capped;
  if (*nstrikes == 0) {
    *wcst = alpha;
    *wspot = beta;
  } else if (fabs(beta) < 1e-16) {
    *nstrikes = 0;
    *wcst = alpha;
    *wspot = 0.0;
    if (floored && *wcst < floor)
      *wcst = floor;
    if (capped && *wcst > cap)
      *wcst = cap;
  } else if (*nstrikes == 1) {
    strikes[0] = ((floored ? floor : cap) - alpha) / beta;
    if ((floored && beta < 0.0) || (capped && beta > 0.0)) {
      *wcst = alpha;
      *wspot = beta;
      weights[0] = -beta;
    } else // (floored && beta > 0.0) || (capped && beta < 0.0)
    {
      *wcst = (floored ? floor : cap);
      *wspot = 0.0;
      weights[0] = beta;
    }
  } else // *nstrikes == 2
  {
    *wspot = 0.0;
    *wcst = (beta < 0.0 ? cap : floor);
    strikes[0] = ((beta < 0.0 ? cap : floor) - alpha) / beta;
    strikes[1] = ((beta < 0.0 ? floor : cap) - alpha) / beta;
    weights[0] = beta;
    weights[1] = -beta;
  }
}

Err cpd_caller(
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
    double *correl_dom_fx, double *correl_for_fx, long corr_n_times,
    CPDBETADLMPARAMS cpd_dlm_params,
    Err (*get_ir_cash_vol)(/*	Function to get IR cash vol from the markets */
                           char *vol_curve_name, double start_date,
                           double end_date, double cash_strike, int zero,
                           char *ref_rate_name, double *vol, double *power),
    /*	Fx vol from the market */
    long *fx_mkt_vol_date, double *fx_mkt_vol, int num_fx_mkt_vol,
    /*	Fx SABR parameters from the market */
    double *fx_mkt_smile_alpha, double *fx_mkt_smile_beta,
    double *fx_mkt_smile_rho, double *fx_mkt_smile_pi,
    int use_sabr,      /*	0: no smile adj-t        , 1: only for underlying      , 2:
                    fees      adj-t for the call/KO */
    int use_3F_interp, /*	0: interpolates the ATM market vols linearly ,
                    1: uses the 3F for the interpolation */
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
    double *pd_not_ref_fix_fx,   //	fx fixings FOR EACH CALL DATE + in the
                                 // end if relevant
    /*		calls */
    int *call_type,  /*	0: call        , 1: KO */
    int ncall,       /*	Number of calls */
    int pay_rec,     /*	0: rec pd        , 1: pay pd */
    long *ex_date,   /*	Call dates */
    long *set_date,  /*	Settlement dates */
    double *barrier, /*	in case of a pure KO or a Callable KO */
    int *bar_type,   /*	0: up and in        , 1: down and in */
    double *fees,    /*  fees if deal is called in domestic currency */
    /* TARN */
    int TARN_Do, double TARN_Floor,
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
    /*	Vega */
    double dom_vol_shift, double for_vol_shift, double fx_vol_shift,
    /*	Exercised flag */
    int exercised,    /*	Flag */
    long ex_date_ex,  /*	Date when exercised */
    long ex_date_set, /*	Corresponding settlement date */
    /*	Prune calls */
    int prune_calls,    /*	Flag: 0 no prune        , 1 do prune */
    int no_prune_years, /*	Number of years from the next call to the first
                     prune */
    int prune_factor,   /*	Ex: 2 -> one call out of two */

    /* Disable call dates */
    int disable_calls_method, /* Name of the disable call method used : 0 None
                                 , 1 All in range 2 list of call ( Parameter1
                           and Parameter 2 not used */
    int first_call_index,     /* Parameter 1 of the disable method */
    int last_call_index,      /* Parameter 2 of the disable method */

    int iNbIndex_List, /* Number of erased call in the list */
    int *Index_List,   /* Array of index of calls to be erased */

    int erasing_call_done, // always to 0        , switch to 1 after haveng
                           // disabled calls
    int *erased_call_list, // list of 0 & 1 used only if erasing_call_done = 1
                           // and if use_GMA = TRUE

    /* Alpha Beta Smile Model */
    int use_smile, /* Flag: 0 no use alphabeta model        , 1 use of the
                alphabeta model */
    double alpha,  /* Alpha        , only used if use_smile =1 */
    double beta,   /* Beta        ,  only used if use_smile =1 */
    /*For the smile impact*/
    int use_GMA,
    int use_3f_optim_barrier, // Optimize the boundary with the 3F and then
                              // price the smile adjustment
    double SSLstd, int SSLniter, int CopNsimul, int CopDoPecs,
    int CummulNpoints, int CummulNstd, double CummulPrecision, int CummulLinear,
    int PayoffFunction, double fwdSmileVisuNstd, int smileOtc, int smileFee,
    int fwdVolMethod, int smileModel, double BMpi, int otc, int FundingSpeedUp,
    /* Do not calculate all OTC*/
    int nStart, int oneOutOfN,
    /*Fast MC*/
    int FMC_do, double FMC_precision, int FMC_min_paths,
    /*	Results */
    double *fund_val,         /*	Value of the funding leg */
    double *pd_val,           /*	Value of the Power Dual leg */
    double *call_val,         /*	Value of the callable feature */
    double *call_stdev,       /*	Standard deviation of the call if applicable */
    double *smile_adjustment, /*	GMA smile adjustment */
    double ***optim_bar,      /*	Contains the value of the optimal KO for
                           corresponding calls */
    double **GMA_Results,     /*	Contains the OTC/KO info */
    int export_ts,            /*	1: Export TS        , 0: don't */
    CPD_UND und_exp)          /*	TS to be exported */
{
  cpd_str cpd_, *cpd = &cpd_;
  cpd_und und_, *und = &und_;
  cpd_tree_arg tree_arg_, *tree_arg = &tree_arg_;
  cpd_mc_arg mc_arg_, *mc_arg = &mc_arg_;

  int call_feat;
  double fund_leg, partial_fund_leg;
  char *fund_yc;
  double fund_mult;
  SrtBasisCode bas;

  double pd_leg, partial_pd_leg;
  int type;
  double df, fwd, adj, std, half_std, str, floor, cap, fixing;
  PD_EXO_CPN cpn;

  double call, call_std, opt_string, smile_opt_string;
  int i, j, str_idx, ex_idx;
  int pd_fut_ncpn, fund_fut_ncpn;
  int free_struct = 0;
  int skip_fill = 0;

  int partial_pv, partial_fund_idx, partial_pd_idx;
  double temp;

  int for_fund;
  long fund_start_date, fin_not_date;
  double eq_final_ex, eq_init_ex;

  int is_ko;
  long next_call_date;
  long prune_date;
  int new_ncall;
  int *new_call_type = NULL;
  long *new_ex_date = NULL;
  long *new_set_date = NULL;
  double *new_barrier = NULL;
  int *new_bar_type = NULL;
  double *new_fees = NULL;

  double *adj_fees = NULL;
  double adj_fee_not = 0.0;
  double cum_adj_fee;
  double smile_std, half_smile_std, interp_coef;
  double smile_floor, smile_cap;
  CPD_PAY_ARG cpd_prm;

  double dlm_fwd, dlm_floor, dlm_cap, dlm_temp;
  double dlm_pd_leg, dlm_partial_pd_leg, saved_pd_leg, saved_partial_pd_leg;
  FxBetaDLM_FxOptInst *Inst = NULL;
  FxBetaDLM_InstPrecalc *InstConst = NULL;
  SrtCallPutType InstCallType;
  double InstStrike;

  Err err = NULL;

  double thirdCcyfundNot;
  char Copyfund_ccy_yc[255];

  /* for alpha fudge */
  int UpOrDown;
  double call_val_alpha;     /*	Value of the callable feature */
  double call_val_alpha_std; /*	Value of the callable features std*/

  double forward_at_smile, vol_at_smile, alpha_at_smile, beta_at_smile,
      rho_at_smile;

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

  // for pruning when interpolating coupon:
  int *new_pd_not_ref_nfxpts = NULL;
  double **new_pd_not_ref_fxpts = NULL;
  double **new_pd_not_ref_cpn_at_pts = NULL;
  int *new_pd_not_ref_lin_xtrpl_l = NULL;
  int *new_pd_not_ref_lin_xtrpl_r = NULL;
  double *new_pd_not_ref_fix_fx = NULL;
  long *new_pd_not_ref_fix = NULL;

  /* For Erasing calls */
  int iErasedIndex;

  // VOL MKT STRUCTURE
  SMILE_VOL_MARKET smile_mkt = NULL;
  SMILE_PARAMETERS smile_params = NULL;
  double *fwd_at_vol_dates = NULL;

  /* ===================================
  ======================================
          FOR THE SMILE ADJUSTMENT
  ======================================
  =================================== */
  /*	Results */
  double **Results = NULL;
  double *mergeTimes = NULL, *sigDom = NULL, *sigFor = NULL, *sigFx = NULL,
         *corDF = NULL, *corDFx = NULL, *corFFx = NULL;
  int mergeNtimes;
  int Height = ncall + 4;
  int Width = 12;
  int Smile_model, Smile_OTC, Smile_Fee_OTC, Smile_fwd_vol;

  otcpd_params *OTC_params = NULL;
  otcpd_precalc *OTC_precalc = NULL;

  if ((use_cpn_opt_str || use_not_opt_str) &&
      (cpd_dlm_params->use_beta_dlm || use_smile)) {
    err = serror("Interpolated coupon is temporarily not supported with "
                 "BetaDLM and AlphaBeta models");
    goto FREE_RETURN;
  }
  /*
          if (ncall > 0 && ex_date[0] < start_date)
          {
                  err = serror("First exercise date is before start date: not
     allowed"); goto FREE_RETURN;
          }
  */
  if (use_GMA && use_sabr == 2) {
    err = serror("Smile impact with 3F fees and GMA should not be calculated "
                 "simultaneously");
    goto FREE_RETURN;
  }

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

  // Calculate max_interp_pts:

  if (use_cpn_opt_str)
    for (i = 0; i < pd_ncpn; i++)
      if (pd_nfxpts[i] > max_interp_pts)
        max_interp_pts = pd_nfxpts[i];
  if (use_not_opt_str)
    for (i = 0; i <= ncall; i++)
      if (pd_not_ref_nfxpts[i] > max_interp_pts)
        max_interp_pts = pd_not_ref_nfxpts[i];

  // Allocate memory for new_*** used for pruning:

  new_call_type = (int *)calloc(alloc_ncall, sizeof(int));
  new_ex_date = (long *)calloc(alloc_ncall, sizeof(long));
  new_set_date = (long *)calloc(alloc_ncall, sizeof(long));
  new_barrier = (double *)calloc(alloc_ncall, sizeof(double));
  new_bar_type = (int *)calloc(alloc_ncall, sizeof(int));
  new_fees = (double *)calloc(alloc_ncall, sizeof(double));

  new_pd_not_ref_fix = (long *)calloc(alloc_ncall + 1, sizeof(double));
  new_pd_not_ref_fix_fx = (double *)calloc(alloc_ncall + 1, sizeof(double));

  if (use_not_opt_str) {
    new_pd_not_ref_nfxpts = (int *)calloc(alloc_ncall + 1, sizeof(int));
    new_pd_not_ref_fxpts = dmatrix(0, alloc_ncall, 0, max_interp_pts - 1);
    new_pd_not_ref_cpn_at_pts = dmatrix(0, alloc_ncall, 0, max_interp_pts - 1);
    new_pd_not_ref_lin_xtrpl_l = (int *)calloc(alloc_ncall + 1, sizeof(int));
    new_pd_not_ref_lin_xtrpl_r = (int *)calloc(alloc_ncall + 1, sizeof(int));
  }
  if (!new_call_type || !new_ex_date || !new_set_date || !new_barrier ||
      !new_bar_type || !new_fees || !new_pd_not_ref_fix ||
      !new_pd_not_ref_fix_fx ||
      (use_not_opt_str &&
       (!new_pd_not_ref_nfxpts || !new_pd_not_ref_fxpts ||
        !new_pd_not_ref_cpn_at_pts || !new_pd_not_ref_lin_xtrpl_l ||
        !new_pd_not_ref_lin_xtrpl_r))) {
    err = serror("Memory allocation failure in cpd_caller");
    goto FREE_RETURN;
  }

  /*	If pruned and at least 1 call after today */
  /*	Automatic pruning */
  if (prune_calls && ncall > 0 && ex_date[ncall - 1] > today + eod_ex_flag &&
      prune_factor > 0) {
    //	Check prune factor
    if (prune_factor < 1)
      prune_factor = 1;

    /*	Calc next call date */
    i = 0;
    while (ex_date[i] < today + eod_ex_flag)
      i++;
    next_call_date = ex_date[i];

    /*	Check for pure KO feature */
    is_ko = 0;
    for (j = i; j < ncall; j++)
      if (call_type[j] != 0)
        is_ko = 1;

    /*	Do not prune if structure is a full ko */
    if (!is_ko) {
      /*	Calc prune date */
      prune_date = add_unit(next_call_date, no_prune_years, SRT_YEAR,
                            MODIFIED_SUCCEEDING);

      /*	Make sure there is at least 1 call after prune date */
      if (ex_date[ncall - 1] > prune_date) {
        i = 0;
        while (ex_date[i] <= prune_date) {
          new_call_type[i] = call_type[i];
          new_ex_date[i] = ex_date[i];
          new_set_date[i] = set_date[i];
          new_barrier[i] = barrier[i];
          new_bar_type[i] = bar_type[i];
          new_fees[i] = fees[i];
          new_pd_not_ref_fix[i] = pd_not_ref_fix[i];
          new_pd_not_ref_fix_fx[i] = pd_not_ref_fix_fx[i];
          if (use_not_opt_str) {
            new_pd_not_ref_nfxpts[i] = pd_not_ref_nfxpts[i];
            memcpy(new_pd_not_ref_fxpts[i], pd_not_ref_fxpts[i],
                   pd_not_ref_nfxpts[i] * sizeof(double));
            memcpy(new_pd_not_ref_cpn_at_pts[i], pd_not_ref_cpn_at_pts[i],
                   pd_not_ref_nfxpts[i] * sizeof(double));
            new_pd_not_ref_lin_xtrpl_l[i] = pd_not_ref_lin_xtrpl_l[i];
            new_pd_not_ref_lin_xtrpl_r[i] = pd_not_ref_lin_xtrpl_r[i];
          }
          i++;
        }
        j = i - 1 + prune_factor;
        while (j < ncall) {
          new_call_type[i] = call_type[j];
          new_ex_date[i] = ex_date[j];
          new_set_date[i] = set_date[j];
          new_barrier[i] = barrier[j];
          new_bar_type[i] = bar_type[j];
          new_fees[i] = fees[j];
          new_pd_not_ref_fix[i] = pd_not_ref_fix[j];
          new_pd_not_ref_fix_fx[i] = pd_not_ref_fix_fx[j];
          if (use_not_opt_str) {
            new_pd_not_ref_nfxpts[i] = pd_not_ref_nfxpts[j];
            memcpy(new_pd_not_ref_fxpts[i], pd_not_ref_fxpts[j],
                   pd_not_ref_nfxpts[j] * sizeof(double));
            memcpy(new_pd_not_ref_cpn_at_pts[i], pd_not_ref_cpn_at_pts[j],
                   pd_not_ref_nfxpts[j] * sizeof(double));
            new_pd_not_ref_lin_xtrpl_l[i] = pd_not_ref_lin_xtrpl_l[j];
            new_pd_not_ref_lin_xtrpl_r[i] = pd_not_ref_lin_xtrpl_r[j];
          }
          i++;
          j += prune_factor;
        }
        new_ncall = i;

        new_pd_not_ref_fix[new_ncall] = pd_not_ref_fix[ncall];
        new_pd_not_ref_fix_fx[new_ncall] = pd_not_ref_fix_fx[ncall];
        if (use_not_opt_str &
            (new_ncall != ncall)) // copy final redemption information:
        {
          new_pd_not_ref_nfxpts[new_ncall] = pd_not_ref_nfxpts[ncall];
          memcpy(new_pd_not_ref_fxpts[new_ncall], pd_not_ref_fxpts[ncall],
                 pd_not_ref_nfxpts[ncall] * sizeof(double));
          memcpy(new_pd_not_ref_cpn_at_pts[new_ncall],
                 pd_not_ref_cpn_at_pts[ncall],
                 pd_not_ref_nfxpts[ncall] * sizeof(double));
          new_pd_not_ref_lin_xtrpl_l[new_ncall] = pd_not_ref_lin_xtrpl_l[ncall];
          new_pd_not_ref_lin_xtrpl_r[new_ncall] = pd_not_ref_lin_xtrpl_r[ncall];
        }

        prune_calls = 0;
        no_prune_years = 999;
        prune_factor = 1;
        err = cpd_caller(
            today, use_calib, fx_spot, fx_spot_date, dom_calib, dom_und, dom_yc,
            dom_vc, dom_ref, dom_swap_freq, dom_swap_basis, dom_lam, for_calib,
            for_und, for_yc, for_vc, for_ref, for_swap_freq, for_swap_basis,
            for_lam, min_fact, max_fact, use_jumps, corr_times, correl_dom_for,
            correl_dom_fx, correl_for_fx, corr_n_times, cpd_dlm_params,
            get_ir_cash_vol, fx_mkt_vol_date, fx_mkt_vol, num_fx_mkt_vol,
            fx_mkt_smile_alpha, fx_mkt_smile_beta, fx_mkt_smile_rho,
            fx_mkt_smile_pi, use_sabr, use_3F_interp, smile_spec_type, fx3dund,
            start_date, fund_not, fund_ccy, fund_ccy_yc, fx_fund_dom,
            fx_fund_dom_spot_date, fund_ncpn, fund_fix, fund_start, fund_pay,
            fund_basis, fund_spr, fund_mrg, fund_fix_cpn, pd_not, pd_ncpn,
            pd_fix, pd_start, pd_pay, pd_basis, pd_alpha, pd_beta, pd_floored,
            pd_floor, pd_capped, pd_cap, pd_nfxpts, pd_fxpts, pd_cpn_at_pts,
            pd_lin_xtrpl_l, pd_lin_xtrpl_r, pd_fix_fx, new_pd_not_ref_fix,
            pd_not_ref_alpha, pd_not_ref_beta, pd_not_ref_floored,
            pd_not_ref_floor, pd_not_ref_capped, pd_not_ref_cap,
            new_pd_not_ref_nfxpts, new_pd_not_ref_fxpts,
            new_pd_not_ref_cpn_at_pts, new_pd_not_ref_lin_xtrpl_l,
            new_pd_not_ref_lin_xtrpl_r, new_pd_not_ref_fix_fx, new_call_type,
            new_ncall, pay_rec, new_ex_date, new_set_date, new_barrier,
            new_bar_type, new_fees, TARN_Do, 0.0, req_stp, req_pth, bar_smooth,
            do_pecs, forcetree, do_optim, force_optim, fx_bound, use_bound,
            do_infos, eod_fix_flag, eod_pay_flag, eod_ex_flag, dom_vol_shift,
            for_vol_shift, fx_vol_shift, exercised, ex_date_ex, ex_date_set,
            prune_calls, no_prune_years, prune_factor,

            disable_calls_method, first_call_index, last_call_index,

            iNbIndex_List, Index_List,

            erasing_call_done, erased_call_list,

            use_smile, alpha, beta,

            /*For the smile impact*/
            use_GMA, use_3f_optim_barrier, SSLstd, SSLniter, CopNsimul,
            CopDoPecs, CummulNpoints, CummulNstd, CummulPrecision, CummulLinear,
            PayoffFunction, fwdSmileVisuNstd, smileOtc, smileFee, fwdVolMethod,
            smileModel, BMpi, otc, FundingSpeedUp,
            /* Do not calculate all OTC*/
            nStart, oneOutOfN,
            // Fast MC
            FMC_do, FMC_precision, FMC_min_paths, fund_val, pd_val, call_val,
            call_stdev, smile_adjustment, optim_bar, GMA_Results, export_ts,
            und_exp);

        goto FREE_RETURN;
      } else {
        smessage("Prune date is on or after last call - no prune");
      }
    } else {
      smessage("Structure has KO feature - no prune");
    }
  }

  /*	Manual pruning */
  if (prune_calls && ncall > 0 && ex_date[ncall - 1] > today + eod_ex_flag &&
      prune_factor == -1) {
    j = 0;
    for (i = 0; i < ncall; i++) {
      if (call_type[i] != -1) {
        new_call_type[j] = call_type[i];
        new_ex_date[j] = ex_date[i];
        new_set_date[j] = set_date[i];
        new_barrier[j] = barrier[i];
        new_bar_type[j] = bar_type[i];
        new_fees[j] = fees[i];
        new_pd_not_ref_fix[j] = pd_not_ref_fix[i];
        new_pd_not_ref_fix_fx[j] = pd_not_ref_fix_fx[i];
        if (use_not_opt_str) {
          new_pd_not_ref_nfxpts[j] = pd_not_ref_nfxpts[i];
          memcpy(new_pd_not_ref_fxpts[j], pd_not_ref_fxpts[i],
                 pd_not_ref_nfxpts[i] * sizeof(double));
          memcpy(new_pd_not_ref_cpn_at_pts[j], pd_not_ref_cpn_at_pts[i],
                 pd_not_ref_nfxpts[i] * sizeof(double));
          new_pd_not_ref_lin_xtrpl_l[j] = pd_not_ref_lin_xtrpl_l[i];
          new_pd_not_ref_lin_xtrpl_r[j] = pd_not_ref_lin_xtrpl_r[i];
        }
        j++;
      }
    }

    new_ncall = j;

    new_pd_not_ref_fix[new_ncall] = pd_not_ref_fix[ncall];
    new_pd_not_ref_fix_fx[new_ncall] = pd_not_ref_fix_fx[ncall];
    if (use_not_opt_str &
        (new_ncall != ncall)) // copy final redemption information:
    {
      new_pd_not_ref_nfxpts[new_ncall] = pd_not_ref_nfxpts[ncall];
      memcpy(new_pd_not_ref_fxpts[new_ncall], pd_not_ref_fxpts[ncall],
             pd_not_ref_nfxpts[ncall] * sizeof(double));
      memcpy(new_pd_not_ref_cpn_at_pts[new_ncall], pd_not_ref_cpn_at_pts[ncall],
             pd_not_ref_nfxpts[ncall] * sizeof(double));
      new_pd_not_ref_lin_xtrpl_l[new_ncall] = pd_not_ref_lin_xtrpl_l[ncall];
      new_pd_not_ref_lin_xtrpl_r[new_ncall] = pd_not_ref_lin_xtrpl_r[ncall];
    }

    prune_calls = 0;
    no_prune_years = 999;
    prune_factor = 1;
    err = cpd_caller(
        today, use_calib, fx_spot, fx_spot_date, dom_calib, dom_und, dom_yc,
        dom_vc, dom_ref, dom_swap_freq, dom_swap_basis, dom_lam, for_calib,
        for_und, for_yc, for_vc, for_ref, for_swap_freq, for_swap_basis,
        for_lam, min_fact, max_fact, use_jumps, corr_times, correl_dom_for,
        correl_dom_fx, correl_for_fx, corr_n_times, cpd_dlm_params,
        get_ir_cash_vol, fx_mkt_vol_date, fx_mkt_vol, num_fx_mkt_vol,
        fx_mkt_smile_alpha, fx_mkt_smile_beta, fx_mkt_smile_rho,
        fx_mkt_smile_pi, use_sabr, use_3F_interp, smile_spec_type, fx3dund,
        start_date, fund_not, fund_ccy, fund_ccy_yc, fx_fund_dom,
        fx_fund_dom_spot_date, fund_ncpn, fund_fix, fund_start, fund_pay,
        fund_basis, fund_spr, fund_mrg, fund_fix_cpn, pd_not, pd_ncpn, pd_fix,
        pd_start, pd_pay, pd_basis, pd_alpha, pd_beta, pd_floored, pd_floor,
        pd_capped, pd_cap, pd_nfxpts, pd_fxpts, pd_cpn_at_pts, pd_lin_xtrpl_l,
        pd_lin_xtrpl_r, pd_fix_fx, new_pd_not_ref_fix, pd_not_ref_alpha,
        pd_not_ref_beta, pd_not_ref_floored, pd_not_ref_floor,
        pd_not_ref_capped, pd_not_ref_cap, new_pd_not_ref_nfxpts,
        new_pd_not_ref_fxpts, new_pd_not_ref_cpn_at_pts,
        new_pd_not_ref_lin_xtrpl_l, new_pd_not_ref_lin_xtrpl_r,
        new_pd_not_ref_fix_fx, new_call_type, new_ncall, pay_rec, new_ex_date,
        new_set_date, new_barrier, new_bar_type, new_fees,
        /* TARN */
        TARN_Do, 0.0, req_stp, req_pth, bar_smooth, do_pecs, forcetree,
        do_optim, force_optim, fx_bound, use_bound, do_infos, eod_fix_flag,
        eod_pay_flag, eod_ex_flag, dom_vol_shift, for_vol_shift, fx_vol_shift,
        exercised, ex_date_ex, ex_date_set, prune_calls, no_prune_years,
        prune_factor,

        disable_calls_method, first_call_index, last_call_index,

        iNbIndex_List, Index_List,

        erasing_call_done, erased_call_list,

        use_smile, alpha, beta,

        /*For the smile impact*/
        use_GMA, use_3f_optim_barrier, SSLstd, SSLniter, CopNsimul, CopDoPecs,
        CummulNpoints, CummulNstd, CummulPrecision, CummulLinear,
        PayoffFunction, fwdSmileVisuNstd, smileOtc, smileFee, fwdVolMethod,
        smileModel, BMpi, otc, FundingSpeedUp,
        /* Do not calculate all OTC*/
        nStart, oneOutOfN,
        // Fast MC
        FMC_do, FMC_precision, FMC_min_paths,

        fund_val, pd_val, call_val, call_stdev, smile_adjustment, optim_bar,
        GMA_Results, export_ts, und_exp);

    goto FREE_RETURN;
  }

  /*	Disable Calls Method */
  if (disable_calls_method && ncall > 0 &&
      ex_date[ncall - 1] > today + eod_ex_flag) {
    j = 0;

    if (use_GMA)
      memset(erased_call_list, 0, ncall * sizeof(int));

    switch (disable_calls_method) {
    case 1:

      /* All in Range Method */
      /* Disable any Call or any KO which index is outside [first_call_index ,
       * last_call_index] */
      for (i = 0; i < ncall; i++) {
        /* Copy the arrays */
        new_call_type[i] = call_type[i];
        new_ex_date[i] = ex_date[i];
        new_set_date[i] = set_date[i];
        new_barrier[i] = barrier[i];
        new_bar_type[i] = bar_type[i];
        new_fees[i] = fees[i];
        new_pd_not_ref_fix[i] = pd_not_ref_fix[i];
        new_pd_not_ref_fix_fx[i] = pd_not_ref_fix_fx[i];
        if (use_not_opt_str) {
          new_pd_not_ref_nfxpts[i] = pd_not_ref_nfxpts[i];
          memcpy(new_pd_not_ref_fxpts[i], pd_not_ref_fxpts[i],
                 pd_not_ref_nfxpts[i] * sizeof(double));
          memcpy(new_pd_not_ref_cpn_at_pts[i], pd_not_ref_cpn_at_pts[i],
                 pd_not_ref_nfxpts[i] * sizeof(double));
          new_pd_not_ref_lin_xtrpl_l[i] = pd_not_ref_lin_xtrpl_l[i];
          new_pd_not_ref_lin_xtrpl_r[i] = pd_not_ref_lin_xtrpl_r[i];
        }

        if ((i < first_call_index) || (i > last_call_index)) {

          if (use_GMA && ex_date[i] >= today + eod_ex_flag) {
            erased_call_list[j] = 1;
            j++;
          }

          /* Disable the calls and KO */
          if (call_type[i] == 0) {
            /* Call Case */
            /* Disable by putting a huge fee */
            new_fees[i] += CPD_HUGE_FEES;
          } else if (call_type[i] == 1) {
            /* KO Case */
            /* Disable by putting a huge Barrier */
            if (bar_type[i]) {
              /* DO Case */
              new_barrier[i] = 0.0;
            } else {
              /* UO Case */
              new_barrier[i] += CPD_HUGE_BARRIER;
            }
          }
        }
      }
      break;

    case 2:

      /* All in Range Method 2 */
      /* Disable any Call or any KO which index belogs to the index list array
       */
      for (i = 0; i < ncall; i++) {
        /* Copy the arrays */
        new_call_type[i] = call_type[i];
        new_ex_date[i] = ex_date[i];
        new_set_date[i] = set_date[i];
        new_barrier[i] = barrier[i];
        new_bar_type[i] = bar_type[i];
        new_fees[i] = fees[i];
        new_pd_not_ref_fix[i] = pd_not_ref_fix[i];
        new_pd_not_ref_fix_fx[i] = pd_not_ref_fix_fx[i];
        if (use_not_opt_str) {
          new_pd_not_ref_nfxpts[i] = pd_not_ref_nfxpts[i];
          memcpy(new_pd_not_ref_fxpts[i], pd_not_ref_fxpts[i],
                 pd_not_ref_nfxpts[i] * sizeof(double));
          memcpy(new_pd_not_ref_cpn_at_pts[i], pd_not_ref_cpn_at_pts[i],
                 pd_not_ref_nfxpts[i] * sizeof(double));
          new_pd_not_ref_lin_xtrpl_l[i] = pd_not_ref_lin_xtrpl_l[i];
          new_pd_not_ref_lin_xtrpl_r[i] = pd_not_ref_lin_xtrpl_r[i];
        }
      }

      for (i = 0; i < iNbIndex_List; i++) {
        iErasedIndex = Index_List[i];

        if (use_GMA)
          erased_call_list[iErasedIndex] = 1;

        if ((iErasedIndex >= 0) && (iErasedIndex < ncall)) {
          /* Disable the calls and KO */
          if (call_type[iErasedIndex] == 0) {
            /* Call Case */
            /* Disable by putting a huge fee */
            new_fees[iErasedIndex] += CPD_HUGE_FEES;
          } else if (call_type[iErasedIndex] == 1) {
            /* KO Case */
            /* Disable by putting a huge Barrier */
            if (bar_type[iErasedIndex]) {
              /* DO Case */
              new_barrier[iErasedIndex] = 0.0;
            } else {
              /* UO Case */
              new_barrier[iErasedIndex] += CPD_HUGE_BARRIER;
            }
          }
        }
      }

      break;

    default:
      err = "disable_calls_method not recognized";
      goto FREE_RETURN;
      break;
    }

    new_ncall = ncall;

    new_pd_not_ref_fix[new_ncall] = pd_not_ref_fix[ncall];
    new_pd_not_ref_fix_fx[new_ncall] = pd_not_ref_fix_fx[ncall];
    if (use_not_opt_str &
        (new_ncall != ncall)) // copy final redemption information:
    {
      new_pd_not_ref_nfxpts[new_ncall] = pd_not_ref_nfxpts[ncall];
      memcpy(new_pd_not_ref_fxpts[new_ncall], pd_not_ref_fxpts[ncall],
             pd_not_ref_nfxpts[ncall] * sizeof(double));
      memcpy(new_pd_not_ref_cpn_at_pts[new_ncall], pd_not_ref_cpn_at_pts[ncall],
             pd_not_ref_nfxpts[ncall] * sizeof(double));
      new_pd_not_ref_lin_xtrpl_l[new_ncall] = pd_not_ref_lin_xtrpl_l[ncall];
      new_pd_not_ref_lin_xtrpl_r[new_ncall] = pd_not_ref_lin_xtrpl_r[ncall];
    }

    disable_calls_method = 0;
    erasing_call_done = 1;

    err = cpd_caller(
        today, use_calib, fx_spot, fx_spot_date, dom_calib, dom_und, dom_yc,
        dom_vc, dom_ref, dom_swap_freq, dom_swap_basis, dom_lam, for_calib,
        for_und, for_yc, for_vc, for_ref, for_swap_freq, for_swap_basis,
        for_lam, min_fact, max_fact, use_jumps, corr_times, correl_dom_for,
        correl_dom_fx, correl_for_fx, corr_n_times, cpd_dlm_params,
        get_ir_cash_vol, fx_mkt_vol_date, fx_mkt_vol, num_fx_mkt_vol,
        fx_mkt_smile_alpha, fx_mkt_smile_beta, fx_mkt_smile_rho,
        fx_mkt_smile_pi, use_sabr, use_3F_interp, smile_spec_type, fx3dund,
        start_date, fund_not, fund_ccy, fund_ccy_yc, fx_fund_dom,
        fx_fund_dom_spot_date, fund_ncpn, fund_fix, fund_start, fund_pay,
        fund_basis, fund_spr, fund_mrg, fund_fix_cpn, pd_not, pd_ncpn, pd_fix,
        pd_start, pd_pay, pd_basis, pd_alpha, pd_beta, pd_floored, pd_floor,
        pd_capped, pd_cap, pd_nfxpts, pd_fxpts, pd_cpn_at_pts, pd_lin_xtrpl_l,
        pd_lin_xtrpl_r, pd_fix_fx, new_pd_not_ref_fix, pd_not_ref_alpha,
        pd_not_ref_beta, pd_not_ref_floored, pd_not_ref_floor,
        pd_not_ref_capped, pd_not_ref_cap, new_pd_not_ref_nfxpts,
        new_pd_not_ref_fxpts, new_pd_not_ref_cpn_at_pts,
        new_pd_not_ref_lin_xtrpl_l, new_pd_not_ref_lin_xtrpl_r,
        new_pd_not_ref_fix_fx, new_call_type, new_ncall, pay_rec, new_ex_date,
        new_set_date, new_barrier, new_bar_type, new_fees,
        /* TARN */
        TARN_Do, 0.0, req_stp, req_pth, bar_smooth, do_pecs, forcetree,
        do_optim, force_optim, fx_bound, use_bound, do_infos, eod_fix_flag,
        eod_pay_flag, eod_ex_flag, dom_vol_shift, for_vol_shift, fx_vol_shift,
        exercised, ex_date_ex, ex_date_set, prune_calls, no_prune_years,
        prune_factor,

        disable_calls_method, first_call_index, last_call_index,

        iNbIndex_List, Index_List,

        erasing_call_done, erased_call_list,

        use_smile, alpha, beta,

        /*For the smile impact*/
        use_GMA, use_3f_optim_barrier, SSLstd, SSLniter, CopNsimul, CopDoPecs,
        CummulNpoints, CummulNstd, CummulPrecision, CummulLinear,
        PayoffFunction, fwdSmileVisuNstd, smileOtc, smileFee, fwdVolMethod,
        smileModel, BMpi, otc, FundingSpeedUp,
        /* Do not calculate all OTC*/
        nStart, oneOutOfN,
        // Fast MC
        FMC_do, FMC_precision, FMC_min_paths, fund_val, pd_val, call_val,
        call_stdev, smile_adjustment, optim_bar, GMA_Results, export_ts,
        und_exp);

    goto FREE_RETURN;
  }

  /*	If tierce funding */
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
    if (err)
      goto FREE_RETURN;
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

  // Prepare beta DLM smile calib:

  if (cpd_dlm_params->use_beta_dlm && cpd_dlm_params->calib_smile &&
      use_calib) {
    cpd_dlm_params->smile_settlmt_date = (long)(add_unit(
        fx_spot_date, (int)(12 * cpd_dlm_params->calib_smile_maturity + 0.5),
        SRT_MONTH, NO_BUSDAY_CONVENTION));

    i = 0;
    while (cpd_dlm_params->smile_settlmt_date > fx_mkt_vol_date[i] &&
           i < (num_fx_mkt_vol - 1))
      i++;

    if (!i || i == num_fx_mkt_vol - 1) {
      if (fabs(cpd_dlm_params->smile_settlmt_date - fx_mkt_vol_date[i]) <
          5) /* 5 days */
        cpd_dlm_params->smile_settlmt_date = fx_mkt_vol_date[i];
      else {
        err =
            "CPDcaller : smile calibration not included in the term structure";
        goto FREE_RETURN;
      }
    } else {
      if (fabs(cpd_dlm_params->smile_settlmt_date - fx_mkt_vol_date[i]) <
          fabs(cpd_dlm_params->smile_settlmt_date - fx_mkt_vol_date[i - 1])) {
        if (fabs(cpd_dlm_params->smile_settlmt_date - fx_mkt_vol_date[i]) <
            5) /* 5 days */
          cpd_dlm_params->smile_settlmt_date = fx_mkt_vol_date[i];
        else {
          err = "CPDcaller : smile calibration not included in the term "
                "structure";
          goto FREE_RETURN;
        }
      } else {
        if (fabs(cpd_dlm_params->smile_settlmt_date - fx_mkt_vol_date[i - 1]) <
            5) /* 5 days */
        {
          cpd_dlm_params->smile_settlmt_date = fx_mkt_vol_date[i - 1];
          i--;
        } else {
          err = "CPDcaller : smile calibration not included in the term "
                "structure";
          goto FREE_RETURN;
        }
      }
    }

    cpd_dlm_params->calib_smile_maturity =
        YEARS_IN_DAY * (add_unit(cpd_dlm_params->smile_settlmt_date, -2,
                                 SRT_BDAY, MODIFIED_SUCCEEDING) -
                        today);

    vol_at_smile = fx_mkt_vol[i];
    alpha_at_smile = fx_mkt_smile_alpha[i];
    beta_at_smile = fx_mkt_smile_beta[i];
    rho_at_smile = fx_mkt_smile_rho[i] + cpd_dlm_params->calib_smile_shift;

    forward_at_smile =
        fx_spot * swp_f_df(today, cpd_dlm_params->smile_settlmt_date, for_yc) /
        swp_f_df(today, cpd_dlm_params->smile_settlmt_date, dom_yc) *
        swp_f_df(today, fx_spot_date, dom_yc) /
        swp_f_df(today, fx_spot_date, for_yc);

    cpd_dlm_params->smile_strikes[0] =
        forward_at_smile *
        exp(cpd_dlm_params->calib_smile_std_up * vol_at_smile *
            sqrt(cpd_dlm_params->calib_smile_maturity));
    cpd_dlm_params->smile_strikes[1] =
        forward_at_smile *
        exp(cpd_dlm_params->calib_smile_std_down * vol_at_smile *
            sqrt(cpd_dlm_params->calib_smile_maturity));

    srt_f_optsarbvol(forward_at_smile, cpd_dlm_params->smile_strikes[0],
                     cpd_dlm_params->calib_smile_maturity, vol_at_smile,
                     alpha_at_smile, beta_at_smile, rho_at_smile, SRT_LOGNORMAL,
                     SRT_LOGNORMAL, &(cpd_dlm_params->smile_bs_vols[0]));

    srt_f_optsarbvol(forward_at_smile, cpd_dlm_params->smile_strikes[1],
                     cpd_dlm_params->calib_smile_maturity, vol_at_smile,
                     alpha_at_smile, beta_at_smile, rho_at_smile, SRT_LOGNORMAL,
                     SRT_LOGNORMAL, &(cpd_dlm_params->smile_bs_vols[1]));
  }

  /*	Initialise structures */

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
  free_struct = 0;
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
      ncall, pay_rec, ex_date, set_date, barrier, bar_type, fees, TARN_Do,
      req_stp, req_pth, do_pecs, forcetree, do_optim, force_optim, fx_bound,
      use_bound, bar_smooth, eod_fix_flag, eod_ex_flag, cpd, und, &call_feat,
      tree_arg, mc_arg, dom_vol_shift, for_vol_shift, fx_vol_shift, skip_fill);

  if (err)
    goto FREE_RETURN;

  pd_fut_ncpn = cpd->pd_leg->num_cpn;
  fund_fut_ncpn = cpd->fund_leg->num_cpn;
  ex_idx = ncall;

  //	If exercised
  if (exercised) {
    // Consider only the coupons which start before the ex date

    for (i = 0; i < cpd->pd_leg->num_cpn &&
                cpd->pd_leg->cpn[i].start_date < ex_date_ex;
         i++)
      ;
    pd_fut_ncpn = i;

    for (i = 0; i < cpd->fund_leg->num_cpn &&
                cpd->fund_leg->cpn[i].start_date < ex_date_ex;
         i++)
      ;
    fund_fut_ncpn = i;

    for (ex_idx = 0; ex_idx < ncall && set_date[ex_idx] < ex_date_ex; ex_idx++)
      ;

    if (ex_idx >= ncall) {
      err = serror("Deal cannot be exercised after the last exercise date");
      goto FREE_RETURN;
    }

    call_feat = 0;
    ncall = 0;

    // Modify redemption coupon for the exercised case:
    cpn = &(cpd->pd_leg->not_ref);

    // Dates:
    cpn->start_date = cpn->pay_date = ex_date_set;
    cpn->fx_fix_date = pd_not_ref_fix[ex_idx];
    cpn->fx_val_date =
        add_unit(cpn->fx_fix_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    // Times:
    cpn->start_time = cpn->pay_time = (ex_date_set - today) * YEARS_IN_DAY;
    cpn->fx_fix_time = (cpn->fx_fix_date - today) * YEARS_IN_DAY;
    cpn->fx_val_time = (cpn->fx_val_date - today) * YEARS_IN_DAY;

    // Coupon:
    if (!use_not_opt_str) {
      cpn->alpha = pd_not;
      cpn->beta = 0.0;
      cpn->floored = cpn->capped = 0;
    } else {
      cpn->nstrikes = pd_not_num_strikes[ex_idx];
      cpn->wcst = pd_not * pd_not_wcst[ex_idx];
      cpn->wspot = pd_not * pd_not_wspot[ex_idx];
      if (pd_not_num_strikes[ex_idx] > 0) {
        if (cpn->strikes)
          free(cpn->strikes);
        if (cpn->weights)
          free(cpn->weights);
        cpn->strikes = calloc(pd_not_num_strikes[ex_idx], sizeof(double));
        cpn->weights = calloc(pd_not_num_strikes[ex_idx], sizeof(double));
        if (!cpn->strikes || !cpn->weights) {
          err = serror("Memory failure in cpd_caller");
          goto FREE_RETURN;
        }
        memcpy(cpn->strikes, pd_not_strikes[ex_idx],
               pd_not_num_strikes[ex_idx] * sizeof(double));
        for (i = 0; i < pd_not_num_strikes[ex_idx]; i++)
          cpn->weights[i] = pd_not * pd_not_weights[ex_idx][i];
      }
    }
  }

  /****************************/
  /* Switch to the smile case */
  /****************************/
  if (use_smile) {
    /* Check if exercised */
    if (exercised) {
      err = "Exercised flag temporary not allowed with UseSmile ";
      goto FREE_RETURN;
    }

    err = call_GRFN_alphabetaModel(
        /* Input */
        und, /* The CPD underlying */
        cpd, /* The CPD structure */

        use_calib,          /*	0: use fx3dund        , 1: calibrate */
        fx3dund, dom_calib, /*	Calibrate domestic underlying */
        dom_und,   /*	If no        , domestic underlying name to be used */
        dom_yc,    /*	Domestic yc */
        for_calib, /*	Same for foreign */
        for_und, for_yc,

        /*	The structure */
        start_date,

        /*		funding */
        fund_not,    /*	Notional */
        for_fund,    /*	1: if funding has been transform into domestic */
        eq_final_ex, /*	equivalent final exchange use only if for_fund = 1 */
        eq_init_ex,  /*	equivalent init exchange use only if for_fund = 1 */
        fund_start_date, fund_ncpn, /*	Number of coupons */
        fund_fix,                   /*	Fixing dates */
        fund_start,                 /*	Start dates */
        fund_pay,                   /*	Pay dates */
        fund_basis,                 /*	Basis */
        fund_mrg,                   /*	Margins */
        fund_fix_cpn,               /*	Past coupon fixing if relevant        ,
                                                                              includes
                                 spr        , but not mrg        , cvg and notional */
        /*		pd */
        pd_not,   /*	Notional */
        pd_ncpn,  /*	Number of coupons */
        pd_fix,   /*	Fx fixing dates */
        pd_start, /*	Start dates */
        pd_pay,   /*	Pay dates */
        pd_basis, /*	Basis */
        pd_alpha, /*	Coupon = alpha + beta * fx [capped        , floored] */
        pd_beta, pd_floored, pd_floor, pd_capped, pd_cap,
        pd_fix_fx, /*	Past Fx fixing if relevant */

        /* Call */
        ncall,   /*	Number of calls */
        ex_date, /*	Call dates */
        barrier, /*	in case of a pure KO or a Callable KO */

        /*	EOD Fixing Flag */
        eod_fix_flag, /*	0: I        , 1: E */
        /*	EOD Payment Flag */
        eod_pay_flag, /*	0: I        , 1: E */

        alpha, beta,

        req_stp,

        /* Fx vol from the market */
        fx_mkt_vol_date, fx_mkt_vol,

        num_fx_mkt_vol, fx_spot, /* Spot Fx */

        /* OutPut */
        fund_val, /*	Value of the funding leg */
        pd_val,   /*	Value of the Power Dual leg */
        call_val, /*	Value of the callable feature */
        call_stdev /*	Standard deviation of the call if applicable */);

    /* export the ts if asked */
    if (export_ts) {
      cpd_copy_und(und, und_exp);
    }

    goto FREE_RETURN;
  }

  /****************************/
  /* End Switch to smile Case */
  /****************************/
  if (call_feat == 1 || (call_feat == 2 && cpd_dlm_params->use_beta_dlm)) {
    partial_pv = 1;
    partial_fund_idx = cpd->call[0].fund_idx;
    partial_pd_idx = cpd->call[0].pd_idx;
  } else {
    partial_pv = 0;
  }

  /*	1)	Value funding leg */

  /*	Funding: domestic/foreign */
  if (cpd->fund_leg->dom_for == 0) {
    fund_yc = (char *)(und->dom_yc);
    fund_mult = 1.0;
  } else {
    fund_yc = (char *)(und->for_yc);
    fund_mult = und->spot_fx;
  }

  fund_leg = partial_fund_leg = 0.0;

  /*	Coupons: spread + margin */
  for (i = 0; i < fund_fut_ncpn; i++) {
    temp = swp_f_df(und->today, cpd->fund_leg->cpn[i].pay_date, fund_yc) *
           cpd->fund_leg->cpn[i].cxxpn;
    fund_leg += temp;

    if (partial_pv && i >= partial_fund_idx) {
      partial_fund_leg += temp;
    }
  }

  /*	Cash libor */
  if (fund_fut_ncpn > 0) {
    if (for_fund) {
      fund_leg += swp_f_df(und->today, cpd->fund_leg->cpn[0].start_date,
                           Copyfund_ccy_yc) *
                  thirdCcyfundNot * fx_fund_dom *
                  swp_f_df(und->today, fx_fund_dom_spot_date, dom_yc) /
                  swp_f_df(und->today, fx_fund_dom_spot_date, Copyfund_ccy_yc);

      // Exerciced Case
      // Force the Notional to be paid at ex_settle_date
      if (exercised &&
          cpd->fund_leg->cpn[fund_fut_ncpn - 1].pay_date != ex_date_set) {
        // Third currency case:
        if (cpd->fund_leg->cpn[fund_fut_ncpn - 1].pay_date >=
            today + eod_pay_flag) {
          fund_leg -=
              swp_f_df(und->today,
                       cpd->fund_leg->cpn[fund_fut_ncpn - 1].pay_date,
                       Copyfund_ccy_yc) *
              thirdCcyfundNot * fx_fund_dom *
              swp_f_df(und->today, fx_fund_dom_spot_date, dom_yc) /
              swp_f_df(und->today, fx_fund_dom_spot_date, Copyfund_ccy_yc);
        }

        if (ex_date_set >= today + eod_pay_flag) {
          fund_leg +=
              swp_f_df(und->today, ex_date_set, Copyfund_ccy_yc) *
              thirdCcyfundNot * fx_fund_dom *
              swp_f_df(und->today, fx_fund_dom_spot_date, dom_yc) /
              swp_f_df(und->today, fx_fund_dom_spot_date, Copyfund_ccy_yc);
        }
      }
    } else {
      fund_leg +=
          swp_f_df(und->today, cpd->fund_leg->cpn[0].start_date, fund_yc) *
          cpd->fund_leg->notional;

      if (exercised &&
          cpd->fund_leg->cpn[fund_fut_ncpn - 1].pay_date != ex_date_set) {
        if (cpd->fund_leg->cpn[fund_fut_ncpn - 1].pay_date >=
            today + eod_pay_flag) {
          fund_leg -= swp_f_df(und->today,
                               cpd->fund_leg->cpn[fund_fut_ncpn - 1].pay_date,
                               fund_yc) *
                      cpd->fund_leg->notional;
        }

        if (ex_date_set >= today + eod_pay_flag) {
          fund_leg += swp_f_df(und->today, ex_date_set, fund_yc) *
                      cpd->fund_leg->notional;
        }
      }
    }

    if (partial_pv) {
      partial_fund_leg +=
          swp_f_df(und->today, cpd->fund_leg->cpn[partial_fund_idx].start_date,
                   fund_yc) *
          cpd->fund_leg->notional;
    }
  } else {
    fin_not_date = (exercised ? ex_date_set : fund_pay[fund_ncpn - 1]);

    if (fin_not_date >= today + eod_pay_flag) {

      if (for_fund) {
        fund_leg +=
            swp_f_df(und->today, fin_not_date, Copyfund_ccy_yc) *
            thirdCcyfundNot * fx_fund_dom *
            swp_f_df(und->today, fx_fund_dom_spot_date, dom_yc) /
            swp_f_df(und->today, fx_fund_dom_spot_date, Copyfund_ccy_yc);

      } else {
        fund_leg += swp_f_df(und->today, fin_not_date, fund_yc) *
                    cpd->fund_leg->notional;
      }
    }
  }

  /*	Initial notional if applicable */
  if (start_date >= und->today + eod_pay_flag) {
    if (for_fund) {
      fund_leg -= swp_f_df(und->today, start_date, Copyfund_ccy_yc) *
                  thirdCcyfundNot * fx_fund_dom *
                  swp_f_df(und->today, fx_fund_dom_spot_date, dom_yc) /
                  swp_f_df(und->today, fx_fund_dom_spot_date, Copyfund_ccy_yc);

    } else {
      fund_leg -= swp_f_df(und->today, start_date, fund_yc) * fund_not;
    }
  }

  if (partial_pv) {
    partial_fund_leg -= swp_f_df(und->today, cpd->call[0].set_date, fund_yc) *
                        cpd->call[0].fund_not_amt;
  }

  /*	PV of coupons fixed in the past and not yet paid */
  i = 0;
  while (i < fund_ncpn && fund_fix[i] < und->today + eod_fix_flag) {
    if (fund_pay[i] >= und->today + eod_pay_flag &&
        (!exercised || fund_start[i] < ex_date_ex)) {
      err = interp_basis(fund_basis[i], &bas);
      if (err)
        goto FREE_RETURN;

      fund_leg += (fund_fix_cpn[i] + fund_mrg[i]) *
                  coverage(fund_start[i], fund_pay[i], bas) * fund_not *
                  swp_f_df(und->today, fund_pay[i], fund_yc);
    }
    i++;
  }

  /*	In order to get the PV in domestic units */
  fund_leg *= fund_mult;
  partial_fund_leg *= fund_mult;

  /*	2)	Value pd leg */

  pd_leg = partial_pd_leg = 0.0;

  if (use_sabr) {
    adj_fees = (double *)calloc(pd_fut_ncpn, sizeof(double));
    if (!adj_fees) {
      err = serror("Memory failure in cpd_caller");
      goto FREE_RETURN;
    }
  }

  /*	Coupons */

  if (cpd_dlm_params->use_beta_dlm) {
    /* allocate the instrument */
    Inst = calloc(1, sizeof(FxBetaDLM_FxOptInst));
    InstConst = calloc(1, sizeof(FxBetaDLM_InstPrecalc));

    if (!Inst || !InstConst) {
      err = "Memory allocation faillure in cpd_caller";
      goto FREE_RETURN;
    }

    err = FxBetaDLM_Allocate_FxOptInst(0, 1, Inst);

    if (err)
      goto FREE_RETURN;
  }

  dlm_pd_leg = pd_leg;
  dlm_partial_pd_leg = partial_pd_leg;

  for (i = 0; i < pd_fut_ncpn; i++) {
    cpn = cpd->pd_leg->cpn + i;

    /*	Discount */
    df = swp_f_df(und->today, cpn->pay_date, und->dom_yc);

    /*	Fwd fx */

    fwd = und->spot_fx * swp_f_df(und->today, cpn->fx_val_date, und->for_yc) /
          swp_f_df(und->today, cpn->fx_val_date, und->dom_yc);

    err = Fx3DtsFwdPayAdjustment_corr(
        0.0, cpn->fx_val_time, cpn->fx_val_time, cpn->pay_time,
        cpn->fx_fix_time, und->sigma_time_rates, und->sigma_n_rates,
        und->sigma_dom, und->lda_dom, und->sigma_for, und->lda_for,
        und->sigma_time_fx, und->sigma_fx, und->sigma_n_fx, und->corr_times,
        und->correl_dom_for, und->correl_dom_fx, und->correl_for_fx,
        und->corr_n_times, &adj);
    if (err) {
      goto FREE_RETURN;
    }

    fwd *= exp(adj);

    /*	Vol */

    if (!use_3F_interp) {
      if (cpn->fx_val_date <= fx_mkt_vol_date[0])
        std = fx_mkt_vol[0];
      else if (cpn->fx_val_date >= fx_mkt_vol_date[num_fx_mkt_vol - 1])
        std = fx_mkt_vol[num_fx_mkt_vol - 1];
      else {
        for (j = 0; cpn->fx_val_date >= fx_mkt_vol_date[j + 1]; j++)
          ;
        interp_coef = (double)(cpn->fx_val_date - fx_mkt_vol_date[j]) /
                      (double)(fx_mkt_vol_date[j + 1] - fx_mkt_vol_date[j]);

        std = interp_coef * fx_mkt_vol[j + 1] +
              (1.0 - interp_coef) * fx_mkt_vol[j];
      }
    } else {
      err = Fx3DtsImpliedVol_corr(
          cpn->fx_val_time, 0.0, cpn->fx_fix_time, und->sigma_time_rates,
          und->sigma_n_rates, und->sigma_dom, und->lda_dom, und->sigma_for,
          und->lda_for, und->sigma_time_fx, und->sigma_fx, und->sigma_n_fx,
          und->corr_times, und->correl_dom_for, und->correl_dom_fx,
          und->correl_for_fx, und->corr_n_times, &std);
      if (err) {
        goto FREE_RETURN;
      }
    }

    smile_params->sigma = std;

    std *= sqrt(cpn->fx_fix_time);
    half_std = 0.5 * std;

    // get the sabr params at the cpn->fx_val_date
    if (use_sabr)
      err = cpd_vol_get_smile_params(0, cpn->fx_val_date, 0.0, smile_mkt,
                                     smile_params);

    if (cpd_dlm_params->use_beta_dlm) {
      InstStrike = fwd / 100000.0;

      err = FxBetaDLM_Setup_FxOptInst(cpn->fx_val_date, cpn->pay_date, today,
                                      cpn->fx_fix_date, 1, &InstStrike, NULL,
                                      SRT_CALL, und->model, Inst);

      if (err)
        goto FREE_RETURN;

      /* Remove discounting */
      Inst->dDfPayDom = 1.0;

      err = FxBetaDLM_Calculate_AllConst(Inst, und->model, und->num_params,
                                         InstConst);

      if (err)
        goto FREE_RETURN;

      err = FxBetaDLM_Update_InstPrecalc_FromMoment(Inst, und->model,
                                                    und->num_params, InstConst);

      if (err)
        goto FREE_RETURN;

      err = FxBetaDLM_Price_FxOptInst(Inst, und->model, InstConst, und->hermite,
                                      und->num_params, &dlm_fwd);

      if (err)
        goto FREE_RETURN;
    }

    if (use_cpn_opt_str) {
      opt_string = smile_opt_string = 0.0;

      for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++) {
        if (fabs(cpn->weights[str_idx]) > 1e-16) {
          if (std > 1e-16 && cpn->strikes[str_idx] > 1e-16)
            type = 3;
          else
            type = 1;

          opt_string +=
              cpn->weights[str_idx] *
              OPT_VAL_MACRO(type, fwd, cpn->strikes[str_idx], std, half_std);

          if (use_sabr) {
            if (type >= 3) {
              cpd_vol_get_vol(fwd, cpn->fx_fix_time, cpn->strikes[str_idx],
                              smile_params, SRT_LOGNORMAL, &smile_std);

              if (err)
                goto FREE_RETURN;

              smile_std *= sqrt(cpn->fx_fix_time);
              half_smile_std = 0.5 * smile_std;
            }

            smile_opt_string += cpn->weights[str_idx] *
                                OPT_VAL_MACRO(type, fwd, cpn->strikes[str_idx],
                                              smile_std, half_smile_std);
          }
        }
      }
      if (use_sabr)
        adj_fees[i] = df * (smile_opt_string - opt_string);

      /*	Coupon pv */
      temp = df * (cpn->wcst + cpn->wspot * fwd + opt_string);
      pd_leg += temp;
      if (use_sabr)
        pd_leg += adj_fees[i];

      if (partial_pv && i >= partial_pd_idx) {
        partial_pd_leg += temp;
      }

    } else {
      /*	Floor */

      if (cpn->floored && fabs(cpn->beta) > 1.0e-16 && cpn->floor > -1.0e-16) {
        str = (cpn->floor - cpn->alpha) / cpn->beta;

        if (cpn->beta > 0.0) {
          if (std > 1.0e-16 && str > 1.0e-16) {
            type = 4; /*	Put */
          } else {
            type = 2; /*	Put IV */
          }

          InstCallType = SRT_PUT;
        } else {
          if (std > 1.0e-16 && str > 1.0e-16) {
            type = 3; /*	Call */
          } else {
            type = 1; /*	Call IV */
          }

          InstCallType = SRT_CALL;
        }
      } else {
        if (cpn->floor <= -1.0e-16) {
          err = "Coupon floor must be positive";
          goto FREE_RETURN;
        } else {
          type = 0; /*	No floor */
        }
      }

      floor = OPT_VAL_MACRO(type, fwd, str, std, half_std);

      if (use_sabr) {
        if (type >= 3) {
          cpd_vol_get_vol(fwd, cpn->fx_fix_time, str, smile_params,
                          SRT_LOGNORMAL, &smile_std);

          if (err)
            goto FREE_RETURN;

          smile_std *= sqrt(cpn->fx_fix_time);
          half_smile_std = 0.5 * smile_std;
        }

        smile_floor = OPT_VAL_MACRO(type, fwd, str, smile_std, half_smile_std);
      }

      if (cpd_dlm_params->use_beta_dlm) {
        if (type >= 3) {
          /* Setup the instrument */

          InstStrike = max(str, 1.0E-16);
          Inst->dStrike[0] = InstStrike;
          Inst->sCallPutType = InstCallType;

          if (err)
            goto FREE_RETURN;

          err = FxBetaDLM_Price_FxOptInst(Inst, und->model, InstConst,
                                          und->hermite, und->num_params,
                                          &dlm_floor);

          if (err)
            goto FREE_RETURN;
        } else {
          dlm_floor = OPT_VAL_MACRO(type, dlm_fwd, str, std, half_std);
        }
      }

      /*	Cap */

      if (cpn->capped && fabs(cpn->beta) > 1.0e-16 && cpn->cap > -1.0e-16) {
        str = (cpn->cap - cpn->alpha) / cpn->beta;

        if (cpn->beta > 0.0) {
          if (std > 1.0e-16 && str > 1.0e-16) {
            type = 3; /*	Call */
          } else {
            type = 1; /*	Call IV */
          }

          InstCallType = SRT_CALL;
        } else {
          if (std > 1.0e-16 && str > 1.0e-16) {
            type = 4; /*	Put */
          } else {
            type = 2; /*	Put IV */
          }

          InstCallType = SRT_PUT;
        }
      } else {
        if (cpn->cap <= -1.0e-16) {
          err = "Coupon cap must be positive";
          goto FREE_RETURN;
        }

        type = 0; /*	No cap */
      }

      cap = OPT_VAL_MACRO(type, fwd, str, std, half_std);

      if (use_sabr) {
        if (type >= 3) {
          cpd_vol_get_vol(fwd, cpn->fx_fix_time, str, smile_params,
                          SRT_LOGNORMAL, &smile_std);

          if (err)
            goto FREE_RETURN;

          smile_std *= sqrt(cpn->fx_fix_time);
          half_smile_std = 0.5 * smile_std;
        }

        smile_cap = OPT_VAL_MACRO(type, fwd, str, smile_std, half_smile_std);
      }

      if (cpd_dlm_params->use_beta_dlm) {
        if (type >= 3) {
          /* Setup the instrument */

          InstStrike = max(str, 1.0E-16);
          Inst->dStrike[0] = InstStrike;
          Inst->sCallPutType = InstCallType;

          if (err)
            goto FREE_RETURN;

          err = FxBetaDLM_Price_FxOptInst(Inst, und->model, InstConst,
                                          und->hermite, und->num_params,
                                          &dlm_cap);

          if (err)
            goto FREE_RETURN;
        } else {
          dlm_cap = OPT_VAL_MACRO(type, dlm_fwd, str, std, half_std);
        }
      }

      if (use_sabr) {
        if (cpd_dlm_params->use_beta_dlm) {
          adj_fees[i] = df * fabs(cpn->beta) *
                        (smile_floor - dlm_floor - smile_cap + dlm_cap);
        } else {
          adj_fees[i] =
              df * fabs(cpn->beta) * (smile_floor - floor - smile_cap + cap);
        }
      }

      /*	Coupon pv */
      temp =
          df * (cpn->alpha + cpn->beta * fwd + fabs(cpn->beta) * (floor - cap));
      pd_leg += temp;
      if (use_sabr)
        pd_leg += adj_fees[i];

      if (partial_pv && i >= partial_pd_idx) {
        partial_pd_leg += temp;
      }

      if (cpd_dlm_params->use_beta_dlm) {
        dlm_temp = df * (cpn->alpha + cpn->beta * dlm_fwd +
                         fabs(cpn->beta) * (dlm_floor - dlm_cap));
        dlm_pd_leg += dlm_temp;

        if (use_sabr)
          dlm_pd_leg += adj_fees[i];

        if (partial_pv && i >= partial_pd_idx) {
          dlm_partial_pd_leg += dlm_temp;
        }
      }
    }
  }

  /*	Final notional */
  cpn = &(cpd->pd_leg->not_ref);

  if (cpn->pay_date >= und->today + eod_pay_flag) {
    /*	Discount */
    df = swp_f_df(und->today, cpn->pay_date, und->dom_yc);

    if (cpn->fx_fix_date >= und->today + eod_fix_flag) {
      /*	Fwd fx */

      fwd = und->spot_fx * swp_f_df(und->today, cpn->fx_val_date, und->for_yc) /
            swp_f_df(und->today, cpn->fx_val_date, und->dom_yc);

      err = Fx3DtsFwdPayAdjustment_corr(
          0.0, cpn->fx_val_time, cpn->fx_val_time, cpn->pay_time,
          cpn->fx_fix_time, und->sigma_time_rates, und->sigma_n_rates,
          und->sigma_dom, und->lda_dom, und->sigma_for, und->lda_for,
          und->sigma_time_fx, und->sigma_fx, und->sigma_n_fx, und->corr_times,
          und->correl_dom_for, und->correl_dom_fx, und->correl_for_fx,
          und->corr_n_times, &adj);
      if (err) {
        goto FREE_RETURN;
      }

      fwd *= exp(adj);

      /*	Vol */

      if (!use_3F_interp) {
        if (cpn->fx_val_date <= fx_mkt_vol_date[0])
          std = fx_mkt_vol[0];
        else if (cpn->fx_val_date >= fx_mkt_vol_date[num_fx_mkt_vol - 1])
          std = fx_mkt_vol[num_fx_mkt_vol - 1];
        else {
          for (j = 0; cpn->fx_val_date >= fx_mkt_vol_date[j + 1]; j++)
            ;
          interp_coef = (double)(cpn->fx_val_date - fx_mkt_vol_date[j]) /
                        (double)(fx_mkt_vol_date[j + 1] - fx_mkt_vol_date[j]);

          std = interp_coef * fx_mkt_vol[j + 1] +
                (1.0 - interp_coef) * fx_mkt_vol[j];
        }
      } else {
        err = Fx3DtsImpliedVol_corr(
            cpn->fx_val_time, 0.0, cpn->fx_fix_time, und->sigma_time_rates,
            und->sigma_n_rates, und->sigma_dom, und->lda_dom, und->sigma_for,
            und->lda_for, und->sigma_time_fx, und->sigma_fx, und->sigma_n_fx,
            und->corr_times, und->correl_dom_for, und->correl_dom_fx,
            und->correl_for_fx, und->corr_n_times, &std);
        if (err) {
          goto FREE_RETURN;
        }
      }

      smile_params->sigma = std;

      std *= sqrt(cpn->fx_fix_time);
      half_std = 0.5 * std;

      // get the sabr params at the cpn->fx_val_date
      if (use_sabr)
        err = cpd_vol_get_smile_params(0, cpn->fx_val_date, 0.0, smile_mkt,
                                       smile_params);

      if (cpd_dlm_params->use_beta_dlm) {
        InstStrike = fwd / 10000.0;

        err = FxBetaDLM_Setup_FxOptInst(cpn->fx_val_date, cpn->pay_date, today,
                                        cpn->fx_fix_date, 1, &InstStrike, NULL,
                                        SRT_CALL, und->model, Inst);

        if (err)
          goto FREE_RETURN;

        /* Remove discounting */
        Inst->dDfPayDom = 1.0;

        err = FxBetaDLM_Calculate_AllConst(Inst, und->model, und->num_params,
                                           InstConst);

        if (err)
          goto FREE_RETURN;

        err = FxBetaDLM_Update_InstPrecalc_FromMoment(
            Inst, und->model, und->num_params, InstConst);

        if (err)
          goto FREE_RETURN;

        err =
            FxBetaDLM_Price_FxOptInst(Inst, und->model, InstConst, und->hermite,
                                      und->num_params, &dlm_fwd);

        if (err)
          goto FREE_RETURN;
      }

      if (use_not_opt_str) {
        opt_string = smile_opt_string = 0.0;

        for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++) {
          if (fabs(cpn->weights[str_idx]) > 1e-16) {
            if (std > 1e-16 && cpn->strikes[str_idx] > 1e-16)
              type = 3;
            else
              type = 1;

            opt_string +=
                cpn->weights[str_idx] *
                OPT_VAL_MACRO(type, fwd, cpn->strikes[str_idx], std, half_std);

            if (use_sabr) {
              if (type >= 3) {
                cpd_vol_get_vol(fwd, cpn->fx_fix_time, cpn->strikes[str_idx],
                                smile_params, SRT_LOGNORMAL, &smile_std);

                if (err)
                  goto FREE_RETURN;

                smile_std *= sqrt(cpn->fx_fix_time);
                half_smile_std = 0.5 * smile_std;
              }

              smile_opt_string +=
                  cpn->weights[str_idx] *
                  OPT_VAL_MACRO(type, fwd, cpn->strikes[str_idx], smile_std,
                                half_smile_std);
            }
          }
        }
        if (use_sabr)
          adj_fee_not = df * (smile_opt_string - opt_string);

        /*	Coupon pv */
        temp = df * (cpn->wcst + cpn->wspot * fwd + opt_string);
        pd_leg += temp;
        if (use_sabr)
          pd_leg += adj_fee_not;

        if (partial_pv) {
          partial_pd_leg += temp;
        }
      } else {
        /*	Floor */

        if (cpn->floored && fabs(cpn->beta) > 1.0e-16 &&
            cpn->floor > -1.0e-16) {
          str = (cpn->floor - cpn->alpha) / cpn->beta;

          if (cpn->beta > 0.0) {
            if (std > 1.0e-16 && str > 1.0e-16) {
              type = 4; /*	Put */
            } else {
              type = 2; /*	Put IV */
            }

            InstCallType = SRT_PUT;
          } else {
            if (std > 1.0e-16 && str > 1.0e-16) {
              type = 3; /*	Call */
            } else {
              type = 1; /*	Call IV */
            }

            InstCallType = SRT_CALL;
          }
        } else {
          type = 0; /*	No floor */
        }

        floor = OPT_VAL_MACRO(type, fwd, str, std, half_std);

        if (use_sabr) {
          if (type >= 3) {
            cpd_vol_get_vol(fwd, cpn->fx_fix_time, str, smile_params,
                            SRT_LOGNORMAL, &smile_std);

            if (err)
              goto FREE_RETURN;

            smile_std *= sqrt(cpn->fx_fix_time);
            half_smile_std = 0.5 * smile_std;
          }

          smile_floor =
              OPT_VAL_MACRO(type, fwd, str, smile_std, half_smile_std);
        }

        if (cpd_dlm_params->use_beta_dlm) {
          if (type >= 3) {
            /* Setup the instrument */

            InstStrike = max(str, 1.0E-16);
            Inst->dStrike[0] = InstStrike;
            Inst->sCallPutType = InstCallType;

            if (err)
              goto FREE_RETURN;

            err = FxBetaDLM_Price_FxOptInst(Inst, und->model, InstConst,
                                            und->hermite, und->num_params,
                                            &dlm_floor);

            if (err)
              goto FREE_RETURN;
          } else {
            dlm_floor = OPT_VAL_MACRO(type, dlm_fwd, str, std, half_std);
          }
        }

        /*	Cap */

        if (cpn->capped && fabs(cpn->beta) > 1.0e-16 && cpn->cap > -1.0e-16) {
          str = (cpn->cap - cpn->alpha) / cpn->beta;

          if (cpn->beta > 0.0) {
            if (std > 1.0e-16 && str > 1.0e-16) {
              type = 3; /*	Call */
            } else {
              type = 1; /*	Call IV */
            }

            InstCallType = SRT_CALL;
          } else {
            if (std > 1.0e-16 && str > 1.0e-16) {
              type = 4; /*	Put */
            } else {
              type = 2; /*	Put IV */
            }

            InstCallType = SRT_PUT;
          }
        } else {
          type = 0; /*	No cap */
        }

        cap = OPT_VAL_MACRO(type, fwd, str, std, half_std);

        if (use_sabr) {
          if (type >= 3) {
            cpd_vol_get_vol(fwd, cpn->fx_fix_time, str, smile_params,
                            SRT_LOGNORMAL, &smile_std);

            if (err)
              goto FREE_RETURN;

            smile_std *= sqrt(cpn->fx_fix_time);
            half_smile_std = 0.5 * smile_std;
          }

          smile_cap = OPT_VAL_MACRO(type, fwd, str, smile_std, half_smile_std);
        }

        if (cpd_dlm_params->use_beta_dlm) {
          if (type >= 3) {
            /* Setup the instrument */

            InstStrike = max(str, 1.0E-16);
            Inst->dStrike[0] = InstStrike;
            Inst->sCallPutType = InstCallType;

            if (err)
              goto FREE_RETURN;

            err = FxBetaDLM_Price_FxOptInst(Inst, und->model, InstConst,
                                            und->hermite, und->num_params,
                                            &dlm_cap);

            if (err)
              goto FREE_RETURN;
          } else {
            dlm_cap = OPT_VAL_MACRO(type, dlm_fwd, str, std, half_std);
          }
        }

        if (use_sabr) {
          if (cpd_dlm_params->use_beta_dlm) {
            adj_fee_not = df * fabs(cpn->beta) *
                          (smile_floor - dlm_floor - smile_cap + dlm_cap);
          } else {
            adj_fee_not =
                df * fabs(cpn->beta) * (smile_floor - floor - smile_cap + cap);
          }
        }

        /*	Coupon pv */
        temp = df *
               (cpn->alpha + cpn->beta * fwd + fabs(cpn->beta) * (floor - cap));

        pd_leg += temp;
        if (use_sabr)
          pd_leg += adj_fee_not;

        if (partial_pv) {
          partial_pd_leg += temp;
        }

        if (cpd_dlm_params->use_beta_dlm) {
          dlm_temp = df * (cpn->alpha + cpn->beta * dlm_fwd +
                           fabs(cpn->beta) * (dlm_floor - dlm_cap));
          dlm_pd_leg += dlm_temp;

          if (use_sabr)
            dlm_pd_leg += adj_fee_not;

          if (partial_pv) {
            dlm_partial_pd_leg += dlm_temp;
          }
        }
      }
    }    // if (cpn->fx_fix_date >= und->today + eod_fix_flag)
    else // if final notional is already fixed but not yet paid
    {
      if (use_not_opt_str) {
        opt_string = 0.0;

        for (str_idx = 0; str_idx < cpn->nstrikes; str_idx++) {
          if (fabs(cpn->weights[str_idx]) > 1e-16 &&
              pd_not_ref_fix_fx[ex_idx] - cpn->strikes[str_idx] > 0.0)
            opt_string += cpn->weights[str_idx] *
                          (pd_not_ref_fix_fx[ex_idx] - cpn->strikes[str_idx]);
        }

        temp =
            (cpn->wcst + cpn->wspot * pd_not_ref_fix_fx[ex_idx] + opt_string) *
            df;
      } else {
        fixing = cpn->alpha + cpn->beta * pd_not_ref_fix_fx[ex_idx];

        if (cpn->floored && fabs(cpn->beta) > 1.0e-16 &&
            cpn->floor > -1.0e-16) {
          fixing = max(fixing, cpn->floor);
        } else {
          if (cpn->floor <= -1.0e-16) {
            err = serror("Floor must be positive");
            goto FREE_RETURN;
          }
        }

        if (cpn->capped && fabs(cpn->beta) > 1.0e-16 && cpn->cap > -1.0e-16) {
          fixing = min(fixing, cpn->cap);
        } else {
          if (cpn->cap <= -1.0e-16) {
            err = serror("Cap must be positive");
            goto FREE_RETURN;
          }
        }
        temp = fixing * df;
      }

      pd_leg += temp;
      dlm_pd_leg += temp;
    }
  }

  if (cpd_dlm_params->use_beta_dlm) {
    saved_pd_leg = pd_leg;
    saved_partial_pd_leg = partial_pd_leg;
  }

  /*	Initial notional if applicable */
  if (start_date >= und->today + eod_pay_flag) {
    pd_leg -= swp_f_df(und->today, start_date, und->dom_yc) * pd_not;
  }
  if (partial_pv) {
    partial_pd_leg -= swp_f_df(und->today, cpd->call[0].set_date, und->dom_yc) *
                      cpd->call[0].pd_not_amt;
  }

  /*	PV of coupons fixed in the past and not yet paid */
  i = 0;
  while (i < pd_ncpn && pd_fix[i] < und->today + eod_fix_flag) {
    if (pd_pay[i] >= und->today + eod_pay_flag &&
        (!exercised || pd_start[i] < ex_date_ex)) {
      err = interp_basis(pd_basis[i], &bas);
      if (err)
        goto FREE_RETURN;

      if (use_cpn_opt_str) {
        opt_string = 0.0;

        for (str_idx = 0; str_idx < pd_num_strikes[i]; str_idx++) {
          if (fabs(pd_weights[i][str_idx]) > 1e-16 &&
              pd_fix_fx[i] - pd_strikes[i][str_idx] > 0.0)
            opt_string += pd_weights[i][str_idx] *
                          (pd_fix_fx[i] - pd_strikes[i][str_idx]);
        }

        pd_leg += (pd_wcst[i] + pd_wspot[i] * pd_fix_fx[i] + opt_string) *
                  coverage(pd_start[i], pd_pay[i], bas) * pd_not *
                  swp_f_df(und->today, pd_pay[i], und->dom_yc);
      } else {
        fixing = pd_alpha[i] + pd_beta[i] * pd_fix_fx[i];

        if (pd_floored[i] && fabs(pd_beta[i]) > 1.0e-16 &&
            pd_floor[i] > -1.0e-16) {
          fixing = max(fixing, pd_floor[i]);
        } else {
          if (pd_floor[i] <= -1.0e-16) {
            err = "Coupon floor must be positive";
            goto FREE_RETURN;
          }
        }

        if (pd_capped[i] && fabs(pd_beta[i]) > 1.0e-16 &&
            pd_cap[i] > -1.0e-16) {
          fixing = min(fixing, pd_cap[i]);
        } else {
          if (pd_cap[i] <= -1.0e-16) {
            err = "Coupon cap must be positive";
            goto FREE_RETURN;
          }
        }

        pd_leg += fixing * coverage(pd_start[i], pd_pay[i], bas) * pd_not *
                  swp_f_df(und->today, pd_pay[i], und->dom_yc);
      }
    }
    i++;
  }

  if (cpd_dlm_params->use_beta_dlm) {
    dlm_pd_leg += pd_leg - saved_pd_leg;
    dlm_partial_pd_leg += partial_pd_leg - saved_partial_pd_leg;
  }

  if (call_feat > 0 &&
      use_sabr >= 1) /* First update fees in all preinitialized structures */
  {
    /* Update cpd->call */
    for (i = 0; i < cpd->num_calls; i++) {
      cpd->call[i].orig_fee = cpd->call[i].fee;
      cpd->call[i].extra_fee = 0.0;

      cum_adj_fee = 0.0;
      for (j = cpd->call[i].pd_idx; j < cpd->pd_leg->num_cpn; j++)
        cum_adj_fee += adj_fees[j];
      cum_adj_fee += adj_fee_not;

      cum_adj_fee /= swp_f_df(und->today, cpd->call[i].set_date, und->dom_yc);

      cpd->call[i].extra_fee =
          (cpd->call[i].pay_rec ? cum_adj_fee : -cum_adj_fee);

      if (use_sabr == 2) {
        cpd->call[i].fee = cpd->call[i].orig_fee + cpd->call[i].extra_fee;
      }
    }

    /* Save the Fees */
    und->nb_fees = cpd->num_calls;
    if (!und->fees)
      und->fees = calloc(und->nb_fees, sizeof(double));
    if (!und->fees_dates)
      und->fees_dates = calloc(und->nb_fees, sizeof(double));

    if (!und->fees || !und->fees_dates) {
      err = "Memory allocation faillure in cpd_caller";
      goto FREE_RETURN;
    }

    for (i = 0; i < cpd->num_calls; i++) {
      und->fees_dates[i] = cpd->call[i].set_date;
      und->fees[i] = cpd->call[i].extra_fee;
    }
  }

  if (export_ts)
    cpd_copy_und(und, und_exp);

  /***********************/
  /* For the alpha Fudge */
  /***********************/
  free_struct = 1;

  /*	3)	If there is at least one call after today        , value call
   * feature
   */
  if (cpd_dlm_params->use_beta_dlm && fabs(und->model->dAlpha) > TINY) {
    use_calib = 0;
    skip_fill = 1;

    for (UpOrDown = 0; UpOrDown < 2; UpOrDown++) {
      err = cpd_AlphaFudgeUpdateUnd(und, UpOrDown); /* Up = 0 and Down = 1 */

      if (err)
        goto FREE_RETURN;

      /*	Free and then re-Initialise structures */
      if (UpOrDown)
        cpd_free_tree_arg(tree_arg);
      if (UpOrDown)
        cpd_free_mc_arg(mc_arg);

      err = cpd_fill_check_all_struct(
          today, use_calib, fx_spot, fx_spot_date, dom_calib, dom_und, dom_yc,
          dom_vc, dom_ref, dom_swap_freq, dom_swap_basis, dom_lam, for_calib,
          for_und, for_yc, for_vc, for_ref, for_swap_freq, for_swap_basis,
          for_lam, min_fact, max_fact, use_jumps, corr_times, correl_dom_for,
          correl_dom_fx, correl_for_fx, corr_n_times, cpd_dlm_params,
          get_ir_cash_vol, fx_mkt_vol_date, fx_mkt_vol, num_fx_mkt_vol, fx3dund,
          fund_not, fund_ccy, fund_ncpn, fund_fix, fund_start, fund_pay,
          fund_basis, fund_spr, fund_mrg, pd_not, pd_ncpn, pd_fix, pd_start,
          pd_pay, pd_basis, pd_alpha, pd_beta, pd_floored, pd_floor, pd_capped,
          pd_cap, 0, NULL, NULL, NULL, NULL, NULL, pd_not_ref_fix,
          pd_not_ref_alpha, pd_not_ref_beta, pd_not_ref_floored,
          pd_not_ref_floor, pd_not_ref_capped, pd_not_ref_cap, 0, NULL, NULL,
          NULL, NULL, NULL, call_type, ncall, pay_rec, ex_date, set_date,
          barrier, bar_type, fees, TARN_Do, req_stp, req_pth, do_pecs,
          forcetree, do_optim, force_optim, fx_bound, use_bound, bar_smooth,
          eod_fix_flag, eod_ex_flag, cpd, und, &call_feat, tree_arg, mc_arg,
          dom_vol_shift, for_vol_shift, fx_vol_shift, skip_fill);

      if (err)
        goto FREE_RETURN;

      if (call_feat > 0 && use_sabr == 2) {
        if ((call_feat == 1 && cpd->type == -1 && do_optim) ||
            call_feat == 2) /* MC or tree */
        {
          /* Update mc_arg */
          for (i = 0; i < mc_arg->nb_dates; i++) {
            cpd_prm = (CPD_PAY_ARG)mc_arg->void_prm[i];
            if (cpd_prm) {
              cpd_prm->eval_const.extra_fee =
                  cpd->call[cpd_prm->call_idx].extra_fee;
              if (cpd->call[cpd_prm->call_idx].pay_rec)
                cpd_prm->eval_const.extra_fee = -cpd_prm->eval_const.extra_fee;
              cpd_prm->eval_const.orig_fee = cpd_prm->eval_const.fee;
              cpd_prm->eval_const.fee =
                  cpd_prm->eval_const.orig_fee + cpd_prm->eval_const.extra_fee;
            }
          }
        } else {
          /* Update tree_arg */
          for (i = 0; i < tree_arg->nstp; i++) {
            cpd_prm = (CPD_PAY_ARG)tree_arg->void_prm[i];
            if (cpd_prm) {
              cpd_prm->eval_const.extra_fee =
                  cpd->call[cpd_prm->call_idx].extra_fee;
              if (cpd->call[cpd_prm->call_idx].pay_rec)
                cpd_prm->eval_const.extra_fee = -cpd_prm->eval_const.extra_fee;
              cpd_prm->eval_const.orig_fee = cpd_prm->eval_const.fee;
              cpd_prm->eval_const.fee =
                  cpd_prm->eval_const.orig_fee + cpd_prm->eval_const.extra_fee;
            }
          }
        }
      }

      if (call_feat == 1) {
        switch (cpd->type) {
        case -1:
          if (do_optim) {
            smessage("Launching optim MC");
            err = cpd_launch_opti_mc(cpd, und, mc_arg, &call, &call_std,
                                     do_infos, optim_bar);
          } else {
            smessage("Launching tree        , time steps requested: %d        "
                     ", actual: %d",
                     req_stp, tree_arg->nstp);
            err = cpd_launch_tree_ko(cpd, und, tree_arg, &call);
          }

          break;

        case 0:
          if (cpd_dlm_params->use_tree_iv_for_call)
            temp = 2;
          else
            temp = 1;

          smessage("Launching tree        , time steps requested: %d        , "
                   "actual: %d",
                   req_stp, tree_arg->nstp);
          err = cpd_launch_tree(cpd, und, tree_arg, &call, &temp);
          break;
        }

        if (err) {
          goto FREE_RETURN;
        }

        if (!do_optim || cpd->type == 0) {
          if (cpd_dlm_params->use_beta_dlm) {
            if (!cpd_dlm_params->use_tree_iv_for_call || cpd->type == -1)
              temp = dlm_partial_pd_leg - partial_fund_leg;
          } else {
            temp = partial_pd_leg - partial_fund_leg;
          }

          if (cpd->call[0].pay_rec) {
            temp *= -1;
          }
          call += temp;
          call_std = 0.0;
        }
      } else if (call_feat == 2) {
        smessage("Launching MC");
        err = cpd_launch_mc(cpd, und, mc_arg, &call, &call_std);
        if (err) {
          goto FREE_RETURN;
        }

        if (cpd_dlm_params->use_beta_dlm) {
          temp = dlm_partial_pd_leg - partial_fund_leg;
          if (cpd->call[0].pay_rec) {
            temp *= -1;
          }
          call += temp;
        }
      } else {
        call = 0.0;
        call_std = 0.0;
      }

      *fund_val = fund_leg;

      if (use_sabr) {
        if (cpd_dlm_params->use_beta_dlm) {
          *pd_val = dlm_pd_leg;
        } else {
          *pd_val = pd_leg;
        }
      } else {
        if (cpd_dlm_params->use_beta_dlm &&
            cpd_dlm_params->use_beta_dlm_for_iv) {
          *pd_val = dlm_pd_leg;
        } else {
          *pd_val = pd_leg;
        }
      }

      *call_val = call;
      *call_stdev = call_std;

      if (!UpOrDown) {
        call_val_alpha = call;
        call_val_alpha_std = call_std;
      }
    }

    *call_val = 0.5 * (*call_val + call_val_alpha);
    *call_stdev = sqrt(0.5 * (*call_stdev * *call_stdev +
                              call_val_alpha_std * call_val_alpha_std));
  } else {
    if (call_feat > 0 && use_sabr == 2) {
      if ((call_feat == 1 && cpd->type == -1 && do_optim) ||
          call_feat == 2) /* MC or tree */
      {
        /* Update mc_arg */
        for (i = 0; i < mc_arg->nb_dates; i++) {
          cpd_prm = (CPD_PAY_ARG)mc_arg->void_prm[i];
          if (cpd_prm) {
            cpd_prm->eval_const.extra_fee =
                cpd->call[cpd_prm->call_idx].extra_fee;
            if (cpd->call[cpd_prm->call_idx].pay_rec)
              cpd_prm->eval_const.extra_fee = -cpd_prm->eval_const.extra_fee;
            cpd_prm->eval_const.orig_fee = cpd_prm->eval_const.fee;
            cpd_prm->eval_const.fee =
                cpd_prm->eval_const.orig_fee + cpd_prm->eval_const.extra_fee;
          }
        }
      } else {
        /* Update tree_arg */
        for (i = 0; i < tree_arg->nstp; i++) {
          cpd_prm = (CPD_PAY_ARG)tree_arg->void_prm[i];
          if (cpd_prm) {
            cpd_prm->eval_const.extra_fee =
                cpd->call[cpd_prm->call_idx].extra_fee;
            if (cpd->call[cpd_prm->call_idx].pay_rec)
              cpd_prm->eval_const.extra_fee = -cpd_prm->eval_const.extra_fee;
            cpd_prm->eval_const.orig_fee = cpd_prm->eval_const.fee;
            cpd_prm->eval_const.fee =
                cpd_prm->eval_const.orig_fee + cpd_prm->eval_const.extra_fee;
          }
        }
      }
    }

    if (call_feat == 1) {
      switch (cpd->type) {
      case -1:
        if (do_optim) {
          smessage("Launching optim MC");
          err = cpd_launch_opti_mc(cpd, und, mc_arg, &call, &call_std, do_infos,
                                   optim_bar);
        } else {
          // Tree_KO not yet supported with option string coupon or redemption:

          if (use_cpn_opt_str || use_not_opt_str) {
            err = serror("KO in a tree temporarily not supported with option "
                         "string coupon or redemption");
            goto FREE_RETURN;
          }

          smessage("Launching tree        , time steps requested: %d        , "
                   "actual: %d",
                   req_stp, tree_arg->nstp);
          err = cpd_launch_tree_ko(cpd, und, tree_arg, &call);
        }
        break;

      case 0:
        if (cpd_dlm_params->use_tree_iv_for_call)
          temp = 2;
        else
          temp = 1;

        smessage("Launching tree        , time steps requested: %d        , "
                 "actual: %d",
                 req_stp, tree_arg->nstp);
        err = cpd_launch_tree(cpd, und, tree_arg, &call,
                              &temp); /* Calculates the IV in the tree */
        break;
      }

      if (err) {
        goto FREE_RETURN;
      }

      if (!do_optim || cpd->type == 0) {

        if (cpd_dlm_params->use_beta_dlm) {
          if (!cpd_dlm_params->use_tree_iv_for_call || cpd->type == -1)
            temp = dlm_partial_pd_leg - partial_fund_leg;
        } else {
          if (!cpd_dlm_params->use_tree_iv_for_call || cpd->type == -1)
            temp = partial_pd_leg - partial_fund_leg;
        }

        if (cpd->call[0].pay_rec) {
          temp *= -1;
        }
        call += temp;
        call_std = 0.0;
      }
    } else if (call_feat == 2) {
      smessage("Launching MC");
      err = cpd_launch_mc(cpd, und, mc_arg, &call, &call_std);
      if (err) {
        goto FREE_RETURN;
      }

      if (cpd_dlm_params->use_beta_dlm) {
        temp = dlm_partial_pd_leg - partial_fund_leg;
        if (cpd->call[0].pay_rec) {
          temp *= -1;
        }
        call += temp;
      }
    } else {
      call = 0.0;
      call_std = 0.0;
    }

    *fund_val = fund_leg;

    if (use_sabr) {
      if (cpd_dlm_params->use_beta_dlm) {
        *pd_val = dlm_pd_leg;
      } else {
        *pd_val = pd_leg;
      }
    } else {
      if (cpd_dlm_params->use_beta_dlm && cpd_dlm_params->use_beta_dlm_for_iv) {
        *pd_val = dlm_pd_leg;
      } else {
        *pd_val = pd_leg;
      }
    }

    *call_val = call;
    *call_stdev = call_std;
  }

  /* ===================================
  ======================================
                  SMILE ADJUSTMENT
  ======================================
  =================================== */

  if (!exercised && use_GMA && cpd->num_calls > 0) {
    Results = dmatrix(0, Height - 1, 0, Width - 1);

    // Test if the funding change is admissible
    if (FundingSpeedUp && cpd->fund_leg->dom_for) {
      if (cpd->call[0].ex_date >
          cpd->pd_leg->cpn[cpd->call[0].pd_idx].fx_val_date) {
        FundingSpeedUp = 0;
      }
    }

    OTC_params = calloc(1, sizeof(otcpd_params));

    cpd_init_otc_params(
        OTC_params, cpd, ncall, CopNsimul, CopDoPecs, CummulNpoints, CummulNstd,
        CummulPrecision, CummulLinear, SSLstd, SSLniter, PayoffFunction,
        fwdSmileVisuNstd,
        /* otc */
        fwdVolMethod, /*	0=Sliding 3F vol; 1=Converging 3F vol; 2=Sliding
                   Cvg Sbeta; 3=Cvg Cvg Sbeta */
        smileOtc,     /*	use smile parameters for the OTC */
        smileFee,     /*	use smile parameters for the OTC Fees */
        smileModel,   /*	0=SSL; 1=BM; 2=HestonD */
        otc,          /*	0=SSL; 1=Log Mix; 2=Beta Mix */
        BMpi,
        /* Correl */
        1, 1, 1, 1, 30, NULL, NULL, NULL, NULL,
        /* Time Dimension */
        FundingSpeedUp, nStart, oneOutOfN,
        /*Fast MC*/
        FMC_do, FMC_precision, FMC_min_paths, use_GMA, 0);

    /* ==================================
    Merge the TS
    ================================== */
    err = merge_rates_fx_corr_ts(
        und->sigma_time_rates, und->sigma_dom, und->sigma_n_rates,
        und->sigma_time_rates, und->sigma_for, und->sigma_n_rates,
        und->sigma_time_fx, und->sigma_fx, und->sigma_n_fx, und->corr_times,
        und->correl_dom_for, und->correl_dom_fx, und->correl_for_fx,
        und->corr_n_times, &mergeTimes, &sigDom, &sigFor, &sigFx, &corDF,
        &corDFx, &corFFx, &mergeNtimes);

    if (err)
      goto FREE_RETURN;

    /* ==================================
    Precalculations
    ================================== */
    OTC_precalc = calloc(1, sizeof(otcpd_precalc));
    cpd_init_otc_precalc(OTC_precalc);

    /* ==================================
    Launch the pricers
    ================================== */
    if (!cpd->type) {
      err = cpd_otc_precalc(OTC_precalc, OTC_params, cpd, und, smile_mkt,
                            mergeTimes, mergeNtimes, sigDom, sigFor, sigFx,
                            corDF, corDFx, corFFx);

      if (err)
        goto FREE_RETURN;

      err = otc_pricer(cpd, und, OTC_params, OTC_precalc, smile_mkt,
                       smile_params, pd_not, mergeNtimes, mergeTimes, sigDom,
                       sigFor, sigFx, corDF, corDFx, corFFx,
                       start_date, /*	Date at which initial notional exchange
                                occurs */
                       NULL,
                       erasing_call_done, //	Don't pass to GMA all the OTC
                       erased_call_list,  //	List of erased calls
                       /*	Results */
                       Results);

      if (err)
        goto FREE_RETURN;

      if (OTC_params->use_GMA == 1 || OTC_params->use_GMA == 3)
        *smile_adjustment = Results[0][0] - Results[0][2];
      else
        *smile_adjustment = Results[0][1] - Results[0][3];

      //=====================================
      // Fill the OTC/KO info
      //=====================================
      if (GMA_Results) {
        if (OTC_params->use_GMA == 1 || OTC_params->use_GMA == 3) {
          // 3f price
          GMA_Results[0][0] = Results[0][2];
          // smile price
          GMA_Results[0][1] = Results[0][0];
        } else {
          // 3f price
          GMA_Results[0][0] = Results[0][3];
          // smile price
          GMA_Results[0][1] = Results[0][1];
        }

        for (i = 0; i < cpd->num_calls; i++) {
          // time
          GMA_Results[1 + i][0] = Results[1 + i][0];

          // Forward IV Smile
          GMA_Results[1 + i][1] = Results[1 + i][1];
          // Long OTC Smile
          GMA_Results[1 + i][2] = Results[1 + i][2];
          // Short OTC Smile
          GMA_Results[1 + i][3] = Results[1 + i][3];
          // Partial OTC Smile
          GMA_Results[1 + i][4] = Results[1 + i][4];

          // Forward IV 3F
          GMA_Results[1 + i][5] = Results[1 + i][5];
          // Long OTC 3F
          GMA_Results[1 + i][6] = Results[1 + i][6];
          // Short OTC 3F
          GMA_Results[1 + i][7] = Results[1 + i][7];
          // Partial OTC 3F
          GMA_Results[1 + i][8] = Results[1 + i][8];
        }
      }
    } else {
      if (OTC_params->KO_do_optim && use_3f_optim_barrier &&
          (cpd->type == -1 && do_optim)) {
        OTC_params->KO_do_optim = 0;
        for (i = 0; i < cpd->num_calls; i++) {
          if (1 - cpd->call[i].cxxall_type) {
            cpd->call[i].orig_barrier = optim_bar[0][i + 1][0];
          }
        }
      }

      // first pricing with 3F in order to save the barrier in the case where
      // the CALLABLE/KO has been priced in the tree and we want to be
      // conservative on the smile adj. (use the same barrier)
      Smile_model = OTC_params->smileModel;
      Smile_OTC = OTC_params->OTCsmile;
      Smile_Fee_OTC = OTC_params->OTCsmileFee;
      Smile_fwd_vol = OTC_params->FwdVolMethod;

      OTC_params->smileModel = 0;
      OTC_params->OTCsmile = 0;
      OTC_params->OTCsmileFee = 0;
      OTC_params->FwdVolMethod = 4;

      err = cpd_otc_precalc(OTC_precalc, OTC_params, cpd, und, smile_mkt,
                            mergeTimes, mergeNtimes, sigDom, sigFor, sigFx,
                            corDF, corDFx, corFFx);

      if (err)
        goto FREE_RETURN;

      err = ko_pricer(cpd, und, OTC_params, OTC_precalc, smile_mkt,
                      smile_params, pd_not, mergeNtimes, mergeTimes, sigDom,
                      sigFor, sigFx, corDF, corDFx, corFFx,
                      start_date, /*	Date at which initial notional exchange
                               occurs */
                      Results);

      *smile_adjustment = -(Results[0][0] - Results[1][0]);

      //=====================================
      // Fill the OTC/KO info
      //=====================================
      // 3F price
      if (GMA_Results) {
        GMA_Results[0][0] = Results[0][0] - Results[1][0];

        for (i = 0; i < cpd->num_calls; i++) {
          // date
          GMA_Results[1 + i][0] = Results[2 + i][0];

          // 3F Proba of paying a coupon
          GMA_Results[1 + i][5] = Results[2 + i][1];
          // 3F KO price
          GMA_Results[1 + i][6] = Results[2 + i][2];
          // 3F KO IV
          GMA_Results[1 + i][7] = Results[2 + i][5] - Results[2 + i][6];
          // 3F Optimal boundary
          if (OTC_params->KO_do_optim)
            GMA_Results[1 + i][8] = Results[1 + i][7];
        }
      }

      if (OTC_params)
        cpd_free_otc_params(OTC_params);
      if (OTC_precalc)
        cpd_free_otc_precalc(OTC_precalc, cpd->num_calls, cpd->pd_leg->num_cpn,
                             cpd->fund_leg->num_cpn);
      Results[0][0] = 0.0;
      Results[1][0] = 0.0;

      // second pricing with smile
      if (OTC_params->KO_do_optim && use_3f_optim_barrier &&
          (cpd->type == -1 && !do_optim))
        OTC_params->KO_do_optim =
            0; // the barrier has allready been updated in the cpd args

      OTC_params->smileModel = Smile_model;
      OTC_params->OTCsmile = Smile_OTC;
      OTC_params->OTCsmileFee = Smile_Fee_OTC;
      OTC_params->FwdVolMethod = Smile_fwd_vol;

      err = cpd_otc_precalc(OTC_precalc, OTC_params, cpd, und, smile_mkt,
                            mergeTimes, mergeNtimes, sigDom, sigFor, sigFx,
                            corDF, corDFx, corFFx);

      if (err)
        goto FREE_RETURN;

      err = ko_pricer(cpd, und, OTC_params, OTC_precalc, smile_mkt,
                      smile_params, pd_not, mergeNtimes, mergeTimes, sigDom,
                      sigFor, sigFx, corDF, corDFx, corFFx,
                      start_date, /*	Date at which initial notional exchange
                               occurs */
                      Results);

      *smile_adjustment += Results[0][0] - Results[1][0];

      //=====================================
      // Fill the OTC/KO info
      //=====================================
      // Smile price
      if (GMA_Results) {
        GMA_Results[0][1] = Results[0][0] - Results[1][0];

        for (i = 0; i < cpd->num_calls; i++) {
          // Smile Proba of paying a coupon
          GMA_Results[1 + i][1] = Results[2 + i][1];
          // Smile KO price
          GMA_Results[1 + i][2] = Results[2 + i][2];
          // Smile KO IV
          GMA_Results[1 + i][3] = Results[2 + i][5] - Results[2 + i][6];
          // Smile Optimal boundary
          if (OTC_params->KO_do_optim)
            GMA_Results[1 + i][4] = Results[1 + i][7];
        }
      }
    }
  } else {
    *smile_adjustment = 0.0;
  }

FREE_RETURN:

  if (smile_mkt) {
    cpd_free_smile_vol_market(smile_mkt);
    free(smile_mkt);
  }

  if (smile_params)
    free(smile_params);

  free(pd_num_strikes);
  free(pd_wcst);
  free(pd_wspot);
  if (pd_strikes)
    free_dmatrix(pd_strikes, 0, alloc_pd_ncpn - 1, 0, max_interp_pts - 1);
  if (pd_weights)
    free_dmatrix(pd_weights, 0, alloc_pd_ncpn - 1, 0, max_interp_pts - 1);
  free(pd_not_num_strikes);
  free(pd_not_wcst);
  free(pd_not_wspot);
  if (pd_not_strikes)
    free_dmatrix(pd_not_strikes, 0, alloc_ncall, 0, max_interp_pts - 1);
  if (pd_not_weights)
    free_dmatrix(pd_not_weights, 0, alloc_ncall, 0, max_interp_pts - 1);

  if (OTC_params)
    cpd_free_otc_params(OTC_params);
  if (OTC_precalc)
    cpd_free_otc_precalc(OTC_precalc, cpd->num_calls, cpd->pd_leg->num_cpn,
                         cpd->fund_leg->num_cpn);

  if (free_struct)
    cpd_free_all_struct(cpd, und, tree_arg, mc_arg);
  free(new_call_type);
  free(new_ex_date);
  free(new_set_date);
  free(new_barrier);
  free(new_bar_type);
  free(new_fees);
  free(adj_fees);

  free(new_pd_not_ref_nfxpts);
  if (new_pd_not_ref_fxpts)
    free_dmatrix(new_pd_not_ref_fxpts, 0, alloc_ncall, 0, max_interp_pts - 1);
  if (new_pd_not_ref_cpn_at_pts)
    free_dmatrix(new_pd_not_ref_cpn_at_pts, 0, alloc_ncall, 0,
                 max_interp_pts - 1);
  free(new_pd_not_ref_lin_xtrpl_l);
  free(new_pd_not_ref_lin_xtrpl_r);
  free(new_pd_not_ref_fix);
  free(new_pd_not_ref_fix_fx);

  if (Inst) {
    FxBetaDLM_Free_FxOptInst(Inst);
    free(Inst);
  }

  if (InstConst)
    free(InstConst);

  if (Results)
    free_dmatrix(Results, 0, Height - 1, 0, Width - 1);
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
