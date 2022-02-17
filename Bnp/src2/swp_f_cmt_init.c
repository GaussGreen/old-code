/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT  , Fixed Income 2020 Addins */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SWP_F_CMT_INIT	                                      */
/*                                                                            */
/*      PURPOSE:        Initilise all that is necessay from the CMS-CMT       */
/*			spread curve                                          */
/*                                                                            */
/*      AUTHORS:        Olivier VAN EYSEREN                    		      */
/*                                                                            */
/*      DATE:           24th October 1995                                     */
/*                                                                            */
/*      VERSION:        01                                                    */
/*                                                                            */
/*      DESCRIPTION:    Fully initialises a CMT_Obj (that has already been
                        allocated) from the CMT-CMS spread curve (entered
                        with a list of tenors and spreads in bp)
                        The CMT_Obj will contain all the informations
                        necessary to rebuild each instrument used to strip
                        the spread curve:
                                start_date
                                end_date
                                compounding
                                basis
                                first_full_fixing
                        (for the CMT and the CMS leg of the swap)
                                                                              */
/*                                                                            */
/*      FUNCTIONS USED: XXX_X_XXXXXXXXX                                       */
/*              Must include all imported function call made by the module    */
/*                                                                            */
/*      PARAMETERS:     <not applicable>                                      */
/*                                                                            */
/*      RETURNS:                                                              */
/*                                                                            */
/*      DATA ACCESSED:                                                        */
/*      	<none>                                                        */
/*                                                                            */
/******************************************************************************/

#include "math.h"
#include "num_h_allhdr.h"
#include "swp_h_all.h"

static Err CMT_set_Instr(CMT_Obj *cmt, CMTInstr *instr, int i);

static Err string_to_CMTInstr(String iname, CMTInstr *instr, double *spread,
                              CMT_Param *cmt_param, Date spot);

static Err CMT_init(CMT_Obj *cmt, Date cmt_today_date,
                    CMT_Param_Struct *cmt_param, String *tenor_names,
                    double *spread_rates, int num_spreads);

static Err FwdT_init(String fwd_T_name, CMT_Obj *cmt, FwdT_Obj *fwdT,
                     CMT_Param *cmt_param);

static Err FwdTEC_init(String fwd_T_name, CMT_Obj *cmt, FwdT_Obj *fwdT,
                       CMT_Param *cmt_param);

static int putoncurve(FwdT_Obj *fwdT, Ddate date, double cum_time,
                      double fwd_spread, double fwd_R_rate, long index);

static Err FwdT_fill(String fwd_spread_name, double *fwd_spread_date,
                     double *fwd_spread_value, int num_fwd_spreads,
                     FwdT_Obj *fwdT, SrtCurvePtr yc_crv, Date clcn_date,
                     Date spot_date, CMT_Param *cmt_param);

static SrtErr swp_f_addcmtcrv_to_list(SrtCurveList *curve_list,
                                      String fwdT_name, String curve_lbl,
                                      Date clcn_date, String yc_name,
                                      FwdT_Obj *fwdT, CMT_Obj *cmt,
                                      CMT_Param *cmtprm);

#define MAXITER 40
#define PRECISION 1e-7

static Err CMT_init(CMT_Obj *cmt, Date cmt_today_date,
                    CMT_Param_Struct *cmt_param, String *tenor_names,
                    double *spread_rates, int num_spreads) {
  int toalloc;
  int ncash = 0, nswap = 0;
  CMTInstr cmt_instr;
  SrtCurvePtr crv;
  Err err;
  int i;
  Date today, spot_date;

  /*  Check cmt_today_date is bus day */
  if (week_day(cmt_today_date) == SAT || week_day(cmt_today_date) == SUN)
    return serror("TODAY date should be a weekday.");

  /*  Check cmt_today_date vs. today (given by reference yield curve */
  crv = lookup_curve(cmt_param->yc_name);
  if (!crv)
    return serror("Unknown yield curve %s in CMT_init", cmt_param->yc_name);

  today = get_clcndate_from_yldcrv(crv);
  spot_date = get_spotdate_from_yldcrv(crv);
  if (today != cmt_today_date)
    return serror("Today dates do not match (yc != cmt)");

  /* Attach TODAY and SPOT to the CMT_Obj */
  CMT_Dset(cmt, CMT_SPOT, spot_date);
  CMT_Dset(cmt, CMT_TODAY, today);

  /* Attach the CMTCode to the CMT_Obj */
  CMT_Iset(cmt, CMT_CODE, cmt_param->cmt_code);

  /* Changes spreads to percentages (they were entered as bp)	*/
  for (i = 0; i < num_spreads; i++)
    spread_rates[i] *= .0001;

  /* Create space in cmt	*/
  toalloc = num_spreads + 4;
  CMT_field_set_Dlength(cmt, toalloc);
  CMT_field_set_Ilength(cmt, toalloc);

  if ((!strcmp(tenor_names[0], "SPOT")) || (!strcmp(tenor_names[0], "S"))) {
    CMT_Dset(cmt, CMT_SPOT_SPREAD, spread_rates[0]);
    CMT_Dset(cmt, CMT_SPOT_DATE, spot_date);
    tenor_names++;
    spread_rates++;
    num_spreads--;

    /* Flag set to check if there is a spot declared */
    CMT_Iset(cmt, CMT_SPOT_EXISTS, 1);
  }

  /* Interpret spread tenor strings as dates */
  for (i = 0; i < num_spreads; i++) {
    /* Computes the end dates (a full CMTInstr) of the swap  , given the tenor
     */
    if (err = string_to_CMTInstr(tenor_names[i], &cmt_instr, &spread_rates[i],
                                 cmt_param, spot_date)) {
      return err;
    }

    /* Attach the corresponding instruments to the CMT_Obj*/
    switch (cmt_instr.type) {
    case CASHINSTR:
      if (nswap > 0)
        return serror("Cash rates must come before swap rates.");

      if (err = CMT_set_Instr(cmt, &cmt_instr, ncash))
        return err;

      ncash++;
      break;

    case SWAPINSTR:
      if (err = CMT_set_Instr(cmt, &cmt_instr, nswap))
        return err;

      nswap++;
      break;

    default:
      return serror("CMT_INIT: unknown instrument encountered");
    }
  }

  CMT_Iset(cmt, CMT_NUMSWAP, nswap);
  CMT_Iset(cmt, CMT_NUMCASH, ncash);

  return NULL;
}

/* ======================================================================== */

static Err string_to_CMTInstr(String iname, CMTInstr *instr, double *spread,
                              CMT_Param *cmt_param, Date spot) {
  int ny = 0, nm = 0, nd = 0;
  int num_brk_month = 0;
  Err err;

  /* Reads the tenor string and returns the number of days  , months  , years to
     be added */

  if (err = interp_tenor_string(iname, &ny, &nm, &nd))
    return err;

  instr->cms_dp.cashdates.start = spot;
  instr->cmt_dp.cashdates.start = spot;

  /* Checks if this is a cash */
  if (nd > 0) {
    return serror("Not implemented yet");
  } else if (ny == 0 && nm > 0 && nm <= 12) {
    return serror("Not implemented yet");
  } else {
    nm = nm + 12 * ny;

    /* For the CMS leg of a swap  */
    instr->cms_dp.swapdates.basis_code = cmt_param->cms_basis_code;
    instr->cms_dp.swapdates.compd = cmt_param->cms_freq;
    num_brk_month = nm % (12 / (cmt_param->cms_freq));
    instr->cms_dp.swapdates.end =
        add_unit(spot, nm, SRT_MONTH, NO_BUSDAY_CONVENTION);
    instr->cms_dp.swapdates.nfp = nm * (cmt_param->cms_freq) / 12;
    if (num_brk_month)
      instr->cms_dp.swapdates.first_full_fixing =
          add_unit(spot, num_brk_month, SRT_MONTH, cmt_param->cms_bus_day_conv);
    else
      instr->cms_dp.swapdates.first_full_fixing = spot;
    instr->type = SWAPINSTR;
    instr->cms_dp.swapdates.spot_lag = cmt_param->spot_lag;
    instr->cms_dp.swapdates.direction = BKWD;

    /* For the CMT leg of a swap */
    instr->cmt_dp.swapdates.basis_code = cmt_param->cmt_basis_code;
    instr->cmt_dp.swapdates.compd = cmt_param->cmt_freq;
    num_brk_month = nm % (12 / (cmt_param->cmt_freq));
    instr->cmt_dp.swapdates.end =
        add_unit(spot, nm, SRT_MONTH, NO_BUSDAY_CONVENTION);
    instr->cmt_dp.swapdates.nfp = nm * (cmt_param->cmt_freq) / 12;
    if (num_brk_month)
      instr->cmt_dp.swapdates.first_full_fixing =
          add_unit(spot, num_brk_month, SRT_MONTH, cmt_param->cmt_bus_day_conv);
    else
      instr->cmt_dp.swapdates.first_full_fixing = spot;
    instr->type = SWAPINSTR;
    instr->cmt_dp.swapdates.spot_lag = cmt_param->spot_lag;
    instr->cmt_dp.swapdates.direction = BKWD;
  }

  if (spread)
    instr->spread = *spread;

  return NULL;
}

/* ======================================================================== */

static Err CMT_set_Instr(CMT_Obj *cmt, CMTInstr *instr, int i) {
  int max_index;

  if (i < 0)
    return serror("Bad instance index %d", i);

  switch (instr->type) {
  case SWAPINSTR:
    max_index = CMT_field_Dlength(cmt, CMS_CMT_SPREAD_RATE);
    if (i >= max_index)
      return serror("Bad spread instance index %d", i);

    CMT_field_Dset(cmt, CMS_CMT_SPREAD_RATE, i, instr->spread);

    /* Default values */
    CMT_field_Dset(cmt, CMT_SPREAD_LIBOR, i, 0);

    CMT_field_Iset(cmt, CMT_SWAP_BASIS_CODE, i,
                   instr->cmt_dp.swapdates.basis_code);
    CMT_field_Iset(cmt, CMS_SWAP_BASIS_CODE, i,
                   instr->cms_dp.swapdates.basis_code);

    CMT_field_Dateset(cmt, CMT_SWAPEND_DATE, i, instr->cmt_dp.swapdates.end);
    CMT_field_Dateset(cmt, CMS_SWAPEND_DATE, i, instr->cms_dp.swapdates.end);

    CMT_field_Dateset(cmt, CMT_SWAP_FIRST_FULL_FIXING, i,
                      instr->cmt_dp.swapdates.first_full_fixing);
    CMT_field_Dateset(cmt, CMS_SWAP_FIRST_FULL_FIXING, i,
                      instr->cms_dp.swapdates.first_full_fixing);

    CMT_field_Iset(cmt, CMT_SWAP_NUM_FULL_PERIOD, i,
                   instr->cmt_dp.swapdates.nfp);
    CMT_field_Iset(cmt, CMS_SWAP_NUM_FULL_PERIOD, i,
                   instr->cms_dp.swapdates.nfp);

    CMT_field_Iset(cmt, CMT_SWAP_FREQ, i, instr->cmt_dp.swapdates.compd);
    CMT_field_Iset(cmt, CMS_SWAP_FREQ, i, instr->cms_dp.swapdates.compd);

    if (i > 0 && instr->cmt_dp.swapdates.end <=
                     CMT_field_Dateget(cmt, CMT_SWAPEND_DATE, i - 1))
      return serror("CMT_set_instr swaps must be in order.");
    break;

  case CASHINSTR:
    return serror("Cash rates not implemented yet");
    /*
    See swp_f_ycinstr.c : Err YC_set_Instr
    */

    break;

  default:
    return serror("Unknown Instrument %d", instr->type);
  }

  return NULL;
}

/* ======================================================================== */

/* fwd_spread is negative */
static int putoncurve(FwdT_Obj *fwdT, Ddate date, double cum_time,
                      double fwd_spread, double fwd_R_rate, long index) {
  FwdT_field_Dset(fwdT, FwdT_DATE, index, date);
  FwdT_field_Dset(fwdT, FwdT_TIME, index, cum_time);
  FwdT_field_Dset(fwdT, FwdT_SPREAD, index, fwd_spread);
  FwdT_field_Dset(fwdT, FwdT_SWAP_RATE, index, fwd_R_rate);
  FwdT_field_Dset(fwdT, FwdT_RATE, index, fwd_R_rate + fwd_spread);
  FwdT_Iset(fwdT, FwdT_NUMSPREAD, (index + 1));

  return (0);
}

/* ======================================================================== */

static Err FwdT_init(String fwd_T_name, CMT_Obj *cmt, FwdT_Obj *fwdT,
                     CMT_Param *cmt_param)

{
  Err solve_flag = NULL;

  Arg_Obj *cms_arg, *cmt_arg;
  Swap_Obj *cms_swap, *cmt_swap;
  Leg_Obj *cms_leg, *cmt_leg;
  int toalloc;
  int num_full_period;
  int first_swap_index = 0;
  Date swap_end;
  double first_full_fixing;
  double fwd_spread, roll_spread_lock, cms_libor_spread, cmt_libor_spread;
  double cum_time;
  double cms_swap_price, cmt_swap_price, cms_cmt_price;
  Date spot, today, tenor_date, last_fixing, start_date;
  double fwd_T_rate, fwd_R_rate, a[3], b[3], yans;
  int i, index, niter;
  int count, num_swap, swap_len;
  double nstop;
  CMTCode cmt_code;
  SwapDP swapdp;
  DateList datelist;
  SrtCurvePtr
      crv; /* Can be either a cmt_crv or a yc_crv depending on last use*/
  String yc_name;
  String ccy_str;

  crv = lookup_curve(fwd_T_name);
  yc_name = get_ycname_from_cmtcrv(crv);

  /* Get the Ccy of the Market */
  ccy_str = get_curve_ccy(crv);

  /* Stores dates used everywhere in the code */
  spot = (Date)CMT_Dget(cmt, CMT_SPOT);
  today = (Date)CMT_Dget(cmt, CMT_TODAY);
  cmt_code = CMT_Dget(cmt, CMT_CODE);

  /* Number of swaps entered in the curve (with spot if entered)*/
  num_swap = CMT_Iget(cmt, CMT_NUMSWAP);
  index = 0;

  /* Memory allocation for the FwdT_Obj (use  more than necesary)*/
  swap_len = CMT_field_Dlength(cmt, CMS_CMT_SPREAD_RATE);
  toalloc = 2 * swap_len + 2 * (CMT_Iget(cmt, CMT_SPOT_EXISTS));
  FwdT_field_set_Dlength(fwdT, toalloc);
  FwdT_field_set_Ilength(fwdT, toalloc);
  FwdT_Dset(fwdT, FwdT_TODAY, today);
  FwdT_Dset(fwdT, FwdT_SPOT, spot);
  FwdT_Dset(fwdT, FwdT_CMTCODE, cmt_code);

  if (CMT_Iget(cmt, CMT_SPOT_EXISTS == 1)) {
    fwd_spread = -CMT_Dget(cmt, CMT_SPOT_SPREAD);
    start_date = (Date)CMT_Dget(cmt, CMT_SPOT_DATE);
    cum_time = (start_date - spot) * YEARS_IN_DAY;
    tenor_date = bus_date_method(start_date, MODIFIED_SUCCEEDING);

    /* Computes the (forward) swap rate (start date = spot_date of CMT swap) */
    swapdp.start = start_date;
    swapdp.nfp = (int)cmt_param->cmt_mat * cmt_param->swap_compd;
    swapdp.end = swapdp.nfp;
    swapdp.direction = FWD;
    swapdp.compd = cmt_param->swap_compd;
    swapdp.basis_code = cmt_param->swap_basis_code;
    swapdp.spot_lag = cmt_param->spot_lag;
    swapdp.first_full_fixing = swapdp.start;
    swp_f_ForwardRate_SwapDP(&swapdp, yc_name, cmt_param->swap_ref_rate,
                             &fwd_R_rate);

    /* Sets the initial value for the fwd_T_rate in the curve:
            please note that FwdT_NUMSPREAD is set to  index+1 */
    putoncurve(fwdT, tenor_date, cum_time, fwd_spread, fwd_R_rate, index);

    index++;
    /*
                    num_swap --;
    */
    /* As there is a spot rate entered
                                    and first_swap_index refers to first swap
       that does not overlap*/
  }

  /* Index is used to store the number of spreads already in the curve */
  index = FwdT_Iget(fwdT, FwdT_NUMSPREAD);

  /* Initialises in the CMS ARG what is independent of the instrument */
  /*
  Err swp_f_create_GenFullSwap(
                          *rec_leg_sdp  ,
                          *pay_leg_sdp  ,
                          fixed  ,
                          initial_exchange  ,
                          double       final_exchange  ,
                          SwapType     swap_type  ,
                          SrtCurvePtr    crv  ,
                          GenFullSwap  *swap)
  */
  init_Arg(&cms_arg, GENERATE_SWAP);
  Arg_Dset(cms_arg, ARG_NOTIONAL, 1000000000.00);
  Arg_Dset(cms_arg, ARG_INITIAL_NOT, 1);
  Arg_Dset(cms_arg, ARG_FINAL_NOT, 1);
  Arg_Dset(cms_arg, ARG_STRIKE, 0.0);
  Arg_set(cms_arg, ARG_COMPD, cmt_param->cms_freq);
  Arg_set(cms_arg, ARG_BASIS_CODE, cmt_param->cms_basis_code);
  Arg_set(cms_arg, ARG_DATE_DIR, BKWD);
  Arg_Dateset(cms_arg, ARG_START, spot);
  Arg_Dateset(cms_arg, ARG_TODAY, today);
  Arg_Dateset(cms_arg, ARG_SPOT, spot);

  /* Initialises in the CMT ARG what is independent of the instrument */

  init_Arg(&cmt_arg, GENERATE_SWAP);
  Arg_Dset(cmt_arg, ARG_NOTIONAL, 1000000000.00);
  Arg_Dset(cmt_arg, ARG_INITIAL_NOT, 1);
  Arg_Dset(cmt_arg, ARG_FINAL_NOT, 1);
  Arg_Dset(cmt_arg, ARG_STRIKE, 0.0);
  Arg_set(cmt_arg, ARG_COMPD, cmt_param->cmt_freq);
  Arg_set(cmt_arg, ARG_BASIS_CODE, cmt_param->cmt_basis_code);
  Arg_set(cmt_arg, ARG_DATE_DIR, BKWD);
  Arg_Dateset(cmt_arg, ARG_START, spot);
  Arg_Dateset(cmt_arg, ARG_TODAY, today);
  Arg_Dateset(cmt_arg, ARG_SPOT, spot);

  if (num_swap > 0) {
    for (i = first_swap_index; i < num_swap; i++) {
      /* Gets the value of the CMS-CMT spread from the spreadsheet (thru
       * CMT_Obj)*/
      roll_spread_lock = CMT_field_Dget(cmt, CMS_CMT_SPREAD_RATE, i);

      /* THIS CONCERNS THE CMS-SPREAD / LIBOR SWAP */
      /* Sets in the CMS ARG what is instrument dependent*/
      first_full_fixing = CMT_field_Dget(cmt, CMS_SWAP_FIRST_FULL_FIXING, i);
      num_full_period = CMT_field_Iget(cmt, CMS_SWAP_NUM_FULL_PERIOD, i);
      swap_end = (Date)CMT_field_Dget(cmt, CMS_SWAPEND_DATE, i);

      Arg_Iset(cms_arg, ARG_NUM_FULL_PERIOD, num_full_period);
      Arg_Dset(cms_arg, ARG_FIRST_FULL_FIXING, first_full_fixing);
      Arg_Dateset(cms_arg, ARG_END, swap_end);

      /* STARTS ITERATION ON SPREAD TO SOLVE FOR A PAR VALUE OF THE CMS/LIBOR
       * SWAP */

      /* Sets a few useful parameters */
      niter = 5;
      yans = 0.0;
      count = 0;
      nstop = 0.0;

      /* Build the CMS-spread / LIBOR swap */
      if (cmt_code == TEC10)
        cms_swap = generate_swap(cms_arg, CMS_FOR_TEC_MARGIN_FLOATING_SWAP);
      else
        cms_swap = generate_swap(cms_arg, CMS_MARGIN_FLOATING_SWAP);

      /* Extracts the CMS-spread Leg of the CMS - Libor swap */
      cms_leg = Swap_field_get(cms_swap, SWAP_LEG, 0);

      /* Populates the CMS-spread Leg of the CMS-Libor swap : fwds  , vols...*/
      crv = lookup_curve(fwd_T_name);
      leg_cms_populate(cms_leg, crv);

      /* Sets a first guess for the CMS/LIBOR spread */
      cms_libor_spread = 0;
      a[0] = cms_libor_spread;
      Arg_Dset(cms_arg, ARG_SPREAD, cms_libor_spread);
      /* Make sure the spread is changed inside the swap leg:
         This prevents from rebuilding the swap leg at each time */
      Leg_Dset(cms_leg, LEG_SPREAD, cms_libor_spread);

      /* Prices the CMS/LIBOR swap */
      crv = lookup_curve(fwd_T_name);
      cms_swap_price = value_swap(crv, cms_swap);
      b[0] = cms_swap_price;

      /* Does this again for two shifts */
      cms_libor_spread -= 0.0001;
      a[1] = cms_libor_spread;
      Arg_Dset(cms_arg, ARG_SPREAD, cms_libor_spread);
      cms_leg = Swap_field_get(cms_swap, SWAP_LEG, 0);
      Leg_Dset(cms_leg, LEG_SPREAD, cms_libor_spread);
      cms_swap_price = value_swap(crv, cms_swap);
      b[1] = cms_swap_price;

      cms_libor_spread -= 0.0001;
      a[2] = cms_libor_spread;

      /* Calls NEWTON for solving the CMS-LIBOR spread*/
      while ((count < MAXITER) && (nstop < 1.0)) {
        cms_libor_spread = a[2];
        Arg_Dset(cms_arg, ARG_SPREAD, cms_libor_spread);
        cms_leg = Swap_field_get(cms_swap, SWAP_LEG, 0);
        Leg_Dset(cms_leg, LEG_SPREAD, cms_libor_spread);
        cms_swap_price = value_swap(crv, cms_swap);
        b[2] = cms_swap_price;
        newton(yans, niter, a, b, &nstop);
        count++;
      }

      /* Checks convergence */
      if (nstop < 1.0)
        solve_flag = "UNABLE TO SOLVE FOR CMS-LIBOR SPREAD ";

      /* Put the last value in the cms_arg  , the cms_swap */

      cms_libor_spread = a[2];
      Arg_Dset(cms_arg, ARG_SPREAD, cms_libor_spread);
      cms_leg = Swap_field_get(cms_swap, SWAP_LEG, 0);
      Leg_Dset(cms_leg, LEG_SPREAD, cms_libor_spread);

      /* Memory free */
      free_Swap(cms_swap);

      /* THIS CONCERNS THE CMT-SPREAD / LIBOR SWAP */
      /* Sets in the CMT ARG what is instrument dependent*/
      first_full_fixing = CMT_field_Dget(cmt, CMT_SWAP_FIRST_FULL_FIXING, i);
      num_full_period = CMT_field_Iget(cmt, CMT_SWAP_NUM_FULL_PERIOD, i);
      swap_end = (Date)CMT_field_Dget(cmt, CMT_SWAPEND_DATE, i);

      Arg_Iset(cmt_arg, ARG_NUM_FULL_PERIOD, num_full_period);
      Arg_Dset(cmt_arg, ARG_FIRST_FULL_FIXING, first_full_fixing);
      Arg_Dateset(cmt_arg, ARG_END, swap_end);

      /* Now that cms_libor_spread is known  , we can store cmt_libor spread */
      cmt_libor_spread = cms_libor_spread + roll_spread_lock;
      ;
      Arg_Dset(cmt_arg, ARG_SPREAD, cmt_libor_spread);

      /* Computes the last fixing date in the swap
         (last payment is made at swap_end) */
      swapdp.end = swap_end;
      swapdp.nfp = swapdp.end;
      swapdp.direction = BKWD;
      swapdp.compd = cmt_param->cmt_freq;
      swapdp.basis_code = cmt_param->cmt_basis_code;
      swapdp.spot_lag = cmt_param->spot_lag;
      swapdp.start = today;
      swapdp.first_full_fixing = swapdp.start;

      datelist = SwapDP_to_DateList(&swapdp, MODIFIED_SUCCEEDING);
      last_fixing = datelist.date[datelist.len - 2];

      /* Computes the forward swap rate (start date = last_fixing of CMT swap)
       */
      swapdp.nfp = (int)cmt_param->cmt_mat * cmt_param->swap_compd;
      swapdp.end = swapdp.nfp;
      swapdp.direction = FWD;
      swapdp.compd = cmt_param->swap_compd;
      swapdp.basis_code = cmt_param->swap_basis_code;
      swapdp.spot_lag = cmt_param->spot_lag;
      swapdp.start = last_fixing;
      swapdp.first_full_fixing = swapdp.start;
      swp_f_ForwardRate_SwapDP(&swapdp, yc_name, cmt_param->swap_ref_rate,
                               &fwd_R_rate);

      /* Sets tenor_date to be the right date in the fwd_T_curve */
      cum_time = (last_fixing - spot) * YEARS_IN_DAY;
      tenor_date = bus_date_method(last_fixing, MODIFIED_SUCCEEDING);

      /* Sets a few useful parameters */
      niter = 5;
      yans = 0.0;
      count = 0;
      nstop = 0.0;

      /* STARTS ITERATION TO SOLVE FOR A PAR VALUE OF THE CMS_CMT SWAP */

      /* Sets an initial guess for the new fwd treasury rate */
      fwd_spread = -roll_spread_lock;
      fwd_T_rate = fwd_R_rate + fwd_spread;
      a[0] = fwd_spread;

      /* Sets the initial guess for the fwd_T_rate in the curve */
      putoncurve(fwdT, tenor_date, cum_time, fwd_spread, fwd_R_rate, index);

      /* Build the CMT-spread / LIBOR swap (the vols depend on the fwdT rate */
      if (cmt_code == TEC10)
        cmt_swap = generate_swap(cmt_arg, TEC_MARGIN_FLOATING_SWAP);
      else
        cmt_swap = generate_swap(cmt_arg, CMT_MARGIN_FLOATING_SWAP);

      /* Extracts the CMT-spread Leg of the CMS - Libor swap */
      cmt_leg = Swap_field_get(cmt_swap, SWAP_LEG, 0);

      /* Computes and store fwd swap rates  , and swap vols */
      crv = lookup_curve(fwd_T_name);
      leg_cmt_populate(cmt_leg, crv);

      /* Initial price of the CMT leg */
      cmt_swap_price = value_swap(crv, cmt_swap);

      /* Initial price of the full CMS-CMT swap */
      cms_cmt_price = cmt_swap_price - cms_swap_price;
      b[0] = cms_cmt_price;

      /* Does this again for two shifts */
      fwd_spread = fwd_spread - 0.0001;
      fwd_T_rate = fwd_R_rate + fwd_spread;
      a[1] = fwd_spread;
      putoncurve(fwdT, tenor_date, cum_time, fwd_spread, fwd_R_rate, index);
      cmt_swap_price = value_swap(crv, cmt_swap);
      cms_cmt_price = cmt_swap_price - cms_swap_price;
      b[1] = cms_cmt_price;

      fwd_spread = fwd_spread - 0.0001;
      fwd_T_rate = fwd_R_rate + fwd_spread;
      a[2] = fwd_spread;

      /* Calls NEWTON for solving */
      while ((count < MAXITER) && (nstop < 1.0)) {
        fwd_spread = a[2];
        fwd_T_rate = fwd_R_rate + fwd_spread;
        putoncurve(fwdT, tenor_date, cum_time, fwd_spread, fwd_R_rate, index);
        cmt_swap_price = value_swap(crv, cmt_swap);
        cms_cmt_price = cmt_swap_price - cms_swap_price;
        b[2] = cms_cmt_price;
        newton(yans, niter, a, b, &nstop);
        count += 1;
      }

      /* Checks convergence */
      if (nstop < 1.0)
        solve_flag = "UNABLE TO SOLVE";

      /* Put the last value in the FwdT_Obj (i.e. in the fwd_T_curve ) */
      fwd_spread = a[2];
      fwd_T_rate = fwd_R_rate + fwd_spread;
      putoncurve(fwdT, tenor_date, cum_time, fwd_spread, fwd_R_rate, index);

      /* Memory free */
      free_Swap(cmt_swap);

      /* One more spread initialised in the FwdT_Obj (i.e. in the fwd_T_curve )
       */
      index += 1;

    } /* END for (i = first_swap_index; i<num_swap; i++)  */

    free_Arg(cmt_arg);
    free_Arg(cms_arg);

  } /* END if (num_swap >0 ) */
  else {
    solve_flag = "NO SPREAD INPUT";
  }
  return (solve_flag);

} /* END FwdT_init() */

/* ------------------------------------------------------------------------ */

static Err FwdTEC_init(String fwd_T_name, CMT_Obj *tec, FwdT_Obj *fwdT,
                       CMT_Param *cmt_param)

{
  Err solve_flag = NULL;

  Arg_Obj *tec_arg;
  Swap_Obj *tec_swap;
  Leg_Obj *tec_leg;
  int toalloc;
  int num_full_period;
  int first_swap_index = 0;
  Date swap_end;
  double first_full_fixing;
  double fwd_spread, tec_libor_spread;
  double cum_time;
  double tec_swap_price;
  Date spot, today, tenor_date, last_fixing;
  double fwd_T_rate, fwd_R_rate, a[3], b[3], yans;
  int i, index, niter;
  int count, num_swap, swap_len;
  double nstop;
  CMTCode cmt_code;
  SwapDP swapdp;
  DateList datelist;
  SrtCurvePtr
      crv; /* Can be either a cmt_crv or a yc_crv depending on last use*/
  String ccy_str;
  String yc_name;

  crv = lookup_curve(fwd_T_name);
  yc_name = get_ycname_from_cmtcrv(crv);

  /* Get the Ccy */
  ccy_str = get_curve_ccy(crv);

  /* Stores dates used everywhere in the code */
  spot = (Date)CMT_Dget(tec, CMT_SPOT);
  today = (Date)CMT_Dget(tec, CMT_TODAY);
  cmt_code = CMT_Dget(tec, CMT_CODE);
  if (cmt_code != TEC10)
    return serror("Use the TEC stripper for a normal CMT");

  /* Number of swaps entered in the curve (with spot if entered)*/
  num_swap = CMT_Iget(tec, CMT_NUMSWAP);
  index = 0;

  /* Memory allocation for the FwdT_Obj (use  more than necesary)*/
  swap_len = CMT_field_Dlength(tec, CMS_CMT_SPREAD_RATE);
  toalloc = 2 * swap_len + 2 * (CMT_Iget(tec, CMT_SPOT_EXISTS));
  FwdT_field_set_Dlength(fwdT, toalloc);
  FwdT_field_set_Ilength(fwdT, toalloc);
  FwdT_Dset(fwdT, FwdT_TODAY, today);
  FwdT_Dset(fwdT, FwdT_SPOT, spot);
  FwdT_Dset(fwdT, FwdT_CMTCODE, cmt_code);

  /* Index is used to store the number of spreads already in the curve */
  index = FwdT_Iget(fwdT, FwdT_NUMSPREAD);

  /* Initialises in the TEC ARG what is independent of the instrument */

  init_Arg(&tec_arg, GENERATE_SWAP);
  Arg_Dset(tec_arg, ARG_NOTIONAL, 1000000000.00);
  Arg_Dset(tec_arg, ARG_INITIAL_NOT, 1);
  Arg_Dset(tec_arg, ARG_FINAL_NOT, 1);
  Arg_Dset(tec_arg, ARG_STRIKE, 0.0);
  Arg_set(tec_arg, ARG_COMPD, cmt_param->cmt_freq);
  Arg_set(tec_arg, ARG_BASIS_CODE, cmt_param->cmt_basis_code);
  Arg_set(tec_arg, ARG_DATE_DIR, BKWD);
  Arg_Dateset(tec_arg, ARG_START, spot);
  Arg_Dateset(tec_arg, ARG_TODAY, today);
  Arg_Dateset(tec_arg, ARG_SPOT, spot);

  if (num_swap > 0) {
    for (i = first_swap_index; i < num_swap; i++) {
      /* Gets the value of the TEC-LIBOR spread from the spreadsheet (thru
       * CMT_Obj)*/
      tec_libor_spread = CMT_field_Dget(tec, CMS_CMT_SPREAD_RATE, i);

      /* THIS CONCERNS THE CMT-SPREAD / LIBOR SWAP */
      /* Sets in the CMT ARG what is instrument dependent*/
      first_full_fixing = CMT_field_Dget(tec, CMT_SWAP_FIRST_FULL_FIXING, i);
      num_full_period = CMT_field_Iget(tec, CMT_SWAP_NUM_FULL_PERIOD, i);
      swap_end = (Date)CMT_field_Dget(tec, CMT_SWAPEND_DATE, i);

      Arg_Iset(tec_arg, ARG_NUM_FULL_PERIOD, num_full_period);
      Arg_Dset(tec_arg, ARG_FIRST_FULL_FIXING, first_full_fixing);
      Arg_Dateset(tec_arg, ARG_END, swap_end);

      /* Store the cmt_libor spread from the spreadsheet */
      Arg_Dset(tec_arg, ARG_SPREAD, tec_libor_spread);

      /* Computes the last fixing date in the swap
         (last payment is made at swap_end) */
      swapdp.end = swap_end;
      swapdp.nfp = swapdp.end;
      swapdp.direction = BKWD;
      swapdp.compd = cmt_param->cmt_freq;
      swapdp.basis_code = cmt_param->cmt_basis_code;
      swapdp.spot_lag = cmt_param->spot_lag;
      swapdp.start = today;
      swapdp.first_full_fixing = swapdp.start;

      datelist = SwapDP_to_DateList(&swapdp, MODIFIED_SUCCEEDING);
      last_fixing = datelist.date[datelist.len - 2];

      /* Computes the forward swap rate (start date = last_fixing of CMT swap)
       */
      swapdp.nfp = (int)cmt_param->cmt_mat * cmt_param->swap_compd;
      swapdp.end = swapdp.nfp;
      swapdp.direction = FWD;
      swapdp.compd = cmt_param->swap_compd;
      swapdp.basis_code = cmt_param->swap_basis_code;
      swapdp.spot_lag = cmt_param->spot_lag;
      swapdp.start = last_fixing;
      swapdp.first_full_fixing = swapdp.start;
      swp_f_ForwardRate_SwapDP(&swapdp, yc_name, cmt_param->swap_ref_rate,
                               &fwd_R_rate);

      /* Sets tenor_date to be the right date in the fwd_T_curve */
      cum_time = (last_fixing - spot) * YEARS_IN_DAY;
      tenor_date = bus_date_method(last_fixing, MODIFIED_SUCCEEDING);

      /* Sets a few useful parameters */
      niter = 5;
      yans = 0.0;
      count = 0;
      nstop = 0.0;

      /* STARTS ITERATION TO SOLVE FOR A PAR VALUE OF THE CMS_CMT SWAP */

      /* Sets an initial guess for the new fwd treasury rate */
      fwd_spread = -0.0015;
      fwd_T_rate = fwd_R_rate + fwd_spread;
      a[0] = fwd_spread;

      /* Sets the initial guess for the fwd_T_rate in the curve */
      putoncurve(fwdT, tenor_date, cum_time, fwd_spread, fwd_R_rate, index);

      /* Build the CMT-spread / LIBOR swap (the vols depend on the fwdT rate */
      tec_swap = generate_swap(tec_arg, TEC_MARGIN_FLOATING_SWAP);

      /* Extracts the CMT-spread Leg of the CMS - Libor swap */
      tec_leg = Swap_field_get(tec_swap, SWAP_LEG, 0);

      /* Computes and store fwd swap rates  , and swap vols */
      crv = lookup_curve(fwd_T_name);
      leg_cmt_populate(tec_leg, crv);

      /* Initial price of the CMT leg */
      tec_swap_price = value_swap(crv, tec_swap);

      /* Initial price of the full CMS-CMT swap */
      b[0] = tec_swap_price;

      /* Does this again for two shifts */
      fwd_spread = fwd_spread - 0.0001;
      fwd_T_rate = fwd_R_rate + fwd_spread;
      a[1] = fwd_spread;
      putoncurve(fwdT, tenor_date, cum_time, fwd_spread, fwd_R_rate, index);
      tec_swap_price = value_swap(crv, tec_swap);
      b[1] = tec_swap_price;

      fwd_spread = fwd_spread - 0.0001;
      fwd_T_rate = fwd_R_rate + fwd_spread;
      a[2] = fwd_spread;

      /* Calls NEWTON for solving */
      while ((count < MAXITER) && (nstop < 1.0)) {
        fwd_spread = a[2];
        fwd_T_rate = fwd_R_rate + fwd_spread;
        putoncurve(fwdT, tenor_date, cum_time, fwd_spread, fwd_R_rate, index);
        tec_swap_price = value_swap(crv, tec_swap);
        b[2] = tec_swap_price;
        newton(yans, niter, a, b, &nstop);
        count += 1;
      }

      /* Checks convergence */
      if (nstop < 1.0)
        solve_flag = "UNABLE TO SOLVE";

      /* Put the last value in the FwdT_Obj (i.e. in the fwd_T_curve ) */
      fwd_spread = a[2];
      fwd_T_rate = fwd_R_rate + fwd_spread;
      putoncurve(fwdT, tenor_date, cum_time, fwd_spread, fwd_R_rate, index);

      /* Memory free */
      free_Swap(tec_swap);

      /* One more spread initialised in the FwdT_Obj (i.e. in the fwd_T_curve )
       */
      index += 1;

    } /* END for (i = first_swap_index; i<num_swap; i++)  */

    free_Arg(tec_arg);

  } /* END if (num_swap >0 ) */
  else {
    solve_flag = "NO SPREAD INPUT";
  }
  return (solve_flag);

} /* END FwdTEC_init() */

/* ------------------------------------------------------------------------ */

/* ======================================================================== */
/*           Add a CMT market to the market list			    */
/*           Initialises all that has to be initialised: 		    */
/*           	   CMT_Obj						    */
/*           	   FwdT_Obj						    */
/* ======================================================================== */

Err swp_f_strip_CMT(String fwd_T_name, String *tenor_names,
                    double *spread_rates, int num_spreads, Date cmt_today_long,
                    String yc_name, String cmt_index_string, String vc_name,
                    String mkt_name,
                    Err (*GetVol)(long, long, double, double *),
                    CMT_Param **cmt_param) {
  Err err;
  CMT_Obj *cmt;
  FwdT_Obj *fwdT;
  SrtCurvePtr yc_crv;
  SrtCurveList *crv_list;
  static int first_call = 0;
  CMTCode cmt_code;
  int clcn_date;

  /* Gets the yc_crv associated to the yc_name */
  yc_crv = lookup_curve(yc_name);
  if (!yc_crv)
    return serror("Could not find yield curve %s", yc_name);
  /* Sets calculation date to today according to the Yield Curve*/
  clcn_date = get_clcndate_from_yldcrv(yc_crv);

  /* Initiliases and gets default cmt_param if not initialized */
  if (!(*cmt_param)) {
    if (err = interp_cmt_string(cmt_index_string, &cmt_code))
      return (err);

    *cmt_param = init_CMT_Param(yc_name, vc_name, mkt_name, GetVol, cmt_code);
  }

  /* Memory allocation*/
  cmt = new_CMT_Obj();
  fwdT = new_FwdT_Obj();

  /* Init th CMT_Obj using the CMT-CMS spread curve (given by tenors)*/
  if (err = CMT_init(cmt, cmt_today_long, *cmt_param, tenor_names, spread_rates,
                     num_spreads)) {
    return (err);
  }

  /* Sets fwd treasury curve interpolation method: linear in r or r*t */
  FwdT_Iset(fwdT, FwdT_INTERP_METHOD, (*cmt_param)->interp_method);

  /* Pick up the next available market (in the static list) and attach
     the CMT_Obj (initialised) and the FwdT_Obj (allocated) to this market*/
  crv_list = get_curve_list();

  if (err = swp_f_addcmtcrv_to_list(crv_list, fwd_T_name, "CMT", clcn_date,
                                    yc_name, fwdT, cmt, *cmt_param)) {
    return (err);
  }

  /* Strips the CMT-CMS spread curve to get the fwd treasury curve */
  if ((*cmt_param)->cmt_code != TEC10) {
    if (err = FwdT_init(fwd_T_name, cmt, fwdT, *cmt_param)) {
      return (err);
    }
  } else {
    if (err = FwdTEC_init(fwd_T_name, cmt, fwdT, *cmt_param)) {
      return (err);
    }
  }

  /* Return a NULL success string */
  return NULL;
}

/* ------------------------------------------------------------------------ */
/*           Fills in a fwdT object from a full spread curve
 */
/* ======================================================================== */

/* here  , the fwd_spreads are negative  and are in full value */

static Err FwdT_fill(String fwd_spread_name, double *fwd_spread_date,
                     double *fwd_spread_value, int num_fwd_spreads,
                     FwdT_Obj *fwdT, SrtCurvePtr yc_crv, Date clcn_date,
                     Date spot_date, CMT_Param *cmt_param)

{
  int i;
  SwapDP swapdp;
  double fwd_R_rate;
  double cum_time;
  String ccy_str;

  /* Get the Ccy */
  ccy_str = get_curve_ccy(yc_crv);

  /* Memory allocation */
  FwdT_field_set_Dlength(fwdT, num_fwd_spreads + 2);
  FwdT_field_set_Ilength(fwdT, num_fwd_spreads + 2);

  /* Sets a few fields in the FwdT_obj (interp_method  , today  , spot)*/
  FwdT_Iset(fwdT, FwdT_INTERP_METHOD, cmt_param->interp_method);
  FwdT_Dset(fwdT, FwdT_TODAY, clcn_date);
  FwdT_Dset(fwdT, FwdT_SPOT, spot_date);

  for (i = 0; i < num_fwd_spreads; i++) {
    /* Computes the (forward) swap rate (start date = fwd_spread_date) */
    swapdp.start = (long)fwd_spread_date[i];
    swapdp.nfp = (int)cmt_param->cmt_mat * cmt_param->swap_compd;
    swapdp.end = swapdp.nfp;
    swapdp.direction = FWD;
    swapdp.compd = cmt_param->swap_compd;
    swapdp.basis_code = cmt_param->swap_basis_code;
    swapdp.spot_lag = cmt_param->spot_lag;
    swapdp.first_full_fixing = swapdp.start;
    swp_f_ForwardRate_SwapDP(&swapdp, get_curve_name(yc_crv),
                             cmt_param->swap_ref_rate, &fwd_R_rate);

    /* Computes cum_time from spot_date */
    cum_time = (double)(fwd_spread_date[i] - spot_date) * YEARS_IN_DAY;

    /* Puts the element in the FwdT_Obj */
    putoncurve(fwdT, fwd_spread_date[i], cum_time, fwd_spread_value[i],
               fwd_R_rate, i);

  } /* END of for(i=0;i<num_fwd_spreads) loop */

  return NULL;
}

/* ======================================================================== */
/*           Add a fwdT curve to the market list from a spread curve*/
/* ======================================================================== */

Err swp_f_store_fwd_spreads(String fwd_spread_name, double *fwd_spread_dates,
                            double *fwd_spread_values, int num_fwd_spreads,
                            String yc_name, String cmt_index_string,
                            String vc_name, String mkt_name,
                            Err (*GetVol)(long, long, double, double *),
                            CMT_Param *cmt_param) {
  CMTCode cmt_code;
  FwdT_Obj *fwdT = NULL;
  SrtCurveList *crv_list;
  SrtCurvePtr yc_crv;
  int clcn_date;
  int spot_date;
  Err err;

  /* Gets the yc_crv associated to the yc_name */
  yc_crv = lookup_curve(yc_name);
  if (!yc_crv)
    return serror("Could not find yield curve %s", yc_name);

  /* Sets calculation date == today & spot_date according to the Yield Curve*/
  clcn_date = get_clcndate_from_yldcrv(yc_crv);
  spot_date = get_spotdate_from_yldcrv(yc_crv);

  /* Initilizes and gets default cmt_param if not initialized */
  if (!(cmt_param)) {
    if (err = interp_cmt_string(cmt_index_string, &cmt_code))
      return (err);

    cmt_param = init_CMT_Param(yc_name, vc_name, mkt_name, GetVol, cmt_code);
  }

  /* Create a new FwdT_Obj */
  fwdT = new_FwdT_Obj();
  if (!fwdT)
    return serror("Could not allocate space for fwdT_obj");

  /* Fills in the Swap-Treas spread curve */
  if (err = FwdT_fill(fwd_spread_name, fwd_spread_dates, fwd_spread_values,
                      num_fwd_spreads, fwdT, yc_crv, clcn_date, spot_date,
                      cmt_param)) {
    return (err);
  }

  /* Pick up the next available market (in the static list) and attach
     the FwdT_Obj (allocated and filled in) to this market*/
  crv_list = get_curve_list();
  if (err = swp_f_addcmtcrv_to_list(crv_list, fwd_spread_name, "CMT", clcn_date,
                                    yc_name, fwdT, NULL, cmt_param)) {
    return (err);
  }

  /* Returns a success message */
  return NULL;
}

/* -------------------------------------------------------------------------- */

static SrtErr swp_f_addcmtcrv_to_list(SrtCurveList *curve_list,
                                      String fwdT_name, String curve_lbl,
                                      Date clcn_date, String yc_name,
                                      FwdT_Obj *fwdT, CMT_Obj *cmt,
                                      CMT_Param *cmtprm) {
  SrtCurveDesc *curve;
  SrtCMTDesc *cmtcrv;
  SrtErr err;
  SrtCrvPtr yc_crv;

  /* create space for SrtCurveDesc */
  curve = (SrtCurveDesc *)srt_calloc(1, sizeof(SrtCurveDesc));

  /* set up the elements in the underlying object val */
  strcpy(curve->curve_name, fwdT_name);
  /*
          strupper(curve->curve_name);
          strip_white_space(curve->curve_name);
  */

  /* Get the currency from the YC */
  yc_crv = lookup_curve(yc_name);
  if (!yc_crv)
    return serror("Could not find yield curve %s", yc_name);

  strcpy(curve->curve_ccy, get_curve_ccy(yc_crv));

  /* Initialise the label */
  err = srt_f_interp_curve(curve_lbl, &(curve->curve_type));
  if (err)
    return err;

  if (curve->curve_type != CMT_CURVE)
    return (
        serror("Trying to use addCMTcrv_to_list for a %s curve", curve_lbl));

  strcpy(curve->curve_label, curve_lbl);

  /* creation of  CMT_CurveDesc objects */
  cmtcrv = (SrtCMTDesc *)srt_calloc(1, sizeof(SrtCMTDesc));
  cmtcrv->clcn_date = clcn_date;
  strcpy(cmtcrv->fwdT_name, fwdT_name);
  cmtcrv->fwdT = fwdT;
  cmtcrv->cmt = cmt;
  strcpy(cmtcrv->yc_name, yc_name);
  cmtcrv->cmtprm = cmtprm;

  curve->curve_desc = cmtcrv;

  if ((err = srt_f_lstins(curve_list, fwdT_name, 0.0, OBJ_PTR_UND, curve,
                          &srt_f_curvevalfree, &(curve->curve_ticker))) != NULL)
    return (serror("Error in initialising %s curve", curve->curve_label));

  return NULL;
} /* END swp_f_crvaddcmtcrv(...) */

/* -----------------------------------------------------------------------------------
 */

Err swp_f_strip_treas(String fwd_T_name, String *tenor_names,
                      double *spread_rates, int num_spreads,
                      Date cmt_today_long, String yc_name,
                      String cmt_index_string, String vc_name, String mkt_name,
                      Err (*GetVol)(long, long, double, double *),
                      CMT_Param **cmt_param) {
  Err solve_flag = NULL;
  Err err;
  Arg_Obj *cmt_arg;
  Swap_Obj *cmt_swap;
  Leg_Obj *cmt_leg;
  CMT_Obj *cmt;
  FwdT_Obj *fwdTtemp, *fwdTresult;
  SrtCrvPtr yc_crv;
  SrtCurveList *crv_list;
  static int first_call = 0;
  CMTCode cmt_code;
  int clcn_date;
  int toalloc;
  int num_full_period;
  int first_swap_index = 0;
  Date swap_end;
  double first_full_fixing;
  double fwd_spread, roll_spread_lock;
  double cum_time;
  double cmt_swap_temp_price, cmt_swap_price, cmt_cmt_price;
  Date spot, today, tenor_date, last_fixing, start_date;
  double fwd_T_rate, fwd_R_rate, a[3], b[3], yans;
  int i, index, niter;
  int count, num_swap, swap_len;
  double nstop;
  SwapDP swapdp;
  DateList datelist;
  SrtCurvePtr
      crv; /* Can be either a cmt_crv or a yc_crv depending on last use*/
  String ccy_str;

  /* Gets the yc_crv associated to the yc_name */
  yc_crv = lookup_curve(yc_name);
  if (!yc_crv)
    return serror("Could not find yield curve %s", yc_name);
  /* Sets calculation date to today according to the Yield Curve*/
  clcn_date = get_clcndate_from_yldcrv(yc_crv);

  /* Get the Ccy of the Market */
  ccy_str = get_curve_ccy(yc_crv);

  /* Initializes and gets default cmt_param if not initialized */
  if (!(*cmt_param)) {
    if (err = interp_cmt_string(cmt_index_string, &cmt_code))
      return (err);

    *cmt_param = init_CMT_Param(yc_name, vc_name, mkt_name, GetVol, cmt_code);
  }

  /* Memory allocation*/
  cmt = new_CMT_Obj();

  /* Init the CMT_Obj using the CMT-CMT spread curve (given by tenors)
     In fact the spread is a constant fwd spread */
  if (err = CMT_init(cmt, cmt_today_long, *cmt_param, tenor_names, spread_rates,
                     num_spreads)) {
    return (err);
  }

  /* we would need two fwdT curves : one which contains the proper result  ,
  the other one for the computation */
  fwdTtemp = new_FwdT_Obj();   /* the one for computation */
  fwdTresult = new_FwdT_Obj(); /* the one for proper result */

  /* Sets fwd treasury curve interpolation method: linear in r or r*t */
  /* we do not really care about the interpolation method for the temp result */
  /* it is just going to be a constant spread across the trade*/
  FwdT_Iset(fwdTtemp, FwdT_INTERP_METHOD, LIN_RT);
  FwdT_Iset(fwdTresult, FwdT_INTERP_METHOD, (*cmt_param)->interp_method);

  /* Pick up the next available market (in the static list) and attach
  the CMT_Obj (initialised) and the FwdT_Obj (allocated) to this market*/
  crv_list = get_curve_list();

  if (err = swp_f_addcmtcrv_to_list(crv_list, fwd_T_name, "CMT", clcn_date,
                                    yc_name, fwdTresult, cmt, *cmt_param)) {
    return (err);
  }

  crv = lookup_curve(fwd_T_name);

  /* Stores dates used everywhere in the code */
  spot = (Date)CMT_Dget(cmt, CMT_SPOT);
  today = (Date)CMT_Dget(cmt, CMT_TODAY);
  cmt_code = CMT_Dget(cmt, CMT_CODE);

  /* Memory allocation for the FwdT_Objects (use  more than necesary)*/
  /* Initialisation of some fields */
  swap_len = CMT_field_Dlength(cmt, CMS_CMT_SPREAD_RATE);
  toalloc = 2 * swap_len + 2 * (CMT_Iget(cmt, CMT_SPOT_EXISTS));

  FwdT_field_set_Dlength(fwdTresult, toalloc);
  FwdT_field_set_Ilength(fwdTresult, toalloc);
  FwdT_Dset(fwdTresult, FwdT_TODAY, today);
  FwdT_Dset(fwdTresult, FwdT_SPOT, spot);
  FwdT_Dset(fwdTresult, FwdT_CMTCODE, cmt_code);

  FwdT_field_set_Dlength(fwdTtemp, toalloc);
  FwdT_field_set_Ilength(fwdTtemp, toalloc);
  FwdT_Dset(fwdTtemp, FwdT_TODAY, today);
  FwdT_Dset(fwdTtemp, FwdT_SPOT, spot);
  FwdT_Dset(fwdTtemp, FwdT_CMTCODE, cmt_code);

  index = 0;

  /* Set the first one only for the result curve */
  if (CMT_Iget(cmt, CMT_SPOT_EXISTS) == 1) {
    fwd_spread = -CMT_Dget(cmt, CMT_SPOT_SPREAD);
    start_date = (Date)CMT_Dget(cmt, CMT_SPOT_DATE);
    cum_time = (start_date - spot) * YEARS_IN_DAY;
    tenor_date = bus_date_method(start_date, MODIFIED_SUCCEEDING);

    /* Computes the (forward) swap rate (start date = spot_date of CMT swap) */
    swapdp.start = start_date;
    swapdp.nfp = (int)(*cmt_param)->cmt_mat * (*cmt_param)->swap_compd;
    swapdp.end = swapdp.nfp;
    swapdp.direction = FWD;
    swapdp.compd = (*cmt_param)->swap_compd;
    swapdp.basis_code = (*cmt_param)->swap_basis_code;
    swapdp.spot_lag = (*cmt_param)->spot_lag;
    swapdp.first_full_fixing = swapdp.start;
    swp_f_ForwardRate_SwapDP(&swapdp, yc_name, (*cmt_param)->swap_ref_rate,
                             &fwd_R_rate);

    /* Sets the initial value for the fwd_T_rate in the curve:
    please note that FwdT_NUMSPREAD is set to  index+1 */
    putoncurve(fwdTresult, tenor_date, cum_time, fwd_spread, fwd_R_rate, index);
  }

  /* Index is used to store the number of spreads already in the curve */
  index = FwdT_Iget(fwdTresult, FwdT_NUMSPREAD);

  /* Initialises in the CMT ARG what is independent of the instrument */
  init_Arg(&cmt_arg, GENERATE_SWAP);
  Arg_Dset(cmt_arg, ARG_NOTIONAL, 1000000000.00);
  Arg_Dset(cmt_arg, ARG_INITIAL_NOT, 1);
  Arg_Dset(cmt_arg, ARG_FINAL_NOT, 1);
  Arg_Dset(cmt_arg, ARG_STRIKE, 0.0);
  Arg_set(cmt_arg, ARG_COMPD, (*cmt_param)->cmt_freq);
  Arg_set(cmt_arg, ARG_BASIS_CODE, (*cmt_param)->cmt_basis_code);
  Arg_set(cmt_arg, ARG_DATE_DIR, BKWD);
  Arg_Dateset(cmt_arg, ARG_START, spot);
  Arg_Dateset(cmt_arg, ARG_TODAY, today);
  Arg_Dateset(cmt_arg, ARG_SPOT, spot);
  Arg_Dset(cmt_arg, ARG_SPREAD, 0);

  /* Lets do the loop for all the instruments */
  num_swap = CMT_Iget(cmt, CMT_NUMSWAP);
  if (num_swap > 0) {
    for (i = first_swap_index; i < num_swap; i++) {
      first_full_fixing = CMT_field_Dget(cmt, CMT_SWAP_FIRST_FULL_FIXING, i);
      num_full_period = CMT_field_Iget(cmt, CMT_SWAP_NUM_FULL_PERIOD, i);
      swap_end = (Date)CMT_field_Dget(cmt, CMT_SWAPEND_DATE, i);

      Arg_Iset(cmt_arg, ARG_NUM_FULL_PERIOD, num_full_period);
      Arg_Dset(cmt_arg, ARG_FIRST_FULL_FIXING, first_full_fixing);
      Arg_Dateset(cmt_arg, ARG_END, swap_end);

      /* Computes the last fixing date in the swap
      (last payment is made at swap_end) */
      swapdp.end = swap_end;
      swapdp.nfp = swapdp.end;
      swapdp.direction = BKWD;
      swapdp.compd = (*cmt_param)->cmt_freq;
      swapdp.basis_code = (*cmt_param)->cmt_basis_code;
      swapdp.spot_lag = (*cmt_param)->spot_lag;
      swapdp.start = today;
      swapdp.first_full_fixing = swapdp.start;

      datelist = SwapDP_to_DateList(&swapdp, MODIFIED_SUCCEEDING);
      last_fixing = datelist.date[datelist.len - 2];
      free(datelist.date);

      /* Computes the forward swap rate (start date = last_fixing of CMT swap)
       */
      swapdp.nfp = (int)(*cmt_param)->cmt_mat * (*cmt_param)->swap_compd;
      swapdp.end = swapdp.nfp;
      swapdp.direction = FWD;
      swapdp.compd = (*cmt_param)->swap_compd;
      swapdp.basis_code = (*cmt_param)->swap_basis_code;
      swapdp.spot_lag = (*cmt_param)->spot_lag;
      swapdp.start = last_fixing;
      swapdp.first_full_fixing = swapdp.start;
      swp_f_ForwardRate_SwapDP(&swapdp, yc_name, (*cmt_param)->swap_ref_rate,
                               &fwd_R_rate);

      /* Gets the value of the CMS-CMT spread from the spreadsheet (thru
       * CMT_Obj)*/
      /* it is in fact the CMT-CMT spread */
      roll_spread_lock = CMT_field_Dget(cmt, CMS_CMT_SPREAD_RATE, i);

      /* Sets tenor_date to be the right date in the fwd_T_curve */
      cum_time = (last_fixing - spot) * YEARS_IN_DAY;
      tenor_date = bus_date_method(last_fixing, MODIFIED_SUCCEEDING);

      /* Sets the initial guess for the fwd_T_rate in the curve temp */
      fwd_spread = -roll_spread_lock;
      fwd_T_rate = fwd_R_rate + fwd_spread;

      /* Be careful index is always 0 */
      /* like that we do not have to reinitialise the curve each time */
      putoncurve(fwdTtemp, tenor_date, cum_time, fwd_spread, fwd_R_rate, 0);

      /* Build the CMT-spread / LIBOR swap (the vols depend on the fwdT rate */
      cmt_swap = generate_swap(cmt_arg, CMT_MARGIN_FLOATING_SWAP);

      /* Extracts the CMT-spread Leg of the CMS - Libor swap */
      cmt_leg = Swap_field_get(cmt_swap, SWAP_LEG, 0);

      /* Get the market and change the fwd curve */
      crv = lookup_curve(fwd_T_name);
      ((SrtCMTDesc *)crv->curve_desc)->fwdT = fwdTtemp;

      /* Computes and store fwd swap rates  , and swap vols */
      leg_cmt_populate(cmt_leg, crv);

      /* Initial price of the CMT leg */
      cmt_swap_temp_price = value_swap(crv, cmt_swap);

      /* STARTS ITERATION TO SOLVE FOR A PAR VALUE OF THE CMT_CMT SWAP */

      /* Sets a few useful parameters */
      niter = 5;
      yans = 0.0;
      count = 0;
      nstop = 0.0;

      /* Sets an initial guess for the new fwd treasury rate */
      fwd_spread = -roll_spread_lock;
      fwd_T_rate = fwd_R_rate + fwd_spread;
      a[0] = fwd_spread;

      /* Sets the initial guess for the fwd_T_rate in the curve */
      putoncurve(fwdTresult, tenor_date, cum_time, fwd_spread, fwd_R_rate,
                 index);

      /* Change the fwd curve */
      ((SrtCMTDesc *)crv->curve_desc)->fwdT = fwdTresult;

      /* Initial price of the CMT leg wit the result curve */
      cmt_swap_price = value_swap(crv, cmt_swap);

      /* Initial price of the CMT-CMT swap */
      cmt_cmt_price = cmt_swap_price - cmt_swap_temp_price;

      /* Test only to speed up the computation when we have a flat curve */
      if (fabs(cmt_cmt_price) > PRECISION) {
        b[0] = cmt_cmt_price;

        /* Does this again for two shifts */
        fwd_spread = fwd_spread - 0.0001;
        fwd_T_rate = fwd_R_rate + fwd_spread;
        a[1] = fwd_spread;
        putoncurve(fwdTresult, tenor_date, cum_time, fwd_spread, fwd_R_rate,
                   index);
        cmt_swap_price = value_swap(crv, cmt_swap);
        cmt_cmt_price = cmt_swap_price - cmt_swap_temp_price;
        b[1] = cmt_cmt_price;

        fwd_spread = fwd_spread - 0.0001;
        fwd_T_rate = fwd_R_rate + fwd_spread;
        a[2] = fwd_spread;

        /* Calls NEWTON for solving */
        while ((count < MAXITER) && (nstop < 1.0)) {
          fwd_spread = a[2];
          fwd_T_rate = fwd_R_rate + fwd_spread;
          putoncurve(fwdTresult, tenor_date, cum_time, fwd_spread, fwd_R_rate,
                     index);
          cmt_swap_price = value_swap(crv, cmt_swap);
          cmt_cmt_price = cmt_swap_price - cmt_swap_temp_price;
          b[2] = cmt_cmt_price;
          newton(yans, niter, a, b, &nstop);
          count += 1;
        }

        /* Checks convergence */
        if (nstop < 1.0)
          solve_flag = "UNABLE TO SOLVE";
      } else {
        /* initialise the spread which is going to be copied into the */
        /* fwd result */
        a[2] = a[0];
      }

      /* Put the last value in the FwdT_Obj (i.e. in the fwd_T_curve ) */
      fwd_spread = a[2];
      fwd_T_rate = fwd_R_rate + fwd_spread;
      putoncurve(fwdTresult, tenor_date, cum_time, fwd_spread, fwd_R_rate,
                 index);

      /* Memory free */
      free_Swap(cmt_swap);

      /* One more spread initialised in the FwdT_Obj (i.e. in the fwd_T_curve )
       */
      index += 1;

    } /* END for (i = first_swap_index; i<num_swap; i++)  */

    free_Arg(cmt_arg);
  }

  free_FwdT_Obj(fwdTtemp);

  return NULL;
}