/* ===========================================================================
   FILENAME:   srt_f_mdlinistp.cxx

   PURPOSE:    All the functions needed to attach at a time step all the
               informaton required for a discretisation that is not SIMULATION
                           dependent (sigma        , taus        , rho        ,
   ...)

   FUNCTIONS:  - srt_f_irministp:      for any interest rate underlying
               - srt_f_lgministp:      specifically for LGM
               - srt_f_betaetainistp:  for the EtaBeta model
               - srt_f_loginistp:      for BlackScholes type underlyings
               - srt_f_basicinistp:    for deterministic underlyings

   DESCRIPTION:  all the functions allocate space for a Srt...TmInf and
                 attach one (populated) to each time step
   ===========================================================================
 */

/* =========================================================================

   IRM function: used for CHE and LGM

   ======================================================================== */

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_cheybetadynamics.h"
#include "srt_h_fwdcurve.h"
#include "srt_h_mdlinistp.h"
#include "srt_h_powermodel.h"
#include "srt_h_repo_obj.h"
#include "srt_h_stpcorata.h"

/* ------------------------------------------------------------------------------
   Sets in every single time step of the discretisation the deterministic
   information relevant for the discrete scheme used (Tree OR MC)
   (und_index corresponds to the underlying und_index        , needed to extract
  the correct tminf from the stp)
  -----------------------------------------------------------------------------
*/

Err srt_f_irministp(
    SrtStpPtr stp, /* The linked list of time steps */
    SrtUndPtr und, /* The interest rate underlying */
    int und_index, /* The underlying index (0 if the domestic) */
    SrtUndPtr
        numeraire_und, /* THe Numeraire Used for Discounting (quanto's...) */
    SrtUndInfo *und_info) {
  SrtStpPtr top;
  Err err;
  SrtIRMTmInf *tminf, *prvtminf;
  TermStruct *ts;
  Ddate today;
  String yc_name, und_name;
  SrtCurvePtr yldcrv;
  SrtMdlDim mdl_dim;
  SrtMdlType mdl_type;

  /* Get the relevant information from the underlying */
  prvtminf = NULL;
  und_name = get_underlying_name(und);

  err = get_underlying_mdldim(und, &mdl_dim);
  if (err)
    return err;
  err = get_underlying_mdltype(und, &mdl_type);
  if (err)
    return err;

  yc_name = get_discname_from_underlying(und);
  err = get_underlying_ts(und, &ts);
  if (err)
    return err;

  /* Go to the top of the steps linked list (i.e. today) */
  top = stp = gototop(stp);

  /* Allocate space for a time info pointer on each element of the step list */
  err = srtstpalloc(top, sizeof(SrtIRMTmInf), 0, und_index);
  if (err)
    return err;

  /* We need the Yield Curve for the Y_t_at_T call */
  yldcrv = lookup_curve(yc_name);
  today = get_clcndate_from_yldcrv(yldcrv);

  while (stp) {
    /* Get the time info */
    tminf = (SrtIRMTmInf *)stp->tminf[und_index];

    /* Discount factor from today to time step */
    if (mdl_type != VASICEK)
      tminf->df = swp_f_df(today, stp->ddate, yc_name);
    else {
      err = srt_f_vasicek_discount_factor(stp->date, und, &tminf->df);
      if (err)
        return err;

      tminf->mean_rev_level = find_mean_rev_level(stp->time, ts);
      tminf->vasicek_init_cond = find_vasicek_init_cond(stp->time, ts);
    }

    /* Sets the instantaneous forward rate for that date */
    if (stp->next) {
      if (mdl_type != VASICEK)
        sam_get(tminf->fwd_sam, und_index, F_0_t) =
            swp_f_zr(stp->ddate, stp->next->ddate, yc_name);
      else {
        err = srt_f_vasicek_cont_fwd_zr(
            stp->ddate, stp->next->ddate, und,
            &sam_get(tminf->fwd_sam, und_index, F_0_t));
        if (err)
          return err;
      }
    } else {
      if (mdl_type != VASICEK)
        sam_get(tminf->fwd_sam, und_index, F_0_t) =
            swp_f_zr(stp->ddate, stp->ddate + 14.0, yc_name);
      else {
        err = srt_f_vasicek_cont_fwd_zr(
            stp->ddate, stp->ddate + 14.0, und,
            &sam_get(tminf->fwd_sam, und_index, F_0_t));
        if (err)
          return err;
      }
    }

    switch (mdl_dim) {
      /* One factor models: put the vols        , mean reversion        , ... */
    case ONE_FAC:
      if (stp->next)
        tminf->ev.onef.sig2 = find_sig2_interp(stp->time, stp->next->time, ts);
      else
        tminf->ev.onef.sig2 =
            find_sig2_interp(stp->time, stp->time + 14.0 / 365.0, ts);

      tminf->ev.onef.sig = sqrt(tminf->ev.onef.sig2);
      tminf->ev.onef.tau = find_tau(stp->time, ts);
      tminf->ev.onef.lambda = 1.0 / tminf->ev.onef.tau;

      if (mdl_type == CHEY_BETA || mdl_type == CHEY_BETA_STOCH_VOL) {
        tminf->ev.onef.beta = find_beta(stp->time, ts);
      }

      if (mdl_type == LGM_STOCH_VOL || mdl_type == CHEY_STOCH_VOL ||
          mdl_type == CHEY_BETA_STOCH_VOL) {
        tminf->vovol = find_vovol(stp->time, ts);
        tminf->vovol_sqr = tminf->vovol * tminf->vovol;
        tminf->rho = find_rho(stp->time, ts);
      }

      break;
      /* Two factor models: put the vols        , mean reversion        ,... */
    case TWO_FAC:
      if (stp->next)
        err = get_2f_sig2_rho_interp(stp->time, stp->next->time, ts,
                                     &tminf->ev.twof[0].sig2,
                                     &tminf->ev.twof[1].sig2, &tminf->correl_x);
      else
        err = get_2f_sig2_rho_interp(stp->time, stp->time + 14.0 / 365.0, ts,
                                     &tminf->ev.twof[0].sig2,
                                     &tminf->ev.twof[1].sig2, &tminf->correl_x);

      tminf->ev.twof[0].sig = sqrt(tminf->ev.twof[0].sig2);
      tminf->ev.twof[1].sig = sqrt(tminf->ev.twof[1].sig2);
      err = find_2f_tau(stp->time, ts, &tminf->ev.twof[0].tau,
                        &tminf->ev.twof[1].tau);
      tminf->ev.twof[0].lambda = 1.0 / tminf->ev.twof[0].tau;
      tminf->ev.twof[1].lambda = 1.0 / tminf->ev.twof[1].tau;

      find_tf_beta(stp->time, ts, &(tminf->ev.twof[0].beta),
                   &(tminf->ev.twof[1].beta));

      break;
    default:
      break;
    }
    /* For Cheyette type of models        , evolves a "deterministic" Phi:
            very important for tree geometry */
    if (prvtminf && mdl_dim == ONE_FAC &&
        ((mdl_type == CHEY) || (mdl_type == CHEY_BETA))) {
      srt_f_CheyBeta_drift_at_sam(mdl_type, stp->prev, &(prvtminf->fwd_sam),
                                  &(tminf->fwd_sam), und_index);
    } else {
      switch (mdl_dim) {
      case ONE_FAC:
        sam_get(tminf->fwd_sam, und_index, PHI) = 0.0;
        break;
      case TWO_FAC:
        sam_get(tminf->fwd_sam, und_index, PHI1) = 0.0;
        sam_get(tminf->fwd_sam, und_index, PHI2) = 0.0;
        sam_get(tminf->fwd_sam, und_index, CROSSPHI) = 0.0;
        break;
      default:
        break;
      }
    }

    /* For DF from date to next date in discretisation (reconstruction formula)
     */
    und = lookup_und(und_name);
    if (is_model_Cheyette_type(mdl_type)) {
      if (stp->next)
        err = Y_T_at_t_param(stp->ddate, &stp->next->ddate, 1, und, &tminf->yp);
    }
    und = lookup_und(und_name);

    /* Move on to the next step */
    prvtminf = tminf;
    stp = stp->next;

  } /* END of wile loop on stp */

  /* Reset market pointer to interest rate und for consistency */
  und = lookup_und(und_name);

  /* If this is the linear gauss markov model (which we regard as a
          subset of the cheyette model) then we can precompute phi and
          several other variables        , which we now do */

  /* Model specific functions */
  if ((mdl_type == LGM) || (mdl_type == ETABETA) || (mdl_type == VASICEK)) {
    err = srt_f_lgministp(top, ts, und_index, mdl_dim, numeraire_und, und_info);
    if (err)
      return err;
  }

  if (mdl_type == ETABETA) {
    err = srt_f_betaetainistp(top, ts, und_index, mdl_dim);
    if (err)
      return err;
  }

  /* Compute the Quanto adjustment and stores it in the step */
  err = srt_f_steps_set_quanto_adjustment(top, und_info, und_index, und,
                                          numeraire_und);
  if (err)
    return err;

  /* Return a success message */
  return err;

} /* Err srt_f_irministp(..) */

/* ----------------------------------------------------------------------- */

/* ======================================================================

        LGM specific functions

   ====================================================================== */

/* ------------------------------------------------------------------------
   FUNCNAME        :srt_f_lgministp (was calc_Phi_LGM)

   AUTHOR          :G.Amblard

   DESCRIPTION     :puts phi (which is deterministic)        ,F        ,G ,Psi
on nodes on step for LGM model; these values are computed from the ts on mdl;
     really the ts should be located on a SrtUndPtr.  Note that we regard the
LGM as a subset of the Cheyette model (with beta = 0)        , justifying the
use of Cheyette structures; this function should be called after srt_f_cheinistp

(old:)
 modified yet again by K Chau for the use of new G_H functions

(new)
  AMENDMENTS      :
  Reference       :
  Author          :E.Auld
  Date            :sep 94
  Description     :modified to use SrtCheTmInf structures

  Reference       :can't use find sig interp
  Author          :E.Auld
  Date            :1 Nov 94
  Description     :changed from using firsttminf->sig to sig0 for computing F
  ----------------------------------------------------------------------------
*/

Err srt_f_lgministp(SrtStpPtr stp, TermStruct *ts, int und_index,
                    SrtMdlDim mdl_dim,
                    SrtUndPtr numeraire_und, /* THe Numeraire Used for
                                          Discounting (quanto's...) */
                    SrtUndInfo *und_info) {
  SrtStpPtr cur, first;
  SrtIRMTmInf *tminf, *prvtminf;
  SrtTFTSVec F, Psi;
  SrtTFTSMat G, H;
  Err err = NULL;

  double dt;

  cur = first = gototop(stp);
  prvtminf = NULL;

  /* init. of the integrale of PHI stochastic rate & volatility gamma smile
   * model  */
  tminf = cur->tminf[und_index];

  while (cur) {
    tminf = cur->tminf[und_index];
    switch (mdl_dim) {
    case ONE_FAC:
      tminf->ev.onef.F = F_func(cur->time, ts);
      tminf->ev.onef.Psi = Psi_func(cur->time, ts);
      /* Commented out for the moment (OVE)
                      tminf->ev.onef.J = J_func(cur->time        ,ts);
      */
      G_H_func(cur->time, ts, &tminf->rf.onef.G, &tminf->rf.onef.H);

      sam_get(tminf->fwd_sam, und_index, PHI) =
          tminf->rf.onef.G * tminf->ev.onef.F * tminf->ev.onef.F;

      /* compute the integral of PHI for stochastic rate & volatility gamma
       * smile model */
      if (cur->prev)
        dt = (cur->time - cur->prev->time);
      else
        dt = cur->time;

      if (prvtminf)
        tminf->rf.onef.stdev_x =
            tminf->ev.onef.F * sqrt(tminf->rf.onef.G - prvtminf->rf.onef.G);

      break;
    case TWO_FAC:
      err = get_2f_F_funcs(cur->time, ts, &F);
      tminf->ev.twof[0].F = F[0];
      tminf->ev.twof[1].F = F[1];
      err = get_2f_Psi_funcs(cur->time, ts, &Psi);
      tminf->ev.twof[0].Psi = Psi[0];
      tminf->ev.twof[1].Psi = Psi[1];

      err = get_2f_G_funcs(cur->time, ts, &G);
      tminf->rf.twof[0][0].G = G[0][0];
      tminf->rf.twof[0][1].G = G[0][1];
      tminf->rf.twof[1][0].G = G[1][0];
      tminf->rf.twof[1][1].G = G[1][1];

      err = get_2f_H_funcs(cur->time, ts, &H);
      tminf->rf.twof[0][0].H = H[0][0];
      tminf->rf.twof[0][1].H = H[0][1];
      tminf->rf.twof[1][0].H = H[1][0];
      tminf->rf.twof[1][1].H = H[1][1];

      sam_get(tminf->fwd_sam, und_index, PHI1) = G[0][0] * F[0] * F[0];
      sam_get(tminf->fwd_sam, und_index, PHI2) = G[1][1] * F[1] * F[1];
      sam_get(tminf->fwd_sam, und_index, CROSSPHI) = G[1][0] * F[1] * F[0];

      if (prvtminf) {
        tminf->rf.twof[0][0].stdev_x =
            F[0] * sqrt(G[0][0] - prvtminf->rf.twof[0][0].G);
        tminf->rf.twof[1][1].stdev_x =
            F[1] * sqrt(G[1][1] - prvtminf->rf.twof[1][1].G);
        /* Here        , [1][0].stedev_r is the correlation        , no the
         * covariance...*/
        tminf->rf.twof[1][0].stdev_x =
            (G[0][1] - prvtminf->rf.twof[0][1].G) /
            sqrt((G[0][0] - prvtminf->rf.twof[0][0].G) *
                 (G[1][1] - prvtminf->rf.twof[1][1].G));
        tminf->rf.twof[0][1].stdev_x = tminf->rf.twof[1][0].stdev_x;
      }
      break;
    default:
      break;
    } /* END switch (mdl_dim) */

    cur = cur->next;
    prvtminf = tminf;

  } /* END  cur loop on the SrtStpPtr */

  /* Returns a success message */
  return NULL;

} /* END function Err srt_f_lgministp */

/* ========================================================================= */
/* NAB for Beta-Eta model THIS IS NOT CORRECT */

/* -------------------------------------------------------------------------
   THe integrals needed for Beta-Eta discretisations */
Err srt_f_betaetainistp(SrtStpPtr stp, TermStruct *ts, int und_index,
                        SrtMdlDim mdl_dim) {
  SrtStpPtr cur, first;
  SrtIRMTmInf *tminf, *prvtminf;
  Err err = NULL;
  double s, theta;
  double **M;

  cur = first = gototop(stp);
  prvtminf = NULL;

  while (cur) {
    tminf = cur->tminf[und_index];

    tminf->ev.onef.F = F_func(cur->time, ts);

    G_H_func(cur->time, ts, &tminf->rf.onef.G, &tminf->rf.onef.H);

    tminf->lambda_t = Psi_func(cur->time, ts);
    tminf->zeta_t = Zeta_func(cur->time, ts);
    tminf->the_beta = find_beta(cur->time, ts);
    tminf->eta = find_eta(cur->time, ts);

    s = tminf->the_beta * tminf->the_beta * tminf->zeta_t;
    theta = tminf->lambda_t / tminf->the_beta;
    M = find_M_eta_beta(cur->time, ts);

    tminf->A_t = M_eta_beta_func(s, theta, tminf->eta, M);

    cur = cur->next;
    prvtminf = tminf;

  } /* END  cur loop on the SrtStpPtr */

  return NULL;
} /* end function srt_f_betaetainistp */

/* -----------------------------------------------------------------------------
 */

/* ======================================================================

        Black Scholes specific functions

   ====================================================================== */

Err srt_f_loginistp(SrtStpPtr stp, SrtUndPtr und, int und_index,
                    SrtUndPtr numeraire_und, SrtUndInfo *und_info) {
  SrtStpPtr top;
  Err err;
  SrtLogTmInf *tminf, *nxttminf;
  TermStruct *ts;
  String crv3 = NULL, crv2 = NULL, crv1 = NULL, und_name;
  String dom_und_name;
  String for_und_name;
  SrtMdlType mdl_type;
  SrtUndPtr dom_und, for_und;
  SrtUnderlyingType und_type, obj_type1, obj_type2, obj_type3;
  SrtDvdObj *dvd_obj, *repo_obj;
  double spot;
  Ddate today, spotdate;
  SrtCurvePtr yldcrv;
  SrtCurvePtr dvdcrv;
  SrtCurvePtr repocrv;
  double repo, nxt_repo;

  und_name = get_underlying_name(und);
  und_type = get_underlying_type(und);
  if (und_type == EQUITY_UND)
    spot = get_spot_from_eqund(und);
  else if (und_type == FOREX_UND)
    spot = get_spot_from_fxund(und);

  /* Gets the term structure */
  err = get_underlying_ts(und, &ts);
  if (err)
    return err;

  /* Gets the model type (for stochastic drfit cases) */
  err = get_underlying_mdltype(und, &mdl_type);
  if (err)
    return (err);

  /* Gets the discount and growth curves names */
  if (mdl_type == FX_STOCH_RATES) {
    /* Get the domestic underlying name and pointer to */
    dom_und_name = get_domname_from_fxund(und);
    dom_und = lookup_und(dom_und_name);
    if (!dom_und)
      return serror("Could not find %s underlying", dom_und_name);

    /* Get the discount curve from the domestic underlying */
    crv1 = get_discname_from_underlying(dom_und);

    /* Get the foreign underlying name and pointer to */
    for_und_name = get_forname_from_fxund(und);
    for_und = lookup_und(for_und_name);
    if (!for_und)
      return serror("Could not find %s underlying", for_und_name);

    /* Get the discount curve for the foreign underlying */
    crv2 = get_discname_from_underlying(for_und);

    /* Get the spread curve */
    crv3 = get_repo_name_from_underlying(und);

  } else if ((mdl_type == EQ_STOCH_RATES) ||
             (mdl_type == EQ_STOCH_RATES_SRVGS)) {
    /* Get the domestic underlying name and pointer to */
    dom_und_name = get_discname_from_underlying(und);
    dom_und = lookup_und(dom_und_name);
    if (!dom_und)
      return serror("Could not find %s underlying", dom_und_name);

    /* Get the discount curve from the domestic underlying */
    crv1 = get_discname_from_underlying(dom_und);
    crv2 = get_dividend_name_from_underlying(und);
    crv3 = get_repo_name_from_underlying(und);
  } else if ((mdl_type == BLACK_SCHOLES)) {

    /* Get the curves names form the underlying */
    crv1 = get_discname_from_underlying(und);
    crv2 = get_dividend_name_from_underlying(und);
    crv3 = get_repo_name_from_underlying(und);
  } else
    return serror("Unknown FX type in srt_f_loginistp");

  /* Get the discount curve */
  yldcrv = lookup_curve(crv1);
  obj_type1 = get_curve_type(yldcrv);
  today = get_today_from_curve(yldcrv);

  /* Get the growth curve */
  dvdcrv = lookup_curve(crv2);
  obj_type2 = get_curve_type(dvdcrv);

  /* If second curve is a Fwd_Obj        , insert the spot value to it */
  if (obj_type2 == DVD_CURVE) {
    dvd_obj = get_dvdobj_from_dvdcrv(dvdcrv);
    spotdate = get_clcndate_from_dvdcrv(dvdcrv);
    err = srt_f_dvdobj_insertpoint(dvd_obj, spotdate, 1.0);
    if (err)
      return err;
  }

  if (crv3 && strcmp(crv3, "")) {
    repocrv = lookup_curve(crv3);
    obj_type3 = get_curve_type(repocrv);

    if (obj_type3 == REPO_CURVE) {
      repo_obj = get_repoobj_from_repocrv(repocrv);
      spotdate = get_clcndate_from_repocrv(repocrv);

      err = srt_f_repo_obj_insert_point(repo_obj, spotdate, 1.0);
      if (err)
        return err;
    }
  }

  /* Go to the first time step */
  top = stp = gototop(stp);

  /* Allocate space */
  err = srtstpalloc(top, sizeof(SrtLogTmInf), 0, und_index);
  if (err)
    return err;

  /* Populate srttimestps */

  /* Sets the volatility */
  while (stp->next) {
    tminf = stp->tminf[und_index];

    tminf->quanto_adjustment = 0.0;
    if (und_type == EQUITY_UND) {
      tminf->int_sig_dt = eq_cum_vol_func(stp->next->time, SRT_NO, ts) -
                          eq_cum_vol_func(stp->time, SRT_NO, ts);
      tminf->int_sig2_dt = eq_cum_vol_func(stp->next->time, SRT_YES, ts) -
                           eq_cum_vol_func(stp->time, SRT_YES, ts);

      if ((mdl_type == EQ_STOCH_RATES_SRVGS)) {
        tminf->time = stp->time;

        tminf->sig = find_eq_sig(stp->time, ts);
        tminf->omega = find_eq_omega(stp->time, ts);
        tminf->beta = find_eq_beta(stp->time, ts);
        tminf->gamma = find_eq_gamma(stp->time, ts);

        tminf->basevol = find_eq_basevol(stp->time, ts);
        tminf->voldrift = find_eq_voldrift(stp->time, ts);

        tminf->vovol = find_eq_vovol(stp->time, ts);
        tminf->rho = find_eq_rho(stp->time, ts);
      }

    } else if (und_type == FOREX_UND) {
      tminf->int_sig_dt = fx_cum_vol_func(stp->next->time, SRT_NO, ts) -
                          fx_cum_vol_func(stp->time, SRT_NO, ts);
      tminf->int_sig2_dt = fx_cum_vol_func(stp->next->time, SRT_YES, ts) -
                           fx_cum_vol_func(stp->time, SRT_YES, ts);

      if (crv3 && strcmp(crv3, "")) {
        /* compute tminf->inv_exp_int_spread_dt */
        err = srt_f_repo_obj_repo(repo_obj, stp->next->date, &nxt_repo);
        if (err)
          return err;

        err = srt_f_repo_obj_repo(repo_obj, stp->date, &repo);
        if (err)
          return err;

        tminf->inv_exp_int_spread_dt = nxt_repo / repo;
      } else
        tminf->inv_exp_int_spread_dt = 1.0;
    }

    tminf->sqrt_int_sig2_dt = sqrt(tminf->int_sig2_dt);

    stp = stp->next;
  }

  /* Last step: no volatility should be needed */
  tminf->quanto_adjustment = 0.0;
  tminf = stp->tminf[und_index];
  tminf->int_sig_dt = 0.0;
  tminf->int_sig2_dt = 0.0;
  tminf->sqrt_int_sig2_dt = 0.0;

  /* Sets the discount factors from the first underlying */
  stp = top;
  while (stp) {
    tminf = stp->tminf[und_index];

    tminf->df = swp_f_df(today, stp->ddate, crv1);
    if (tminf->df == SRT_DF_ERROR)
      return serror("Could not compute df for %s curve", crv1);

    stp = stp->next;
  }

  /* Sets the forward (as of today) using the second growth curve */
  dvdcrv = lookup_curve(crv2);
  stp = top;
  while (stp) {
    tminf = stp->tminf[und_index];

    if (obj_type2 == DVD_CURVE) {
      tminf->init_fwd_val =
          spot * srt_f_forward_from_fwdcrv(stp->time, dvdcrv, repocrv) /
          tminf->df;

      tminf->vol2_cum = eq_cum_vol_func(stp->time, SRT_YES, ts);
      tminf->ajust_init_fwd_val =
          tminf->init_fwd_val * exp(-0.5 * tminf->vol2_cum);
    } else

        /* For Equity        , the forward is computed using ONLY the second
           yield curve input */
        if (und_type == EQUITY_UND) {
      tminf->init_fwd_val = spot / tminf->df;
      tminf->vol2_cum = eq_cum_vol_func(stp->time, SRT_YES, ts);
      tminf->ajust_init_fwd_val =
          tminf->init_fwd_val * exp(-0.5 * tminf->vol2_cum);
    } else
        /* For FX        , the forward is computed using the foreign yield curve
           input and the discount one */
        if ((und_type == FOREX_UND) && !(mdl_type == FX_STOCH_RATES)) {
      tminf->init_fwd_val =
          spot * swp_f_df(today, stp->ddate, crv2) / tminf->df;
      tminf->vol2_cum = fx_cum_vol_func(stp->time, SRT_YES, ts);
      tminf->ajust_init_fwd_val =
          tminf->init_fwd_val * exp(-0.5 * tminf->vol2_cum);
    } else
        /* For FX_STOCH_RATES init_fwd_val = spot */
        if (mdl_type == FX_STOCH_RATES) {
      tminf->init_fwd_val = spot;
      tminf->ajust_init_fwd_val =
          tminf->init_fwd_val * exp(-0.5 * tminf->vol2_cum);
    }

    stp = stp->next;
  }

  /* Sets the DETERMINSITIC drift (from todays forward values) */
  stp = top;
  while (stp->next) {
    tminf = stp->tminf[und_index];
    nxttminf = stp->next->tminf[und_index];
    tminf->drift =
        log(nxttminf->init_fwd_val / tminf->init_fwd_val) / stp->delta_t;
    stp = stp->next;
  }

  return NULL;
}

/* ======================================================================== */

/* ======================================================================

        Deterministic Model  specific functions

   ====================================================================== */

Err srt_f_basicinistp(SrtStpPtr stp, SrtUndPtr und, int und_index) {
  SrtCurvePtr yldcrv;
  SrtStpPtr top;
  Err err;
  SrtBasicTmInf *tminf;
  Ddate today;
  String und_name, yc_name;

  top = stp = gototop(stp);

  und_name = get_underlying_name(und);
  yc_name = get_ycname_from_irund(und);

  yldcrv = lookup_curve(yc_name);
  today = get_clcndate_from_yldcrv(yldcrv);

  /* allocate space */
  err = srtstpalloc(top, sizeof(SrtBasicTmInf), 0, und_index);
  if (err)
    return err;

  /* populate srttimestps */

  while (stp->next) {
    tminf = stp->tminf[und_index];
    tminf->df = swp_f_df(today, stp->ddate, yc_name);
    stp = stp->next;
  }

  /* last stp */
  tminf = stp->tminf[und_index];
  tminf->df = swp_f_df(today, stp->ddate, yc_name);

  und = lookup_und(und_name);
  return err;

} /* END Err srt_f_basicinistp(...) */

/* ========================================================================= */

/* -------------------------------------------------------------------------
   When Dealing with Multi-Currencies        , an extra initialisation is
   required: the QUANTO ADJUSTMENT in the drift The Name of the FX underlying
   should be FOR/DOM The Qunato Adjustement is stored in the TmInf attached to
   the SrtStp
   ------------------------------------------------------------------------- */

Err srt_f_steps_set_quanto_adjustment(SrtStpPtr stp, SrtUndInfo *und_info,
                                      int und_index, SrtUndPtr und,
                                      SrtUndPtr numeraire_und) {
  TermStruct *ts;    /* Term structure for for_und*/
  TermStruct *fx_ts; /* Term structure for FX_und*/
  TermStruct
      *dom_ts; /* Term structure for dom_und (need it in Jumping Numeraire)*/
  Err err;
  double correlation;
  SrtUnderlyingType und_type;
  char *und_name;
  double time;
  double fx_vol;
  char FX_name[32];
  double correlation_sign;
  double M_dt, Ofd_dt, Rfd_dt;
  double for_G, for_H, for_G_dt, for_S_dt, next_for_G, next_for_H;
  SrtCorrLstPtr cls;
  String und_ccy;
  String numeraire_ccy;
  String dom_und_name;
  SrtMdlType eModelType;
  SrtMdlDim eModelDim;
  SrtUndPtr fx_und;
  SrtUndPtr dom_und;

  /* Gets the Currency of this underlying (not necessarily different from
   * domestic...) */
  und_ccy = get_underlying_ccy(und);
  numeraire_ccy = get_underlying_ccy(numeraire_und);

  /* If the Currencies are the same        , the adjustment is zero */
  if (!strcmp(und_ccy, numeraire_ccy)) {
    return NULL;
  }

  /* For a normal quanto        , the adjustemnt is - rho(fx        ,for) *
   * fx_vol */
  correlation_sign = -1.0;

  /* Build the first possible FX name: "FOR/DOM" */
  sprintf(FX_name, "%s/%s", und_ccy, numeraire_ccy);

  /* Check the FX underlying is defined */
  fx_und = lookup_und(FX_name);
  if (fx_und == NULL) {
    /* The correlation will be the wrong way round */
    correlation_sign = 1.0;

    /* If not found: Build the second possible FX name: "DOM/FOR" */
    sprintf(FX_name, "%s/%s", numeraire_ccy, und_ccy);
    fx_und = lookup_und(FX_name);
    if (fx_und == NULL)
      return serror("Can't find the FX for %s/%s or even %s/%S ", und_ccy,
                    numeraire_ccy, numeraire_ccy, und_ccy);
  }

  /* Get the FX underlying Term Structure */
  err = get_underlying_ts(fx_und, &fx_ts);
  if (err) {
    return (err);
  }

  /* Get the Foreign underlying name (the one that needs a quanto adj) */
  und_name = get_underlying_name(und);

  /* Get the type of the (foreign) underlying */
  und_type = get_underlying_type(und);

  /* Get the foreign undelying Term Structure */
  err = get_underlying_ts(und, &ts);
  if (err) {
    return (err);
  }

  /* Get the Model type of the (foreign) underlying */
  get_underlying_mdltype(und, &eModelType);
  get_underlying_mdldim(und, &eModelDim);

  /* Gets the Global correlation list (not only the deal one: the FX could be
   * defined but not used) */
  cls = srt_f_GetTheCorrelationList();

  /* Extends the Undelying Term Structure if LGM */
  if (((eModelType == LGM) && (eModelDim == ONE_FAC)) &&
      (und_info->jumping == SRT_YES)) {
    /* Get the domestic underlying Term Structure */
    dom_und_name = get_domname_from_fxund(fx_und);
    dom_und = lookup_und(dom_und_name);
    if (!dom_und)
      return serror("Could not find %s underlying", dom_und_name);

    err = get_underlying_ts(dom_und, &dom_ts);
    if (err) {
      return (err);
    }
    /* Computes the M function for the LGm Term Structure */
    err = srt_f_extend_lgm_jumping_ts_for_quanto(ts, fx_ts, und_name, FX_name);
    if (err)
      return err;
  }

  /* Populate all the time steps of the dicretisation */
  stp = gototop(stp);
  while (stp->next) {
    time = stp->time;

    /* If LGM and Jumping: special quanto adjustment */
    if ((und_type == INTEREST_RATE_UND) && (eModelType == LGM) &&
        (eModelDim == ONE_FAC) && (und_info->jumping == SRT_YES)) {
      /* Computes the Quanto Adjustment Term : - M(Ti)  */

      G_H_func(stp->time, ts, &for_G, &for_H);
      G_H_func(stp->next->time, ts, &next_for_G, &next_for_H);

      M_dt = M_func(stp->next->time, ts) - M_func(stp->time, ts);
      for_S_dt = S_func(stp->next->time, ts) - S_func(stp->time, ts);

      for_G_dt = next_for_G - for_G;
      Ofd_dt = O_fd_func(stp->next->time, fx_ts) - O_fd_func(stp->time, fx_ts);
      Rfd_dt = R_fd_func(stp->next->time, fx_ts) - R_fd_func(stp->time, fx_ts);

      ((SrtIRMTmInf *)(stp->tminf[und_index]))->fxquanto_adjustment =
          correlation_sign * F_func(stp->next->time, ts) * M_dt;

      ((SrtIRMTmInf *)(stp->tminf[und_index]))->ffdquanto_adjustment =
          -correlation_sign * F_func(stp->next->time, dom_ts) *
          (Rfd_dt - Psi_func(stp->next->time, ts) * Ofd_dt);

      ((SrtIRMTmInf *)(stp->tminf[und_index]))->sfdquanto_adjustment =
          -correlation_sign * F_func(stp->next->time, ts) *
          (Psi_func(stp->next->time, ts) * for_G_dt - for_S_dt);

    } else {
      /* Get the correlation between the FX underlying and the foreign one (from
       * the Global Correlation Curve) */
      err = srt_f_get_corr_from_CorrList(cls, FX_name, und_name, time,
                                         &correlation);
      if (err)
        return (err);

      /* Get the FX volatility */
      fx_vol = find_fx_sig(time, fx_ts);

      /* The Quanto Adjustment itself */
      if ((und_type == FOREX_UND) || (und_type == EQUITY_UND)) {
        ((SrtLogTmInf *)(stp->tminf[und_index]))->quanto_adjustment =
            correlation_sign * correlation * fx_vol;
      } else {
        ((SrtIRMTmInf *)(stp->tminf[und_index]))->quanto_adjustment =
            correlation_sign * correlation * fx_vol;
      }
    }

    /* Move to the next time step */
    stp = stp->next;
  }

  /* Return a success mesage */
  return (err);

} /* Err srt_f_steps_set_quanto_adjustment(...) */

/* ============================================================================
 */

/* ===========================================================================================

                                                        FX STOCHASTIC RATES WITH
JUMPING NUMERAIRE

==============================================================================================*/

Err srt_f_fx_initstp(SrtStpPtr stp, SrtUndInfo *und_info, SrtUndPtr und,
                     SrtUndPtr dom_und, SrtUndPtr for_und, int und_index,
                     SrtUndPtr numeraire_und)

{

  SrtStpPtr top;
  Err err = NULL;
  SrtFXTmInf *tminf, *nexttminf;
  TermStruct *ts, *dom_ts, *for_ts, *gen_ts;
  String szCrv1, szCrv2, und_name;
  String dom_und_name;
  String for_und_name;
  String *used_und_names = NULL;
  double spot;
  int i;
  SrtCurvePtr yldcrv;
  SrtMdlType dom_mdl_type, for_mdl_type;
  SrtUnderlyingType sObjType1;
  Ddate today;
  double dom_Phi, dom_Psi, dom_F, dom_G, dom_H, dom_lambda_dt, next_dom_G,
      next_dom_H, dom_G_dt, next_dom_Psi, next_dom_F;
  double for_Psi, for_F, for_G, for_H, next_for_G, next_for_H, for_G_dt,
      next_for_Psi, next_for_F, for_lambda_dt;
  double dom_S_dt, dom_K_dt, dom_L_dt, dom_O_dt, dom_Psi_dt;
  double for_S_dt, for_K_dt, for_L_dt, for_O_dt, for_Psi_dt;
  double Mfd_dt, Ofd_dt, Pfd_dt, Qfd_dt, Rfd_dt;
  double Vdx_dt, Wdx_dt;
  double Vfx_dt, Wfx_dt;
  double mean_int_dom_sr_dt, mean_int_for_sr_dt, mean_int_fx_dt, mean_ln_fx;
  double var_int_dom_sr, var_int_for_sr, var_int_fx, var_ln_fx;
  double std_ln_fx;
  double cov_int_dom_fx, cov_int_for_fx, cov_dom_for_sr;
  double **new_corr_matrix;
  double cov_int_dom_for;
  double var_dom_sr, var_for_sr, stdev_dom_sr, stdev_for_sr, cov_ln_fx_dom_sr,
      cov_ln_fx_for_sr;
  String **ppszUndNames;
  long num_time_points;
  long step_index;
  double *pdDates;
  SrtCorrLstPtr corrlist;
  SrtCorrLstPtr cls = NULL;

  /* Gets the underlying name	and the spot value */

  und_name = get_underlying_name(und);
  spot = get_spot_from_fxund(und);

  /* Gets the term structure */
  err = get_underlying_ts(und, &ts);
  if (err)
    return err;

  /* Gets the domestic underlying name and pointer to */
  dom_und_name = get_domname_from_fxund(und);
  dom_und = lookup_und(dom_und_name);
  if (!dom_und)
    return serror("Could not find %s underlying", dom_und_name);

  /* Get the discount curve from the domestic underlying */
  szCrv1 = get_discname_from_underlying(dom_und);

  /* Get the domestic mdl type */
  err = get_underlying_mdltype(dom_und, &dom_mdl_type);
  if (err)
    return (err);

  /* Gets the foreign underlying name and pointer to  */

  for_und_name = get_forname_from_fxund(und);
  for_und = lookup_und(for_und_name);
  if (!for_und)
    return serror("Could not find %s underlying", for_und_name);

  /* Get the discount curve for the foreign underlying */
  szCrv2 = get_discname_from_underlying(for_und);

  /* Get the foreign mdl type */
  err = get_underlying_mdltype(for_und, &for_mdl_type);
  if (err)
    return (err);

  if ((dom_mdl_type != LGM) && (for_mdl_type != LGM))
    return serror("dom & for  mdl type should be LGM");

  /* Get the discount curve  */

  yldcrv = lookup_curve(szCrv1);
  sObjType1 = get_curve_type(yldcrv);
  today = get_today_from_curve(yldcrv);

  /* Go to the first time step */
  top = stp = gototop(stp);

  /* Allocate space for the Time Info */
  err = srtstpalloc(top, sizeof(SrtFXTmInf), 0, und_index);
  if (err)
    return err;

  /*----------------------	Populate srttimestps
   * -----------------------------*/

  stp = top;
  /* Set the volatility */
  while (stp->next) {
    tminf = stp->tminf[und_index];
    tminf->int_sigx_dt = fx_cum_vol_func(stp->next->time, SRT_NO, ts) -
                         fx_cum_vol_func(stp->time, SRT_NO, ts);
    tminf->int_sigx2_dt = fx_cum_vol_func(stp->next->time, SRT_YES, ts) -
                          fx_cum_vol_func(stp->time, SRT_YES, ts);
    tminf->sqrt_int_sigx2_dt = sqrt(tminf->int_sigx2_dt);

    stp = stp->next;
  }

  stp = top;
  while (stp) {
    tminf = stp->tminf[und_index];

    /* Sets the discount factors from the first underlying */
    tminf->df = swp_f_df(today, stp->ddate, szCrv1);
    if (tminf->df == SRT_DF_ERROR)
      return serror("Could not compute df for %s curve", szCrv1);

    /* Sets the forward value of the exchange rate */
    tminf->StochRatesVal[0].df = tminf->df;
    tminf->StochRatesVal[1].df = swp_f_df(today, stp->ddate, szCrv2);
    tminf->init_spot = spot;

    stp = stp->next;
  }

  /* Get 	the values of the discretised functions   */
  stp = top;
  while (stp) {
    tminf = (SrtFXTmInf *)stp->tminf[und_index];

    for (i = 0; i < 2; i++) {
      if (i == 0)
        err = get_underlying_ts(dom_und, &gen_ts);
      else if (i == 1)
        err = get_underlying_ts(for_und, &gen_ts);

      tminf->StochRatesVal[i].F = F_func(stp->time, gen_ts);
      G_H_func(stp->time, gen_ts, &tminf->StochRatesVal[i].G,
               &tminf->StochRatesVal[i].H);
      tminf->StochRatesVal[i].I = I_func(stp->time, gen_ts);
      tminf->StochRatesVal[i].K = K_func(stp->time, gen_ts);
      tminf->StochRatesVal[i].L = L_func(stp->time, gen_ts);
      tminf->StochRatesVal[i].O = O_func(stp->time, gen_ts);
      tminf->StochRatesVal[i].Q = Q_func(stp->time, gen_ts);
      tminf->StochRatesVal[i].Psi = Psi_func(stp->time, gen_ts);
      tminf->StochRatesVal[i].Phi = Phi_func(stp->time, gen_ts);
    }

    stp = stp->next;
  }

  /* Allocate space for a new correlation matrix with values at each time step
   * (now we look at the Fx Forward) */
  num_time_points = 0;
  num_time_points = create_index(top) + 1;
  new_corr_matrix = dmatrix(0, 2, 0, num_time_points - 1);
  pdDates = dvector(0, num_time_points - 1);

  stp = top;
  step_index = 0;
  while (stp->next) {
    tminf = (SrtFXTmInf *)stp->tminf[und_index];
    nexttminf = (SrtFXTmInf *)stp->next->tminf[und_index];
    err = get_underlying_ts(dom_und, &dom_ts);
    if (err)
      return err;

    err = get_underlying_ts(for_und, &for_ts);
    if (err)
      return err;

    /* computation of the integral between step and next step */
    Vdx_dt = V_dx_func(stp->next->time, ts) - V_dx_func(stp->time, ts);
    Wdx_dt = W_dx_func(stp->next->time, ts) - W_dx_func(stp->time, ts);
    Vfx_dt = V_fx_func(stp->next->time, ts) - V_fx_func(stp->time, ts);
    Wfx_dt = W_fx_func(stp->next->time, ts) - W_fx_func(stp->time, ts);

    Qfd_dt = Q_fd_func(stp->next->time, ts) - Q_fd_func(stp->time, ts);
    Ofd_dt = O_fd_func(stp->next->time, ts) - O_fd_func(stp->time, ts);
    Pfd_dt = P_fd_func(stp->next->time, ts) - P_fd_func(stp->time, ts);
    Rfd_dt = R_fd_func(stp->next->time, ts) - R_fd_func(stp->time, ts);
    Mfd_dt = M_fd_func(stp->next->time, ts) - M_fd_func(stp->time, ts);

    dom_K_dt = K_func(stp->next->time, dom_ts) - K_func(stp->time, dom_ts);
    dom_Psi_dt =
        Psi_func(stp->next->time, dom_ts) - Psi_func(stp->time, dom_ts);
    dom_L_dt = L_func(stp->next->time, dom_ts) - L_func(stp->time, dom_ts);
    dom_O_dt = O_func(stp->next->time, dom_ts) - O_func(stp->time, dom_ts);
    dom_lambda_dt = Lambda_func(stp->time, stp->next->time, dom_ts);

    G_H_func(stp->time, dom_ts, &dom_G, &dom_H);
    G_H_func(stp->next->time, dom_ts, &next_dom_G, &next_dom_H);
    dom_G_dt = next_dom_G - dom_G;

    for_L_dt = L_func(stp->next->time, for_ts) - L_func(stp->time, for_ts);
    for_O_dt = O_func(stp->next->time, for_ts) - O_func(stp->time, for_ts);
    dom_S_dt = S_func(stp->next->time, dom_ts) - S_func(stp->time, dom_ts);
    for_S_dt = S_func(stp->next->time, for_ts) - S_func(stp->time, for_ts);
    for_K_dt = K_func(stp->next->time, for_ts) - K_func(stp->time, for_ts);
    for_Psi_dt =
        Psi_func(stp->next->time, for_ts) - Psi_func(stp->time, for_ts);
    for_lambda_dt = Lambda_func(stp->time, stp->next->time, for_ts);

    G_H_func(stp->time, for_ts, &for_G, &for_H);
    G_H_func(stp->next->time, for_ts, &next_for_G, &next_for_H);
    for_G_dt = next_for_G - for_G;

    /* computation of the step information */
    dom_F = F_func(stp->time, dom_ts);
    dom_Psi = Psi_func(stp->time, dom_ts);
    dom_Phi = Phi_func(stp->time, dom_ts);
    for_F = F_func(stp->time, for_ts);
    for_Psi = Psi_func(stp->time, for_ts);

    /* computation of the next step information */
    next_dom_Psi = Psi_func(stp->next->time, dom_ts);
    next_dom_F = F_func(stp->next->time, dom_ts);
    next_for_Psi = Psi_func(stp->next->time, for_ts);
    next_for_F = F_func(stp->next->time, for_ts);

    mean_int_dom_sr_dt =
        -log(nexttminf->StochRatesVal[0].df) + log(tminf->StochRatesVal[0].df) -
        next_dom_Psi * dom_L_dt + dom_O_dt + pow(dom_lambda_dt, 2) * dom_Phi;

    cov_int_dom_fx = next_dom_Psi * Vdx_dt - Wdx_dt;
    cov_int_for_fx = next_for_Psi * Vfx_dt - Wfx_dt;
    cov_int_dom_for = next_dom_Psi * next_for_Psi * Ofd_dt -
                      next_for_Psi * Rfd_dt - next_dom_Psi * Pfd_dt + Mfd_dt;

    mean_int_for_sr_dt =
        -log(nexttminf->StochRatesVal[1].df) + log(tminf->StochRatesVal[1].df) +
        next_for_Psi * for_L_dt - for_O_dt - cov_int_for_fx - cov_int_dom_for;

    mean_int_fx_dt = Wdx_dt - next_dom_Psi * Vdx_dt;

    mean_ln_fx = mean_int_dom_sr_dt - mean_int_for_sr_dt + mean_int_fx_dt;

    tminf->mean_ln_fx = mean_ln_fx;

    /* compute the variance */
    var_int_dom_sr =
        next_dom_Psi * next_dom_Psi * dom_G_dt -
        2 * next_dom_Psi *
            (next_dom_Psi * next_dom_G - dom_Psi * dom_G - dom_L_dt) +
        next_dom_G * next_dom_Psi * next_dom_Psi - dom_G * dom_Psi * dom_Psi -
        2 * dom_O_dt;

    var_int_for_sr =
        next_for_Psi * next_for_Psi * for_G_dt -
        2 * next_for_Psi *
            (next_for_Psi * next_for_G - for_Psi * for_G - for_L_dt) +
        next_for_G * next_for_Psi * next_for_Psi - for_G * for_Psi * for_Psi -
        2 * for_O_dt;

    var_int_fx = tminf->int_sigx2_dt;

    var_ln_fx = var_int_dom_sr + var_int_for_sr + var_int_fx -
                2 * cov_int_dom_for + 2 * cov_int_dom_fx - 2 * cov_int_for_fx;
    std_ln_fx = sqrt(var_ln_fx);

    tminf->var_ln_fx = var_ln_fx;

    tminf->dom_lambda_dt = dom_lambda_dt;
    tminf->for_lambda_dt = for_lambda_dt;

    /* TEST VARIABLES: TO SUPRESS LATER */
    tminf->mean_int_dom_sr_dt = mean_int_dom_sr_dt;
    tminf->mean_int_for_sr_dt = mean_int_for_sr_dt;
    tminf->var_int_dom_sr = var_int_dom_sr;
    tminf->var_int_for_sr = var_int_for_sr;
    tminf->cov_int_dom_for = cov_int_dom_for;
    tminf->cov_int_dom_fx = cov_int_dom_fx;
    tminf->cov_int_for_fx = cov_int_for_fx;

    /* computation of the new correlation : ( 0: FX/DOM | 1 : FX/FOR | 2 :
     * DOM/FOR ) */

    var_dom_sr = next_dom_F * next_dom_F * dom_G_dt;
    var_for_sr = next_for_F * next_for_F * for_G_dt;
    stdev_dom_sr = sqrt(var_dom_sr);
    stdev_for_sr = sqrt(var_for_sr);

    /* Covariance between rd and Z */
    cov_ln_fx_dom_sr = next_dom_F * (next_dom_Psi * dom_G_dt - dom_S_dt) +
                       next_dom_F * (Pfd_dt - next_for_Psi * Ofd_dt) +
                       next_dom_F * Vdx_dt;

    /* Covariance between rf and Z */
    cov_ln_fx_for_sr = -next_for_F * (next_for_Psi * for_G_dt - for_S_dt) -
                       next_for_F * (Rfd_dt - next_dom_Psi * Ofd_dt) +
                       next_for_F * Vfx_dt;

    cov_dom_for_sr = next_dom_F * next_for_F * Ofd_dt;

    new_corr_matrix[0][step_index] =
        cov_ln_fx_dom_sr / (stdev_dom_sr * std_ln_fx);

    new_corr_matrix[1][step_index] =
        cov_ln_fx_for_sr / (stdev_for_sr * std_ln_fx);

    new_corr_matrix[2][step_index] =
        cov_dom_for_sr / (stdev_for_sr * stdev_dom_sr);

    if (err) {
      free_dmatrix(new_corr_matrix, 0, 2, 0, num_time_points - 1);
      free_dvector(pdDates, 0, num_time_points - 1);
      return err;
    }

    pdDates[step_index] = stp->ddate;

    /* Moves on to the next step */
    stp = stp->next;
    step_index++;
  }

  /* Sets the Underlying Names in a [3][2] array */
  ppszUndNames = smatrix_size(0, 2, 0, 1, 128);
  strncpy(ppszUndNames[0][0], und_name, strlen(und_name));
  strncpy(ppszUndNames[0][1], dom_und_name, strlen(dom_und_name));
  strncpy(ppszUndNames[1][0], und_name, strlen(und_name));
  strncpy(ppszUndNames[1][1], for_und_name, strlen(for_und_name));
  strncpy(ppszUndNames[2][0], dom_und_name, strlen(dom_und_name));
  strncpy(ppszUndNames[2][1], for_und_name, strlen(for_und_name));

  /* Replaces the current correlation matrix by the new one (with coeffcients
   * done inside) in UndInfo */
  err = srt_f_init_Corr_TermStruct(step_index, 3, new_corr_matrix, pdDates,
                                   ppszUndNames, &und_info->corr_ts);

  /* Just stacks the underlying names in one single vector */
  used_und_names = svector(0, und_info->no_of_underlyings - 1);
  for (i = 0; i < und_info->no_of_underlyings; i++)
    used_und_names[i] = und_info->und_data[i].und_name;

  /* Gets the Global Correlation list */
  corrlist = srt_f_GetTheCorrelationList();

  /* Creates and fills the local SrtCorrLst with relevent underlying correlation
   * matrix*/
  err = srt_f_make_deal_corrlist(used_und_names, und_info->no_of_underlyings,
                                 "Local deal correlation list", corrlist, &cls);

  if (err) {
    free_svector(used_und_names, 0, und_info->no_of_underlyings - 1);
    return err;
  }

  /* Attach the local SrtCorrLst to the und_info */
  und_info->corr_ts = cls;

  err = srt_f_attach_correl_to_stp(top, und_info->corr_ts);
  if (err)
    return err;

  /* Frees the Memory allocated */
  free_svector(used_und_names, 0, und_info->no_of_underlyings - 1);
  free_dmatrix(new_corr_matrix, 0, 2, 0, num_time_points - 1);
  free_smatrix_size(ppszUndNames, 0, 2, 0, 1, 128);
  free_dvector(pdDates, 0, num_time_points - 1);

  return err;

} /* END  Err srt_f_fx_initstp() */

/*======================END OF FX STOCHASTIC RATES WITH JUMPING NUMERAIRE ====*/
