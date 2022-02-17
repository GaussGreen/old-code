/*	==========================================================================
        FILE_NAME:	Fx3FBetaDLMCalculations.c

        PURPOSE:	All the calculations in the Fx3DBetaDLM model

        DATE:		10/15/03

        AUTHOR:		L.C.
        ==========================================================================
 */

#include "Fx3FBetaDLMCalculations.h"
#include "Fx3FBetaDLMCalibration.h"
#include "Fx3FBetaDLMUtil.h"
#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "Math.h"
#include "opfnctns.h"

Err FxBetaDLM_FxOptionUnd(char *und, long opt_settlmt_date, long fwd_start_date,
                          long fixing_date, int nb_strike, double *strikes,
                          SrtCallPutType callput,
                          FxBetaDLM_OptNumerParams *NumParams,
                          double *premium) {
  Err err = NULL;
  FxBetaDLM_model *model = NULL;
  FxBetaDLM_Hermite *hermite = NULL;

  model = calloc(1, sizeof(FxBetaDLM_model));
  hermite = calloc(1, sizeof(FxBetaDLM_Hermite));

  if (!model || !hermite) {
    err = "Memory allocation faillure in FxBetaDLM_FxOptionUnd";
    goto FREE_RETURN;
  }

  err = FxBetaDLM_Get_Model(und, model);

  if (err)
    goto FREE_RETURN;

  err = initialise_FxBetaDLM_Hermite(NumParams, hermite);

  if (err)
    goto FREE_RETURN;

  err = FxBetaDLM_FxOptionTS_struct(opt_settlmt_date, fwd_start_date,
                                    fixing_date, nb_strike, strikes, callput,
                                    model, NumParams, hermite, premium);

FREE_RETURN:

  if (model) {
    free_FxBetaDLM_model(model);
    free(model);
  }

  if (hermite) {
    free_FxBetaDLM_Hermite(hermite);
    free(hermite);
  }

  return err;
}

Err FxBetaDLM_Allocate_FxOptInst(int iIndex, int iNbStrike,
                                 FxBetaDLM_FxOptInst *Inst) {
  Err err = NULL;

  Inst->iInstIndex = iIndex;
  Inst->dStrike = calloc(iNbStrike, sizeof(double));
  Inst->dLogStrike = calloc(iNbStrike, sizeof(double));
  Inst->dPrice = calloc(iNbStrike, sizeof(double));
  Inst->dLogVol = calloc(iNbStrike, sizeof(double));
  Inst->dVega = calloc(iNbStrike, sizeof(double));

  if (!Inst->dStrike || !Inst->dLogStrike || !Inst->dPrice || !Inst->dLogVol ||
      !Inst->dVega) {
    err = "Memory allocation faillure in FxBetaDLM_Allocate_FxOptInst";
    goto FREE_RETURN;
  }

FREE_RETURN:

  if (err) {
    FxBetaDLM_Free_FxOptInst(Inst);
  }

  return err;
}

void FxBetaDLM_Free_FxOptInst(FxBetaDLM_FxOptInst *Inst) {
  if (Inst) {
    if (Inst->dStrike)
      free(Inst->dStrike);
    Inst->dStrike = NULL;

    if (Inst->dLogStrike)
      free(Inst->dLogStrike);
    Inst->dLogStrike = NULL;

    if (Inst->dPrice)
      free(Inst->dPrice);
    Inst->dPrice = NULL;

    if (Inst->dLogVol)
      free(Inst->dLogVol);
    Inst->dLogVol = NULL;

    if (Inst->dVega)
      free(Inst->dVega);
    Inst->dVega = NULL;
  }
}

Err FxBetaDLM_Setup_FxOptInst(long opt_settlmt_date, long opt_pay_date,
                              long fwd_start_date, long fixing_date,
                              int nb_strike, double *strikes, double *bs_vols,
                              SrtCallPutType callput, FxBetaDLM_model *model,
                              FxBetaDLM_FxOptInst *Inst) {
  Err err = NULL;
  int i;
  double std_fx, temp;

  Inst->lExeDate = fixing_date;
  Inst->lFwdStartDate = fwd_start_date;
  Inst->lPayDate = opt_pay_date;
  Inst->lSettlementDate = opt_settlmt_date;

  Inst->dExeTime = (fixing_date - model->lToday) * YEARS_IN_DAY;
  Inst->dFwdStartTime = (fwd_start_date - model->lToday) * YEARS_IN_DAY;
  Inst->dPayTime = (opt_pay_date - model->lToday) * YEARS_IN_DAY;
  Inst->dSettlementTime = (opt_settlmt_date - model->lToday) * YEARS_IN_DAY;

  if (Inst->dFwdStartTime > Inst->dExeTime) {
    err = "fixing date must be before fwd_start_date in "
          "FxBetaDLM_Setup_FxOptInst";
    goto FREE_RETURN;
  }

  Inst->dVolTime = Inst->dExeTime - Inst->dFwdStartTime;
  Inst->dSqVolTime = sqrt(Inst->dVolTime);

  if (Inst->dSettlementTime < 1.0E-08) {
    err = "settlmt of the option must be after today";
    goto FREE_RETURN;
  }

  /* Get the informations */
  Inst->dDfDom = swp_f_df(model->lToday, opt_settlmt_date, model->cYcDom);
  Inst->dDfFor = swp_f_df(model->lToday, opt_settlmt_date, model->cYcFor);
  Inst->dDfPayDom = swp_f_df(model->lToday, opt_pay_date, model->cYcDom);
  Inst->dFwdFx = model->dCashFx * Inst->dDfFor / Inst->dDfDom;
  Inst->dLnFwdFx = log(Inst->dFwdFx);

  Inst->sCallPutType = callput;

  Inst->iNbStrike = nb_strike;

  if (fabs(strikes[0]) < 1.0E-08) {
    strikes[0] = Inst->dFwdFx;
  }

  if (bs_vols) {
    std_fx = bs_vols[0] * sqrt(Inst->dExeTime);
  }

  if (nb_strike > 2 && strikes[2] < 0.0 && bs_vols) {
    for (i = 1; i < nb_strike; i++) {
      strikes[i] = Inst->dFwdFx * exp(strikes[i] * std_fx);
      err = srt_f_optsarbvol(Inst->dFwdFx, strikes[i], Inst->dExeTime,
                             bs_vols[0], 0.0, bs_vols[i], 0.0, SRT_LOGNORMAL,
                             SRT_LOGNORMAL, &temp);

      bs_vols[i] = temp;

      if (err)
        goto FREE_RETURN;
    }
  }

  for (i = 0; i < nb_strike; i++) {
    Inst->dStrike[i] = strikes[i];
    Inst->dLogStrike[i] = log(strikes[i]);

    if (bs_vols) {
      Inst->dLogVol[i] = bs_vols[i];
      Inst->dPrice[i] = srt_f_optblksch(
          Inst->dFwdFx, Inst->dStrike[i], Inst->dLogVol[i], Inst->dExeTime,
          Inst->dDfDom, Inst->sCallPutType, PREMIUM);

      Inst->dVega[i] = srt_f_optblksch(Inst->dFwdFx, Inst->dStrike[i],
                                       Inst->dLogVol[i], Inst->dExeTime,
                                       Inst->dDfDom, Inst->sCallPutType, VEGA);
    }
  }

FREE_RETURN:

  return err;
}

Err FxBetaDLM_Price_FxOptInst(FxBetaDLM_FxOptInst *Inst, FxBetaDLM_model *Model,
                              FxBetaDLM_InstPrecalc *Precalc,
                              FxBetaDLM_Hermite *Hermite,
                              FxBetaDLM_OptNumerParams *cNumParams,
                              double *Premium) {
  Err err = NULL;
  int i, j;
  double fwd_fx, Val1, Val2, Val;

  /* Initialisation */
  if (Inst->iNbStrike == 0) {
    *Premium = 0.0;
    return err;
  }

  for (i = 0; i < Inst->iNbStrike; i++) {
    Premium[i] = 0.0;
  }

  /* Pricing With a first conditioning on the quadratic variable*/
  if (cNumParams->iMethod == 1) {
    for (i = 1; i <= Hermite->iNbPoints; i++) {
      fwd_fx = exp(Precalc->dConstCoef +
                   (Precalc->dLinCoef + Precalc->dQuadCoef * Hermite->X[i]) *
                       Hermite->X[i]);

      for (j = 0; j < Inst->iNbStrike; j++) {
        Premium[j] +=
            Hermite->W[i] * srt_f_optblksch(fwd_fx, Inst->dStrike[j],
                                            Precalc->dStdRates, 1.0, 1.0,
                                            Inst->sCallPutType, PREMIUM);
      }
    }
  }
  /* Pricing With a first conditioning on the Linear varible*/
  else if (cNumParams->iMethod == 2 || cNumParams->iMethod == 3 ||
           cNumParams->iMethod == 4) {
    if (Precalc->iIsCPD) {
      for (i = 1; i <= Hermite->iNbPoints; i++) {
        for (j = 0; j < Inst->iNbStrike; j++) {
          err = quadratic_bs_CPD(Precalc, i, Inst->dStrike[j],
                                 Inst->dLogStrike[j], Inst->sCallPutType, &Val);
          if (err)
            return err;

          Premium[j] += Hermite->W[i] * Val;
        }
      }
    } else {
      for (i = 1; i <= Hermite->iNbPoints; i++) {
        for (j = 0; j < Inst->iNbStrike; j++) {
          if (fabs(Model->dAlpha) < TINY) {
            err = quadratic_bs(
                Precalc->dQuadConst + Hermite->X[i] * Precalc->dQuadStdY,
                Precalc->dQuadLin, Precalc->dQuadQuad,
                Precalc->dQuadMean * Hermite->X[i] * Precalc->dQuadStdY,
                Precalc->dQuadVar, Inst->dStrike[j], Inst->sCallPutType, &Val);
          } else {
            err = quadratic_bs(
                Precalc->dQuadConst + Hermite->X[i] * Precalc->dQuadStdY,
                Precalc->dQuadLin, Precalc->dQuadQuad,
                Precalc->dQuadMean * Hermite->X[i] * Precalc->dQuadStdY,
                Precalc->dQuadVar, Inst->dStrike[j], Inst->sCallPutType, &Val1);

            err = quadratic_bs(Precalc->dQuadConst_down +
                                   Hermite->X[i] * Precalc->dQuadStdY_down,
                               Precalc->dQuadLin_down, Precalc->dQuadQuad_down,
                               Precalc->dQuadMean_down * Hermite->X[i] *
                                   Precalc->dQuadStdY_down,
                               Precalc->dQuadVar_down, Inst->dStrike[j],
                               Inst->sCallPutType, &Val2);

            Val = 0.5 * (Val1 + Val2);
          }

          if (err)
            return err;

          Premium[j] += Hermite->W[i] * Val;
        }
      }
    }
  } else {
    err = "Parameter Method is invalid";
    return err;
  }

  /* Discounting */
  for (j = 0; j < Inst->iNbStrike; j++) {
    Premium[j] *= Inst->dDfPayDom;
  }

  return err;
}

Err FxBetaDLM_FxOptionTS_struct(long opt_settlmt_date, long fwd_start_date,
                                long fixing_date, int nb_strike,
                                double *strikes, SrtCallPutType callput,
                                FxBetaDLM_model *model,
                                FxBetaDLM_OptNumerParams *NumParams,
                                FxBetaDLM_Hermite *hermite, double *premium) {
  Err err = NULL;
  FxBetaDLM_FxOptInst *Inst = NULL;
  FxBetaDLM_InstPrecalc *InstConst = NULL;
  double *SigFx;
  double *premium1;
  double *premium2;

  premium1 = calloc(nb_strike, sizeof(double));
  premium2 = calloc(nb_strike, sizeof(double));
  SigFx = calloc(model->iNbPWTime, sizeof(double));
  Inst = calloc(1, sizeof(FxBetaDLM_FxOptInst));
  InstConst = calloc(1, sizeof(FxBetaDLM_InstPrecalc));

  if (!Inst || !InstConst) {
    err = "Memory allocation faillure in FxBetaDLM_FxOptionTS_struct";
    goto FREE_RETURN;
  }

  /* Setup the instrument */
  err = FxBetaDLM_Allocate_FxOptInst(0, nb_strike, Inst);

  if (err)
    goto FREE_RETURN;

  err = FxBetaDLM_Setup_FxOptInst(opt_settlmt_date, opt_settlmt_date,
                                  fwd_start_date, fixing_date, nb_strike,
                                  strikes, NULL, callput, model, Inst);

  if (err)
    goto FREE_RETURN;

  err = FxBetaDLM_Calculate_AllConst(Inst, model, NumParams, InstConst);

  if (err)
    goto FREE_RETURN;

  err = FxBetaDLM_Update_InstPrecalc_FromMoment(Inst, model, NumParams,
                                                InstConst);

  if (err)
    goto FREE_RETURN;

  err = FxBetaDLM_Price_FxOptInst(Inst, model, InstConst, hermite, NumParams,
                                  premium);
  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (Inst) {
    FxBetaDLM_Free_FxOptInst(Inst);
    free(Inst);
  }

  if (InstConst)
    free(InstConst);
  if (SigFx)
    free(SigFx);
  if (premium1)
    free(premium1);
  if (premium2)
    free(premium2);
  return err;
}

Fx3DBetaDLM_ExpectAndVarGrfn(long nstp, double *time, FxBetaDLM_model *model,

                             double *dom_fwd, double *dom_var_std,
                             double *dom_phi, double *dom_beta,
                             double *for_fwd_const, double *for_fwd_lin,
                             double *for_var_std, double *for_phi,
                             double *for_beta, double *fx_fwd,
                             double *fx_var_std, double *ffx_var,

                             int mc_tree, /*0: MC  , 1: Tree */
                             double *dom_for_cov_corr, double *dom_fx_cov_corr,
                             double *for_fx_cov_corr,

                             double min_time) {
  int i, j, k;
  double t1, t2, dt;
  double var_rates_dom, var_rates_for, var_ffx;
  double covar_dom_for, covar_dom_ffx, covar_for_ffx, covar_dom_for_adjust;
  double sig_dom, sig_for, sig_fx, sig_dom2, sig_for2, sig_fx2, sig_domfor,
      sig_domfx, sig_forfx;
  double dom_lam, for_lam, tstar;

  double start_mat, end_mat;
  int start_index, end_index;

  int nb_points;
  double min_dt, old_t, new_t;
  double temp, temp_dom, temp_for, temp_bracket;
  double expect_for, total_var_ffx, var_0_Ti;
  double old_integ, new_integ, old_integ2, new_integ2;
  double last_dom_phi, last_for_phi, last_var_ffx;

  dom_lam = model->dLambdaDom;
  for_lam = model->dLambdaFor;
  tstar = model->dTStar;

  start_mat = time[0];
  start_index = Get_Index(start_mat, model->dPWTime, model->iNbPWTime);

  var_rates_dom = 0.0;
  var_rates_for = 0.0;
  covar_dom_for = 0.0;
  var_ffx = 0.0;
  covar_dom_ffx = 0.0;
  covar_for_ffx = 0.0;
  expect_for = 0.0;
  old_t = 0.0;
  old_integ = 0.0;
  old_integ2 = 0.0;
  total_var_ffx = 0.0;
  covar_dom_for_adjust = 0.0;

  last_dom_phi = 0.0;
  last_for_phi = 0.0;
  last_var_ffx = 0.0;

  for (i = 1; i < nstp; i++) {
    end_mat = time[i];
    end_index = Get_Index(end_mat, model->dPWTime, model->iNbPWTime);
    var_0_Ti = total_var_ffx;

    if (mc_tree == 0) {
      /* initialisation for the next period */

      if (i > 1) {
        last_dom_phi = dom_phi[i - 1];
        last_for_phi = for_phi[i - 1];
        last_var_ffx = ffx_var[i - 1];
      }

      var_rates_dom = 0.0;
      var_rates_for = 0.0;
      covar_dom_for = 0.0;
      var_ffx = 0.0;
      covar_dom_ffx = 0.0;
      covar_for_ffx = 0.0;
      expect_for = 0.0;
      covar_dom_for_adjust = 0.0;

      old_integ2 = 0.0;
    }

    for (j = start_index; j < end_index + 1; j++) {
      if (j > start_index) {
        t1 = model->dPWTime[j - 1];
      } else {
        /* First part */
        t1 = start_mat;
      }

      if (j == end_index || start_index == end_index) {
        /* Last part */
        t2 = end_mat;
      } else {
        t2 = model->dPWTime[j];
      }

      dt = t2 - t1;

      sig_dom = model->dSigmaDom[j];
      sig_for = model->dSigmaFor[j];
      sig_fx = model->dSigmaFx[j];
      sig_dom2 = sig_dom * sig_dom;
      sig_for2 = sig_for * sig_for;
      sig_fx2 = sig_fx * sig_fx;
      sig_domfor = model->dCorrDomFor[j] * sig_dom * sig_for;
      sig_domfx = model->dCorrDomFx[j] * sig_dom * sig_fx;
      sig_forfx = model->dCorrForFx[j] * sig_for * sig_fx;

      /* Expectation and Variance approximation for foreign */
      nb_points = max((int)(dt / min_time + 0.5), 1);
      min_dt = dt / nb_points;

      for (k = 0; k < nb_points; k++) {
        new_t = old_t + min_dt;

        /* update var_ffx */
        total_var_ffx += sig_fx2 * min_dt;
        total_var_ffx +=
            -2.0 * sig_forfx * Etha_Func(for_lam, tstar, old_t, new_t) +
            2.0 * sig_domfx * Etha_Func(dom_lam, tstar, old_t, new_t);
        total_var_ffx +=
            sig_for2 * Psi_Func(for_lam, for_lam, tstar, old_t, new_t) +
            sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, old_t, new_t) -
            2.0 * sig_domfor * Psi_Func(dom_lam, for_lam, tstar, old_t, new_t);

        covar_dom_ffx +=
            sig_domfx * Phi_Func(dom_lam, tstar, old_t, new_t) -
            sig_domfor * Gamma_Func(for_lam, dom_lam, tstar, old_t, new_t) +
            sig_dom2 * Gamma_Func(dom_lam, dom_lam, tstar, old_t, new_t);

        temp_dom = exp(-dom_lam * (model->dTStar - new_t));
        temp_for = exp(-for_lam * (model->dTStar - new_t));

        temp_bracket =
            temp_for * (sig_forfx - sig_for2 * (1.0 - temp_for) / for_lam +
                        sig_domfor * (1.0 - temp_dom) / dom_lam);

        new_integ = temp_bracket / (1.0 + 2.0 * model->dC0 * total_var_ffx);
        new_integ2 = covar_dom_ffx * new_integ;
        new_integ *= total_var_ffx;

        expect_for += (old_integ + new_integ) * min_dt / 2.0;
        covar_dom_for_adjust += (old_integ2 + new_integ2) * min_dt / 2.0;

        old_t = new_t;
        old_integ = new_integ;
        old_integ2 = new_integ2;
      }

      /* Update */
      var_rates_dom += sig_dom2 * Phi_Func(2.0 * dom_lam, tstar, t1, t2);
      var_rates_for += sig_for2 * Phi_Func(2.0 * for_lam, tstar, t1, t2);

      covar_dom_for += sig_domfor * Phi_Func(dom_lam + for_lam, tstar, t1, t2);

      /* Fx component */
      var_ffx += sig_fx2 * dt;

      /* Fx / Rates component */
      var_ffx += -2.0 * sig_forfx * Etha_Func(for_lam, tstar, t1, t2) +
                 2.0 * sig_domfx * Etha_Func(dom_lam, tstar, t1, t2);

      /* Rates / Rates component */
      var_ffx += sig_for2 * Psi_Func(for_lam, for_lam, tstar, t1, t2) +
                 sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, t1, t2) -
                 2.0 * sig_domfor * Psi_Func(dom_lam, for_lam, tstar, t1, t2);

      covar_for_ffx += sig_forfx * Phi_Func(for_lam, tstar, t1, t2) -
                       sig_for2 * Gamma_Func(for_lam, for_lam, tstar, t1, t2) +
                       sig_domfor * Gamma_Func(dom_lam, for_lam, tstar, t1, t2);
    }

    dom_beta[i] = (1.0 - exp(-dom_lam * (end_mat - tstar))) / dom_lam;
    dom_fwd[i] = 0.0;
    dom_var_std[i] = var_rates_dom;
    dom_phi[i] = var_rates_dom + last_dom_phi;

    for_beta[i] = (1.0 - exp(-for_lam * (end_mat - tstar))) / for_lam;

    temp = -covar_for_ffx + 2.0 * model->dC0 * expect_for;
    for_fwd_const[i] = model->dB0 * temp;
    for_var_std[i] = var_rates_for - 2.0 * model->dC0 *
                                         (1.0 + 2.0 * model->dC0 * var_0_Ti) *
                                         temp * temp;
    for_phi[i] = var_rates_for + last_for_phi;

    fx_fwd[i] = 0.0;
    fx_var_std[i] = var_ffx;
    ffx_var[i] = var_ffx + last_var_ffx;

    if (mc_tree == 0) {
      for_fwd_lin[i] = 2.0 * model->dC0 * temp;

      dom_var_std[i] = sqrt(dom_var_std[i]);
      for_var_std[i] = sqrt(for_var_std[i]);
      fx_var_std[i] = sqrt(fx_var_std[i]);

      dom_for_cov_corr[i] =
          covar_dom_for - 2.0 * model->dC0 * covar_dom_for_adjust;
      dom_fx_cov_corr[i] = covar_dom_ffx;
      for_fx_cov_corr[i] = -(1.0 + 2.0 * model->dC0 * var_0_Ti) * temp;
    } else {
      dom_for_cov_corr[i] =
          (covar_dom_for - 2.0 * model->dC0 * covar_dom_for_adjust) /
          sqrt(dom_var_std[i] * for_var_std[i]);
      dom_fx_cov_corr[i] = covar_dom_ffx / sqrt(dom_var_std[i] * fx_var_std[i]);
      for_fx_cov_corr[i] = -temp / sqrt(for_var_std[i] * fx_var_std[i]);
    }

    start_mat = end_mat;
    start_index = end_index;
  }
}

Fx3DBetaDLM_PrecalcGRFNTreeQBeta(
    long nstp, double *time, FxBetaDLM_model *model,

    double *dom_fwd, double *dom_var, double *dom_phi, double *dom_beta,
    double *for_fwd, double *for_var, double *for_phi, double *for_beta,
    double *fx_fwd, double *fx_var, double *ffx_var,

    double *dom_for_corr, double *dom_fx_corr, double *for_fx_corr,

    double min_time) {
  int i, j;
  double t1, t2, dt;
  double var_rates_dom, var_rates_for, var_ffx, var_lnZ;
  double covar_dom_for;
  double covar_dom_lnZ;
  double covar_for_lnZ;
  double sig_dom, sig_for, sig_fx, sig_dom2, sig_for2, sig_fx2, sig_domfor,
      sig_domfx, sig_forfx;
  double dom_lam, for_lam, tstar;

  double start_mat, end_mat;
  int start_index, end_index;

  double expect_for, expect_dom, expect_lnZ;

  dom_lam = model->dLambdaDom;
  for_lam = model->dLambdaFor;
  tstar = model->dTStar;

  start_mat = 0.0;
  start_index = Get_Index(start_mat, model->dPWTime, model->iNbPWTime);

  var_ffx = 0.0;
  var_lnZ = 0.0;
  var_rates_dom = 0.0;
  var_rates_for = 0.0;

  covar_dom_lnZ = 0.0;
  covar_for_lnZ = 0.0;
  covar_dom_for = 0.0;

  expect_for = 0.0;
  expect_dom = 0.0;
  expect_lnZ = 0.0;

  for (i = 1; i < nstp; i++) {
    end_mat = time[i];
    end_index = Get_Index(end_mat, model->dPWTime, model->iNbPWTime);

    var_ffx = 0.0;
    var_lnZ = 0.0;
    var_rates_dom = 0.0;
    var_rates_for = 0.0;

    covar_dom_lnZ = 0.0;
    covar_for_lnZ = 0.0;
    covar_dom_for = 0.0;

    expect_for = 0.0;
    expect_dom = 0.0;
    expect_lnZ = 0.0;

    for (j = start_index; j < end_index + 1; j++) {
      if (j > start_index) {
        t1 = model->dPWTime[j - 1];
      } else {
        /* First part */
        t1 = start_mat;
      }

      if (j == end_index || start_index == end_index) {
        /* Last part */
        t2 = end_mat;
      } else {
        t2 = model->dPWTime[j];
      }

      dt = t2 - t1;

      sig_dom = model->dSigmaDom[j];
      sig_for = model->dSigmaFor[j];
      sig_fx = model->dSigmaFx[j];
      sig_dom2 = sig_dom * sig_dom;
      sig_for2 = sig_for * sig_for;
      sig_fx2 = sig_fx * sig_fx;
      sig_domfor = model->dCorrDomFor[j] * sig_dom * sig_for;
      sig_domfx = model->dCorrDomFx[j] * sig_dom * sig_fx;
      sig_forfx = model->dCorrForFx[j] * sig_for * sig_fx;

      /* Update */
      var_rates_dom += sig_dom2 * Phi_Func(2.0 * dom_lam, end_mat, t1, t2);
      var_rates_for += sig_for2 * Phi_Func(2.0 * for_lam, end_mat, t1, t2);

      var_lnZ += sig_fx2 * dt;
      var_lnZ += -2.0 * sig_forfx * Etha_Func(for_lam, end_mat, t1, t2) +
                 2.0 * sig_domfx * Etha_Func(dom_lam, end_mat, t1, t2);
      var_lnZ += sig_for2 * Psi_Func(for_lam, for_lam, end_mat, t1, t2) +
                 sig_dom2 * Psi_Func(dom_lam, dom_lam, end_mat, t1, t2) -
                 2.0 * sig_domfor * Psi_Func(dom_lam, for_lam, end_mat, t1, t2);

      covar_dom_for +=
          sig_domfor * Phi_Func(dom_lam + for_lam, end_mat, t1, t2);

      covar_dom_lnZ += sig_domfx * Phi_Func(dom_lam, end_mat, t1, t2);
      covar_dom_lnZ -=
          sig_domfor * Gamma_Func(for_lam, dom_lam, end_mat, t1, t2);
      covar_dom_lnZ += sig_dom2 * Gamma_Func(dom_lam, dom_lam, end_mat, t1, t2);

      covar_for_lnZ += sig_forfx * Phi_Func(for_lam, end_mat, t1, t2);
      covar_for_lnZ -= sig_for2 * Gamma_Func(for_lam, for_lam, end_mat, t1, t2);
      covar_dom_lnZ +=
          sig_domfor * Gamma_Func(dom_lam, for_lam, end_mat, t1, t2);

      expect_dom += sig_dom2 * Gamma_Func(dom_lam, dom_lam, end_mat, t1, t2);
      expect_for += sig_for2 * Gamma_Func(dom_lam, dom_lam, end_mat, t1, t2);
      expect_for -= sig_forfx * Phi_Func(for_lam, end_mat, t1, t2);
      // expect_lnZ += sig_forfx * Etha_Func(for_lam  , end_mat  , t1  , t2);
      expect_lnZ += -sig_domfor * Psi_Func(dom_lam, for_lam, end_mat, t1, t2) +
                    sig_domfx * Etha_Func(dom_lam, end_mat, t1, t2) +
                    sig_dom2 * Psi_Func(dom_lam, dom_lam, end_mat, t1, t2);

      /* Fx component */
      var_ffx += sig_fx2 * dt;

      /* Fx / Rates component */
      var_ffx += -2.0 * sig_forfx * Etha_Func(for_lam, tstar, t1, t2) +
                 2.0 * sig_domfx * Etha_Func(dom_lam, tstar, t1, t2);

      /* Rates / Rates component */
      var_ffx += sig_for2 * Psi_Func(for_lam, for_lam, tstar, t1, t2) +
                 sig_dom2 * Psi_Func(dom_lam, dom_lam, tstar, t1, t2) -
                 2.0 * sig_domfor * Psi_Func(dom_lam, for_lam, tstar, t1, t2);
    }

    dom_fwd[i] = expect_dom;
    dom_var[i] = var_rates_dom;
    dom_phi[i] = var_rates_dom;
    dom_beta[i] = (1.0 - exp(-dom_lam * (tstar - end_mat))) / dom_lam;

    for_fwd[i] = expect_for;
    for_var[i] = var_rates_for;
    for_phi[i] = var_rates_for;
    for_beta[i] = (1.0 - exp(-for_lam * (tstar - end_mat))) / for_lam;

    fx_fwd[i] = -0.5 * var_lnZ + expect_lnZ;
    fx_var[i] = var_lnZ;

    ffx_var[i] = var_ffx;

    dom_for_corr[i] = covar_dom_for / sqrt(dom_var[i] * for_var[i]);
    dom_fx_corr[i] = covar_dom_lnZ / sqrt(dom_var[i] * fx_var[i]);
    for_fx_corr[i] = covar_for_lnZ / sqrt(for_var[i] * fx_var[i]);
  }
}

Err quadratic_bs_CPD(FxBetaDLM_InstPrecalc *Precalc, int i, double strike,
                     double LogStrike, SrtCallPutType callput,
                     double *Premium) {
  double root_1, root_2, swap_root, k1, k2, J1, J2, delta;
  double constant_coef, mean;
  Err err = NULL;

  constant_coef = Precalc->dQuadConst + Precalc->coupon_hermite[i];
  mean = Precalc->dQuadMean * Precalc->coupon_hermite[i];

  if (Precalc->dQuadVar == 0.0) {
    if (callput == SRT_CALL) {
      *Premium = exp(constant_coef + Precalc->dQuadLin * mean +
                     Precalc->dQuadQuad * mean * mean) -
                 strike;
      if (*Premium < 0) {
        *Premium = 0;
      }
    } else {
      *Premium = strike - exp(constant_coef + Precalc->dQuadLin * mean +
                              Precalc->dQuadQuad * mean * mean);
      if (*Premium < 0) {
        *Premium = 0;
      }
    }
  } else {
    if (Precalc->dQuadQuad >= 0.5 / Precalc->dQuadVar) {
      err = serror("quadratic_bs  The process has infinite Precalc->dQuadVar  "
                   ", C = %f > %f",
                   Precalc->dQuadQuad, 0.5 / Precalc->dQuadVar);
      goto FREE_RETURN;
    } else if (Precalc->dQuadQuad == 0) {
      root_1 = (LogStrike - constant_coef) / Precalc->dQuadLin;
      if (callput == SRT_CALL) {
        J1 = exp(constant_coef - mean * mean / (2.0 * Precalc->dQuadVar) +
                 (Precalc->dQuadLin * Precalc->dQuadVar + mean) *
                     (Precalc->dQuadLin * Precalc->dQuadVar + mean) /
                     (2.0 * Precalc->dQuadVar));
        J1 *= norm((mean + Precalc->dQuadLin * Precalc->dQuadVar - root_1) /
                   Precalc->dQuadStd);

        J2 = strike * norm((mean - root_1) / Precalc->dQuadStd);

        *Premium = J1 - J2;
      } else {
        J1 = exp(constant_coef - mean * mean / (2.0 * Precalc->dQuadVar) +
                 (Precalc->dQuadLin * Precalc->dQuadVar + mean) *
                     (Precalc->dQuadLin * Precalc->dQuadVar + mean) /
                     (2.0 * Precalc->dQuadVar));
        J1 *= norm((root_1 - mean - Precalc->dQuadLin * Precalc->dQuadVar) /
                   Precalc->dQuadStd);

        J2 = strike * norm((root_1 - mean) / Precalc->dQuadStd);

        *Premium = J2 - J1;
      }
    } else {
      delta = Precalc->dQuadLin * Precalc->dQuadLin -
              4 * Precalc->dQuadQuad * (constant_coef - LogStrike);
      if (delta < 0) {
        if (callput == SRT_CALL) {
          if (Precalc->dQuadQuad < 0) {
            *Premium = 0;
          } else {
            k1 = Precalc->dQuadLin + mean / Precalc->dQuadVar;
            k2 = sqrt(1.0 / Precalc->dQuadVar - 2.0 * Precalc->dQuadQuad);
            *Premium =
                exp(constant_coef - mean * mean / (2.0 * Precalc->dQuadVar) +
                    k1 * k1 / (2.0 * k2 * k2)) /
                Precalc->dQuadStd / k2;
            *Premium += -strike;
          }
        } else {
          if (Precalc->dQuadQuad > 0) {
            *Premium = 0;
          } else {
            k1 = Precalc->dQuadLin + mean / Precalc->dQuadVar;
            k2 = sqrt(1.0 / Precalc->dQuadVar - 2.0 * Precalc->dQuadQuad);
            *Premium = strike;
            *Premium +=
                -exp(constant_coef - mean * mean / (2.0 * Precalc->dQuadVar) +
                     k1 * k1 / (2.0 * k2 * k2)) /
                Precalc->dQuadStd / k2;
          }
        }
      }

      else {
        root_1 =
            (-Precalc->dQuadLin + sqrt(delta)) / (2.0 * Precalc->dQuadQuad);
        root_2 =
            (-Precalc->dQuadLin - sqrt(delta)) / (2.0 * Precalc->dQuadQuad);

        /* put the roots in the right order */
        if (root_1 > root_2) {
          swap_root = root_1;
          root_1 = root_2;
          root_2 = swap_root;
        }
        k1 = Precalc->dQuadLin + mean / Precalc->dQuadVar;
        k2 = sqrt(1.0 / Precalc->dQuadVar - 2.0 * Precalc->dQuadQuad);

        if (Precalc->dQuadQuad > 0) {
          if (callput == SRT_CALL) {
            J1 = norm(k2 * root_1 - k1 / k2);
            J1 += norm(k1 / k2 - k2 * root_2);

            J2 = norm((root_1 - mean) / Precalc->dQuadStd);
            J2 += norm((mean - root_2) / Precalc->dQuadStd);
          } else {
            J1 = norm(k2 * root_2 - k1 / k2);
            J1 += -norm(k2 * root_1 - k1 / k2);

            J2 = norm((root_2 - mean) / Precalc->dQuadStd);
            J2 += -norm((root_1 - mean) / Precalc->dQuadStd);
          }
        } else {
          if (callput == SRT_CALL) {
            J1 = norm(k2 * root_2 - k1 / k2);
            J1 += -norm(k2 * root_1 - k1 / k2);

            J2 = norm((root_2 - mean) / Precalc->dQuadStd);
            J2 += -norm((root_1 - mean) / Precalc->dQuadStd);
          } else {
            J1 = norm(k2 * root_1 - k1 / k2);
            J1 += norm(k1 / k2 - k2 * root_2);

            J2 = norm((root_1 - mean) / Precalc->dQuadStd);
            J2 += norm((mean - root_2) / Precalc->dQuadStd);
          }
        }

        J1 *= exp(constant_coef - mean * mean / (2.0 * Precalc->dQuadVar) +
                  k1 * k1 / (2 * k2 * k2)) /
              Precalc->dQuadStd / k2;
        J2 *= strike;
        if (callput == SRT_CALL) {
          *Premium = J1 - J2;
        } else {
          *Premium = J2 - J1;
        }
      }
    }
  }
FREE_RETURN:

  return err;
}

Err FxBetaDLM_Price_ForwardFx(FxBetaDLM_InstPrecalc *Precalc, double *Premium) {
  Err err = NULL;

  *Premium = Precalc->dFwdBeta *
             exp(Precalc->dQuadConst +
                 Precalc->dFwdQuad2 * (Precalc->dQuadLin + Precalc->dFwdAlpha) *
                     (Precalc->dQuadLin + Precalc->dFwdAlpha));

  return err;
}

Err FxBetaDLM_GetEqui3FactorTS(FxBetaDLM_model *Model,
                               FxBetaDLM_OptNumerParams *NumParams,
                               FxBetaDLM_Hermite *hermite, int nb_fx,
                               double *fx_maturities, double **fx_cal_vols) {
  Err err = NULL;
  int i;
  long settlmt_date, fixing_date;
  double strike, premium;
  double *fx_vols = NULL, *settlmt_time = NULL;
  FxBetaDLM_FxOptInst *Inst = NULL;
  FxBetaDLM_InstPrecalc *InstConst = NULL;

  fx_vols = calloc(nb_fx, sizeof(double));
  settlmt_time = calloc(nb_fx, sizeof(double));

  Inst = calloc(1, sizeof(FxBetaDLM_FxOptInst));
  InstConst = calloc(1, sizeof(FxBetaDLM_InstPrecalc));

  if (!fx_vols || !settlmt_time || !Inst || !InstConst) {
    err = "Memory allocation faillure in FxBetaDLM_GetEqui3FactorTS";
    goto FREE_RETURN;
  }

  /* Setup the instrument */
  err = FxBetaDLM_Allocate_FxOptInst(0, 1, Inst);

  if (err)
    goto FREE_RETURN;

  /* Get the ATM vols */
  for (i = 0; i < nb_fx; i++) {
    fixing_date = (long)(Model->lToday + fx_maturities[i] * DAYS_IN_YEAR + 0.5);
    settlmt_date = add_unit(fixing_date, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
    settlmt_time[i] = (settlmt_date - Model->lToday) * YEARS_IN_DAY;

    strike = 0.0;

    err = FxBetaDLM_Setup_FxOptInst(settlmt_date, settlmt_date, Model->lToday,
                                    fixing_date, 1, &strike, NULL, SRT_CALL,
                                    Model, Inst);

    if (err)
      goto FREE_RETURN;

    err = FxBetaDLM_Calculate_AllConst(Inst, Model, NumParams, InstConst);

    if (err)
      goto FREE_RETURN;

    err = FxBetaDLM_Update_InstPrecalc_FromMoment(Inst, Model, NumParams,
                                                  InstConst);

    if (err)
      goto FREE_RETURN;

    err = FxBetaDLM_Price_FxOptInst(Inst, Model, InstConst, hermite, NumParams,
                                    &premium);

    if (err)
      goto FREE_RETURN;

    err = srt_f_optimpvol(premium, Inst->dFwdFx, Inst->dStrike[0],
                          Inst->dVolTime, Inst->dDfPayDom, SRT_CALL,
                          SRT_LOGNORMAL, &(fx_vols[i]));

    if (err)
      goto FREE_RETURN;
  }

  /* Calibrate the 3 Factor */
  err = Fx3DtsCalibration_corr(
      fx_maturities, settlmt_time, fx_vols, nb_fx, Model->dPWTime,
      Model->iNbPWTime, Model->dSigmaDom, Model->dLambdaDom, Model->dSigmaFor,
      Model->dLambdaFor, Model->dPWTime, Model->dCorrDomFor,
      Model->dCorrDomFxInput, Model->dCorrForFxInput, Model->iNbPWTime,
      fx_cal_vols);

  if (err)
    goto FREE_RETURN;

FREE_RETURN:

  if (Inst) {
    FxBetaDLM_Free_FxOptInst(Inst);
    free(Inst);
  }

  if (InstConst) {
    free(InstConst);
  }

  if (fx_vols)
    free(fx_vols);
  if (settlmt_time)
    free(settlmt_time);

  return err;
}