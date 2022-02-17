#include "LGM2FMC.h"
#include "LGMSVUtil.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"

#include "TargetNote.h"
#include "TargetNoteCalib.h"
#include "TargetNoteProdStruct.h"
#include "TargetNoteProdStructSV.h"
#include "TargetNoteStatic.h"

#define BIG_BIG_NUMBER 1E40
#define ONE_MONTH 0.083333333

/* ------------------------------------------------------------------------------------------------------------------
 */
/* deal definition structures */

char *init_TARN_AUX(TARN_Struct *tarn, TARN_AUX *aux) {
  /* Variable declaration */
  int i = 0;
  int iFund = 0;
  double dDF = 0.0;
  char *err = 0;
  long lFirstFixing = tarn->market.lToday + !tarn->market.eodFixFlag;
  long lFirstPayment = tarn->market.lToday + !tarn->market.eodPayFlag;

  /* The coupon schedule */
  aux->dvCouponCvg = 0;
  aux->dvCumulCvg = 0;
  aux->dvCouponPayDF = 0;
  aux->dvKnockOutFee = 0;
  /* The funding */
  aux->ivFundStartIndex = 0;
  aux->ivFundEndIndex = 0;
  aux->dvFundCvg = 0;
  aux->dvFundDomPayDF = 0;
  aux->dvFundDomStartDF = 0;
  aux->dvFundForStartDF = 0;

  /* Generate the coverages and DF */
  aux->dvCouponCvg = (double *)calloc(tarn->deal.nCouponDates, sizeof(double));
  aux->dvCumulCvg = (double *)calloc(tarn->deal.nCouponDates, sizeof(double));
  aux->dvCouponPayDF =
      (double *)calloc(tarn->deal.nCouponDates, sizeof(double));

  aux->dvKnockOutFee =
      (double *)calloc(tarn->deal.nCouponDates, sizeof(double));
  for (i = 0; i < tarn->deal.nCouponDates; i++) {
    aux->dvKnockOutFee[i] = 0.0;
    if (tarn->deal.dFixedCvg > 0.0001)
      aux->dvCouponCvg[i] = tarn->deal.dFixedCvg;
    else
      aux->dvCouponCvg[i] =
          coverage(tarn->deal.lvCpnStartDates[i], tarn->deal.lvCpnPayDates[i],
                   tarn->deal.basisCpn);
    aux->dvCumulCvg[i] = aux->dvCouponCvg[i];
    if (tarn->deal.lvCpnPayDates[i] >= tarn->market.lToday) {
      if (err =
              tarn->market.getDF(tarn->market.szYieldCurve, tarn->market.lToday,
                                 (double)tarn->deal.lvCpnPayDates[i], &dDF))
        return err;
      aux->dvCouponPayDF[i] = dDF;
    } else
      aux->dvCouponPayDF[i] = 0.0;
  }

  aux->dvFundCvg = (double *)calloc(tarn->deal.nFundDates, sizeof(double));
  aux->dvFundDomPayDF = (double *)calloc(tarn->deal.nFundDates, sizeof(double));
  aux->dvFundDomStartDF =
      (double *)calloc(tarn->deal.nFundDates, sizeof(double));
  aux->dvFundForStartDF =
      (double *)calloc(tarn->deal.nFundDates, sizeof(double));
  for (i = 0; i < tarn->deal.nFundDates; i++) {
    aux->dvFundCvg[i] =
        coverage(tarn->deal.lvFundStartDates[i], tarn->deal.lvFundEndDates[i],
                 tarn->deal.basisFund);

    if (tarn->deal.lvFundStartDates[i] >= tarn->market.lToday) {
      if (err =
              tarn->market.getDF(tarn->market.szYieldCurve, tarn->market.lToday,
                                 (double)tarn->deal.lvFundStartDates[i], &dDF))
        return err;
      aux->dvFundDomStartDF[i] = dDF;
      if (err = tarn->market.getDF(
              tarn->market.szFundYieldCurve, tarn->market.lToday,
              (double)tarn->deal.lvFundStartDates[i], &dDF))
        return err;
      aux->dvFundForStartDF[i] = dDF;
    } else {
      aux->dvFundDomStartDF[i] = 0.0;
      aux->dvFundForStartDF[i] = 0.0;
    }

    if (tarn->deal.lvFundEndDates[i] >= tarn->market.lToday) {
      if (err =
              tarn->market.getDF(tarn->market.szYieldCurve, tarn->market.lToday,
                                 (double)tarn->deal.lvFundEndDates[i], &dDF))
        return err;
      aux->dvFundDomPayDF[i] = dDF;
    } else
      aux->dvFundDomPayDF[i] = 0.0;
  }

  /* Loop over the coupons calculating the floating periods that can be
     calculated . A floating payment belongs to period n if it starts strictly
     before the nth coupon payment date */
  iFund = 0;
  aux->ivFundStartIndex = (int *)calloc(tarn->deal.nCouponDates, sizeof(int));
  aux->ivFundEndIndex = (int *)calloc(tarn->deal.nCouponDates, sizeof(int));
  for (i = 0; i < tarn->deal.nCouponDates - 1; i++) {
    aux->ivFundStartIndex[i] = iFund;
    while (iFund < tarn->deal.nFundDates &&
           tarn->deal.lvFundStartDates[iFund] < tarn->deal.lvCpnPayDates[i])
      iFund++;
    aux->ivFundEndIndex[i] = iFund - 1;
  }

  /* The last period gets the final bit */
  if (tarn->deal.nCouponDates == 1) {
    aux->ivFundStartIndex[0] = 0;
    aux->ivFundEndIndex[0] = tarn->deal.nFundDates - 1;
  } else {
    aux->ivFundStartIndex[tarn->deal.nCouponDates - 1] =
        aux->ivFundEndIndex[tarn->deal.nCouponDates - 2] + 1;
    aux->ivFundEndIndex[tarn->deal.nCouponDates - 1] =
        tarn->deal.nFundDates - 1;
  }

  /* historical variables */
  aux->i1stCpn = 0;
  aux->dHistFundPV = 0.0;
  aux->dHistCpnPV = 0.0;
  aux->dHistCumulCpn = 0.0;

  /* check if today is a fixing date */
  aux->iIsTodayFixing = 0;
  for (i = 0; i < tarn->deal.nCouponDates; i++)
    if (tarn->deal.lvFltrFixingDates[i] == tarn->market.lToday)
      aux->iIsTodayFixing = 1;

  /* Convert the foreign funding to domestic (if necessary) */
  if (tarn->deal.fund_ccy == 1) {
    //		long lStartDate;
    //		double eq_final_ex  , eq_init_ex  , dFundingPV  , dCouponPV;
    double dFundDF, dDomDF;
    double dNotRatio = tarn->deal.dFundingNotional / tarn->deal.dCpnNotional *
                       tarn->market.fx_fund_dom;

    /* Find the first funding coupon not yet fixed */
    i = 0;
    while (i < tarn->deal.nFundDates &&
           tarn->deal.lvFundFixingDates[i] <= lFirstFixing)
      ++i;

    /* add the notional exchange for the first funding date to the Funding PV */
    /*		if ( i < tarn->deal.nFundDates )
                    {
                            if ( err = tarn->market.getDF(
       tarn->market.szFundYieldCurve  , tarn->market.lToday  ,
                                                                                            tarn->deal.lvFundStartDates[i]  , &dFundDF ) )
                                    goto FREE_RETURN;
                            if ( err = tarn->market.getDF(
       tarn->market.szYieldCurve  , tarn->market.lToday  ,
                                                                                            tarn->deal.lvFundStartDates[i]  , &dDomDF ) )
                                    goto FREE_RETURN;
                            aux->dHistFundPV = tarn->deal.dFundingNotional /
       tarn->deal.dCpnNotional * dFundDF * tarn->market.fx_fund_dom
                                                                            -
       dDomDF;
                    }
    */
    /*
                    if ( err = tarn->market.getDF( tarn->market.szFundYieldCurve
       , tarn->market.lToday  ,
                                                                                    tarn->deal.lvFundStartDates[aux->ivFundStartIndex[aux->i1stCpn+1]]  , &dFundDF ) )
                            goto FREE_RETURN;
                    if ( err = tarn->market.getDF( tarn->market.szYieldCurve  ,
       tarn->market.lToday  ,
                                                                                    tarn->deal.lvFundStartDates[aux->ivFundStartIndex[aux->i1stCpn+1]]  , &dDomDF ) )
                            goto FREE_RETURN;
                    aux->dHistFundPV = tarn->deal.dFundingNotional /
       tarn->deal.dCpnNotional * dFundDF * tarn->market.fx_fund_dom
                                                                    - dDomDF;
    */

    /* calculate the potential notional exchange at each of the KO settle dates
     * and store it in the KO fee vector */
    /*		for ( i=aux->i1stCpn; i<tarn->deal.nCouponDates; i++ )
                    {

                            if (
       tarn->deal.lvFundEndDates[aux->ivFundEndIndex[i]] >= tarn->market.lToday
       )
                            {
                                    if ( err = tarn->market.getDF(
       tarn->market.szFundYieldCurve  , tarn->market.lToday  ,
                                                                                                    tarn->deal.lvFundEndDates[aux->ivFundEndIndex[i]]  , &dFundDF ) )
                                            goto FREE_RETURN;
                                    if ( err = tarn->market.getDF(
       tarn->market.szYieldCurve  , tarn->market.lToday  ,
                                                                                                    tarn->deal.lvFundEndDates[aux->ivFundEndIndex[i]]  , &dDomDF ) )
                                            goto FREE_RETURN;
                                    aux->dvKnockOutFee[i] = 1.0 -
       tarn->deal.dFundingNotional / tarn->deal.dCpnNotional * dFundDF / dDomDF
                                                                                            * tarn->market.fx_fund_dom;
                            }
                    }
    */
    /* Calculate the effective margin in the domestic currency */
    for (i = 0; i < tarn->deal.nFundDates; i++) {
      if (tarn->deal.lvFundFixingDates[i] < lFirstFixing) {
        if (tarn->deal.lvFundEndDates[i] >= tarn->market.lToday) {
          if (err = tarn->market.getDF(tarn->market.szFundYieldCurve,
                                       tarn->market.lToday,
                                       tarn->deal.lvFundEndDates[i], &dFundDF))
            goto FREE_RETURN;
          if (err = tarn->market.getDF(tarn->market.szYieldCurve,
                                       tarn->market.lToday,
                                       tarn->deal.lvFundEndDates[i], &dDomDF))
            goto FREE_RETURN;
          tarn->deal.dvFundHistFixings[i] =
              (tarn->deal.dvFundHistFixings[i] + tarn->deal.dvFundMargin[i]) *
              dFundDF / dDomDF * dNotRatio;
          tarn->deal.dvFundMargin[i] = 0.0;
        }
      }

      if (tarn->deal.lvFundFixingDates[i] >= lFirstFixing) {
        if (err = tarn->market.getDF(tarn->market.szFundYieldCurve,
                                     tarn->market.lToday,
                                     tarn->deal.lvFundEndDates[i], &dFundDF))
          goto FREE_RETURN;
        if (err = tarn->market.getDF(tarn->market.szYieldCurve,
                                     tarn->market.lToday,
                                     tarn->deal.lvFundEndDates[i], &dDomDF))
          goto FREE_RETURN;
        tarn->deal.dvFundMargin[i] =
            (tarn->deal.dvFundMargin[i] + tarn->deal.dvFundSpread[i]) *
            dFundDF / dDomDF * dNotRatio;
        tarn->deal.dvFundSpread[i] = 0.0;
      }
    }

    /* Convert into domestic */
    tarn->deal.fund_ccy = 0;
    //		TARGET_NOTE.m_UseNotionalExchange = 1;
    /*		if ( err = convert_funding_to_domestic(
                                                                    tarn->market.lToday
       , tarn->deal.lvFundStartDates[aux->i1stCpn]  , tarn->market.eodFixFlag  ,
                                                                    tarn->market.eodPayFlag
       , tarn->market.fx_fund_dom  , tarn->market.fx_fund_dom_spot_date  ,
                                                                    tarn->deal.dCpnNotional
       , tarn->market.szYieldCurve  , tarn->deal.nFundDates  ,
                                                                    tarn->deal.lvFundFixingDates
       , tarn->deal.lvFundStartDates  , tarn->deal.lvFundEndDates  ,
                                                                    tarn->deal.szvBasisFund
       , tarn->market.szFundYieldCurve  , &tarn->deal.dFundingNotional  ,
                                                                    tarn->deal.dvFundSpread
       , tarn->deal.dvFundMargin  , tarn->deal.dvFundHistFixings  , &lStartDate
       , &eq_final_ex  , &eq_init_ex)  ) goto FREE_RETURN;

                    TARGET_NOTE.m_dFundingNotionalExchange = eq_init_ex -
       eq_final_ex; TARGET_NOTE.m_dCouponNotionalExchange =
       tarn->deal.dCpnNotional;
    *//*
		if ( lStartDate >= tarn->market.lToday + tarn->market.eodPayFlag )
		{
			double dStartDF;
			if ( err = tarn->market.getDF( tarn->market.szYieldCurve  , tarn->market.lToday  , lStartDate  , &dStartDF ) )
				goto FREE_RETURN;
			dFundingPV= dStartDF * TARGET_NOTE.m_dFundingNotionalExchange;
			dCouponPV = - dStartDF * TARGET_NOTE.m_dCouponNotionalExchange;
		}
*/
  }
  //	else
  //		TARGET_NOTE.m_UseNotionalExchange = 0;

FREE_RETURN:
  /* return */
  return err;
}

void free_TARN_AUX(TARN_AUX *aux) {
  /* Free the allocated memory */
  if (aux->dvCouponCvg)
    free_and_zero(&aux->dvCouponCvg);
  if (aux->dvCouponPayDF)
    free_and_zero(&aux->dvCouponPayDF);
  if (aux->dvFundCvg)
    free_and_zero(&aux->dvFundCvg);
  free_and_zero(&aux->dvFundDomPayDF);
  free_and_zero(&aux->dvFundDomStartDF);
  free_and_zero(&aux->dvFundForStartDF);
  if (aux->ivFundStartIndex)
    free_and_zero(&aux->ivFundStartIndex);
  if (aux->ivFundEndIndex)
    free_and_zero(&aux->ivFundEndIndex);
  free_and_zero(&aux->dvKnockOutFee);
}
/* deal definition structures */
/* ------------------------------------------------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------------------------------------------
 */
/* Setup function for the global variable */
void initTargetNote(TARN_Struct *tarn) {
  TARGET_NOTE.m_dTarget = tarn->deal.dTarget;
  TARGET_NOTE.m_couponType = tarn->deal.couponType;
  TARGET_NOTE.m_dSpread = tarn->pricing.dSpread;
  TARGET_NOTE.m_dHistCumulCoupon = tarn->output.dCumulCoupon;
  TARGET_NOTE.tarn = tarn;
  TARGET_NOTE.m_dFinalCoupon = tarn->output.dAccretedCoupon;
  TARGET_NOTE.m_dCouponNotionalRepaid = tarn->output.iIsKnockedOut > 0 ? 1 : 0;
  TARGET_NOTE.m_dPrevAccrual = 1.0;
}

void resetTargetNote() {
  // set the cumulative coupon to zero
  TARGET_NOTE.m_dCumulCoupon = TARGET_NOTE.m_dHistCumulCoupon;
  TARGET_NOTE.m_dCouponNotional = 1.0;
  TARGET_NOTE.m_dFundingNotional = 1.0;
  TARGET_NOTE.m_dFinalCoupon = TARGET_NOTE.tarn->output.dAccretedCoupon;
  TARGET_NOTE.m_dCouponNotionalRepaid =
      TARGET_NOTE.tarn->output.iIsKnockedOut > 0 ? 1 : 0;
  TARGET_NOTE.m_dPrevAccrual = 1.0;
}
/* Setup function for the global variable */
/* -----------------------------------------------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------------------------------------------
 */
/* Memory management for TARN_MC_AUX */

void init_TARN_MC_AUX(TARN_MC_AUX *mc_aux, char *szUnd) {
  mc_aux->date = 0;
  mc_aux->time = 0;
  mc_aux->evt_tms = 0;
  mc_aux->evt_dts = 0;
  mc_aux->is_event = 0;
  mc_aux->ptrModel = 0;
  mc_aux->void_prm = 0;
  mc_aux->covar = 0;
  mc_aux->und_name = 0;
  mc_aux->szUnd = szUnd;
}

void free_TARN_MC_AUX(TARN_MC_AUX *mc_aux) {
  free_and_zero(&mc_aux->date);
  free_and_zero(&mc_aux->time);
  free_and_zero(&mc_aux->evt_tms);
  free_and_zero(&mc_aux->evt_dts);
  free_and_zero(&mc_aux->is_event);
  if (mc_aux->tnModel)
    switch (mc_aux->tnModel) {
    case TN_2F:
      free_TARN_MC_AUX_2F((TARN_MC_AUX_2F *)mc_aux->ptrModel);
      break;
    default:
      free_TARN_MC_AUX_SV((TARN_MC_AUX_SV *)mc_aux->ptrModel);
      break;
    }

  if (mc_aux->void_prm) {
    int i;
    TARN_EVENT *event = 0;
    for (i = 0; i < mc_aux->nStps; i++)
      if (mc_aux->void_prm[i]) {
        event = (TARN_EVENT *)mc_aux->void_prm[i];
        free_TARN_EVENT(event);
        free_and_zero(&event);
      }
    free(mc_aux->void_prm);
  }
}

char *alloc_TARN_MC_AUX(TARN_Struct *tarn, TARN_AUX *aux, TARN_MC_AUX *mc_aux) {
  /* Variable declaration */
  SrtUndPtr und;
  char *err = 0;

  /* look for the underlying */
  if (!(und = lookup_und(mc_aux->szUnd))) {
    err = "cannot find the underlying";
    goto FREE_RETURN;
  }
  mc_aux->und_name = und->underl_name;

  /* check whether we have SV or LGM */
  if (get_mdltype_from_irund(und) == LGM) {
    alloc_TARN_MC_AUX_2F(tarn, aux, mc_aux);
  } else if (get_mdltype_from_irund(und) == LGM_STOCH_VOL) {
    alloc_TARN_MC_AUX_SV(tarn, aux, mc_aux);
  } else {
    err = "Unknown model type";
    goto FREE_RETURN;
  }

FREE_RETURN:
  return err;
}

/* Memory management for TARN_MC_AUX */
/* -----------------------------------------------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------------------------------------------
 */
/* Memory management for TARN_MC_AUX_2F */

char *init_TARN_MC_AUX_2F(TARN_MC_AUX_2F **mc_aux_2f) {
  *mc_aux_2f = (TARN_MC_AUX_2F *)malloc(sizeof(TARN_MC_AUX_2F));
  if (!*mc_aux_2f)
    return "Memory allocation error in TARN_MC_AUX_2F";
  (*mc_aux_2f)->ifr = 0;
  (*mc_aux_2f)->fwd1 = 0;
  (*mc_aux_2f)->fwd2 = 0;
  (*mc_aux_2f)->exp1 = 0;
  (*mc_aux_2f)->exp2 = 0;
  (*mc_aux_2f)->phi1 = 0;
  (*mc_aux_2f)->phi2 = 0;
  (*mc_aux_2f)->phi12 = 0;
  (*mc_aux_2f)->gam1_fwd = 0;
  (*mc_aux_2f)->gam2_fwd = 0;
  (*mc_aux_2f)->bond_pay = 0;
  (*mc_aux_2f)->gam1_pay = 0;
  (*mc_aux_2f)->gam2_pay = 0;

  return 0;
}

void free_TARN_MC_AUX_2F(TARN_MC_AUX_2F *mc_aux_2f) {
  free_and_zero(&mc_aux_2f->ifr);
  free_and_zero(&mc_aux_2f->fwd1);
  free_and_zero(&mc_aux_2f->fwd2);
  free_and_zero(&mc_aux_2f->exp1);
  free_and_zero(&mc_aux_2f->exp2);
  free_and_zero(&mc_aux_2f->phi1);
  free_and_zero(&mc_aux_2f->phi2);
  free_and_zero(&mc_aux_2f->phi12);
  free_and_zero(&mc_aux_2f->gam1_fwd);
  free_and_zero(&mc_aux_2f->gam2_fwd);
  free_and_zero(&mc_aux_2f->bond_pay);
  free_and_zero(&mc_aux_2f->gam1_pay);
  free_and_zero(&mc_aux_2f->gam2_pay);
}

char *alloc_TARN_MC_AUX_2F(TARN_Struct *tarn, TARN_AUX *aux,
                           TARN_MC_AUX *mc_aux) {
  /* Variable declaration */
  double tau, alpha, gamma, rho;
  int num_sig;
  double *sig_date = 0, *sig = 0;
  TARN_MC_AUX_2F *mc_aux_2f = 0;
  char *err = 0;

  /* allocate the memory */
  if (err = init_TARN_MC_AUX_2F(&mc_aux_2f))
    goto FREE_RETURN;

  /*	Next  , set time steps at each of the floater fixing dates */
  if (err = TARN_TimeSteps_2F(tarn, aux, mc_aux))
    goto FREE_RETURN;

  /* Get the term structure */
  if (err = Get_LGM2F_TermStructure(mc_aux->und_name, &sig_date, &sig, &num_sig,
                                    &tau, &alpha, &gamma, &rho))
    goto FREE_RETURN;

  /* set the constant values */
  mc_aux_2f->lam = 1.0 / tau;
  mc_aux_2f->lam2 = mc_aux_2f->lam + gamma;

  /*	Allocate the memory */
  if (!(mc_aux_2f->ifr = (double *)calloc(mc_aux->nStps, sizeof(double))))
    return "Memory allocation error";
  if (!(mc_aux_2f->fwd1 = (double *)calloc(mc_aux->nStps, sizeof(double))))
    return "Memory allocation error";
  if (!(mc_aux_2f->fwd2 = (double *)calloc(mc_aux->nStps, sizeof(double))))
    return "Memory allocation error";
  if (!(mc_aux_2f->exp1 = (double *)calloc(mc_aux->nStps, sizeof(double))))
    return "Memory allocation error";
  if (!(mc_aux_2f->exp2 = (double *)calloc(mc_aux->nStps, sizeof(double))))
    return "Memory allocation error";
  if (!(mc_aux_2f->phi1 = (double *)calloc(mc_aux->nStps, sizeof(double))))
    return "Memory allocation error";
  if (!(mc_aux_2f->phi2 = (double *)calloc(mc_aux->nStps, sizeof(double))))
    return "Memory allocation error";
  if (!(mc_aux_2f->phi12 = (double *)calloc(mc_aux->nStps, sizeof(double))))
    return "Memory allocation error";
  if (!(mc_aux_2f->gam1_fwd = (double *)calloc(mc_aux->nStps, sizeof(double))))
    return "Memory allocation error";
  if (!(mc_aux_2f->gam2_fwd = (double *)calloc(mc_aux->nStps, sizeof(double))))
    return "Memory allocation error";
  if (!(mc_aux_2f->bond_pay = (double *)calloc(mc_aux->nStps, sizeof(double))))
    return "Memory allocation error";
  if (!(mc_aux_2f->gam1_pay = (double *)calloc(mc_aux->nStps, sizeof(double))))
    return "Memory allocation error";
  if (!(mc_aux_2f->gam2_pay = (double *)calloc(mc_aux->nStps, sizeof(double))))
    return "Memory allocation error";

  /* allocate the memory for the covariance matrix */
  if (!(mc_aux->covar = f3tensor(0, mc_aux->nStps - 1, 0, 1, 0, 1))) {
    err = "Memory allocation error";
    goto FREE_RETURN;
  }

  mc_aux_2f->pay_date = (long)mc_aux->date[mc_aux->nStps - 1];
  mc_aux_2f->pay_time = mc_aux->time[mc_aux->nStps - 1];

  /* calculate the model quantities */
  fill_mc_init_lgm2f(
      0, mc_aux_2f->pay_date, mc_aux_2f->pay_time, mc_aux->date, mc_aux->time,
      mc_aux->nStps, sig_date, num_sig, sig, mc_aux_2f->lam, alpha, gamma, rho,
      tarn->market.szYieldCurve, mc_aux_2f->fwd1, mc_aux_2f->fwd2,
      mc_aux_2f->exp1, mc_aux_2f->exp2, mc_aux_2f->phi1, mc_aux_2f->phi2,
      mc_aux_2f->phi12, mc_aux_2f->gam1_fwd, mc_aux_2f->gam2_fwd,
      mc_aux_2f->bond_pay, mc_aux_2f->gam1_pay, mc_aux_2f->gam2_pay,
      mc_aux->covar);

  /* assign the 2F to the mc aux */
  mc_aux->tnModel = TN_2F;
  mc_aux->ptrModel = mc_aux_2f;

/* return safely */
FREE_RETURN:
  return 0;
}

/* Memory management for TARN_MC_AUX_2F */
/* -----------------------------------------------------------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------------------------------------------
 */
/* Memory management for 2F event */

void init_TARN_EVENT_2F(TARN_EVENT_2F *event) {
  memset(event, 0, sizeof(TARN_EVENT_2F));
  /*
          event->fund_start_gam1 = 0;
          event->fund_start_gam2 = 0;
          event->fund_start_gam12 = 0;
          event->fund_end_gam1 = 0;
          event->fund_end_gam2 = 0;
          event->fund_end_gam12 = 0;
  */
}

void free_TARN_EVENT_2F(TARN_EVENT_2F *event) {
  if (event->fund_start_gam1)
    free_dvector(event->fund_start_gam1, 0, event->num_fund_df - 1);
  if (event->fund_start_gam2)
    free_dvector(event->fund_start_gam2, 0, event->num_fund_df - 1);
  if (event->fund_start_gam12)
    free_dvector(event->fund_start_gam12, 0, event->num_fund_df - 1);
  if (event->fund_end_gam1)
    free_dvector(event->fund_end_gam1, 0, event->num_fund_df - 1);
  if (event->fund_end_gam2)
    free_dvector(event->fund_end_gam2, 0, event->num_fund_df - 1);
  if (event->fund_end_gam12)
    free_dvector(event->fund_end_gam12, 0, event->num_fund_df - 1);
}

/* Memory management for 2F structure */
/* ----------------------------------------------------------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------------------------------------------
 */
/* Memory management for SV structure */

void init_TARN_EVENT_SV(TARN_EVENT_SV *event) {
  memset(event, 0, sizeof(TARN_EVENT_SV));
  /*
          event->fund_start_gam1 = 0;
          event->fund_start_gam2 = 0;
          event->fund_start_gam1_2 = 0;
          event->fund_start_gam2_2 = 0;
          event->fund_start_gam12 = 0;
          event->fund_end_gam1 = 0;
          event->fund_end_gam2 = 0;
          event->fund_end_gam1_2 = 0;
          event->fund_end_gam2_2 = 0;
          event->fund_end_gam12 = 0;
  */
}

void free_TARN_EVENT_SV(TARN_EVENT_SV *event) {
  if (event->fund_start_gam1)
    free_dvector(event->fund_start_gam1, 0, event->num_fund_df - 1);
  if (event->fund_start_gam2)
    free_dvector(event->fund_start_gam2, 0, event->num_fund_df - 1);
  if (event->fund_start_gam1_2)
    free_dvector(event->fund_start_gam1_2, 0, event->num_fund_df - 1);
  if (event->fund_start_gam2_2)
    free_dvector(event->fund_start_gam2_2, 0, event->num_fund_df - 1);
  if (event->fund_start_gam12)
    free_dvector(event->fund_start_gam12, 0, event->num_fund_df - 1);
  if (event->fund_end_gam1)
    free_dvector(event->fund_end_gam1, 0, event->num_fund_df - 1);
  if (event->fund_end_gam2)
    free_dvector(event->fund_end_gam2, 0, event->num_fund_df - 1);
  if (event->fund_end_gam1_2)
    free_dvector(event->fund_end_gam1_2, 0, event->num_fund_df - 1);
  if (event->fund_end_gam2_2)
    free_dvector(event->fund_end_gam2_2, 0, event->num_fund_df - 1);
  if (event->fund_end_gam12)
    free_dvector(event->fund_end_gam12, 0, event->num_fund_df - 1);
}

/* Memory management for 2FSV structure */
/* ----------------------------------------------------------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------------------------------------------
 */
/* General Structure */
void init_TARN_EVENT(TARN_EVENT *event) {
  memset(event, 0, sizeof(TARN_EVENT));
  /*
          event->fund_start = 0;
          event->fund_end = 0;
          event->fund_start_tms = 0;
          event->fund_start_dts = 0;
          event->fund_start_logdff = 0;
          event->fund_end_tms = 0;
          event->fund_end_dts = 0;
          event->fund_end_logdff = 0;
          event->fund_cvgMargin = 0;
          event->ptrModel = 0;
  //	event->do_dom = 1;
  //	event->bIsFirst = 0;
          event->iIsFinal = 0;
  */
}

void free_TARN_EVENT(TARN_EVENT *event) {
  if (event->fund_start_tms)
    free_dvector(event->fund_start_tms, 0, event->num_fund_df - 1);
  if (event->fund_start_dts)
    free_dvector(event->fund_start_dts, 0, event->num_fund_df - 1);
  if (event->fund_start_logdff)
    free_dvector(event->fund_start_logdff, 0, event->num_fund_df - 1);
  if (event->fund_start)
    free_dvector(event->fund_start, 0, event->num_fund_df - 1);
  if (event->fund_end_tms)
    free_dvector(event->fund_end_tms, 0, event->num_fund_df - 1);
  if (event->fund_end_dts)
    free_dvector(event->fund_end_dts, 0, event->num_fund_df - 1);
  if (event->fund_end_logdff)
    free_dvector(event->fund_end_logdff, 0, event->num_fund_df - 1);
  if (event->fund_end)
    free_dvector(event->fund_end, 0, event->num_fund_df - 1);
  if (event->fund_cvgMargin)
    free_dvector(event->fund_cvgMargin, 0, event->num_fund_df - 1);
  if (event->ptrModel)
    switch (event->iModelType) {
    case TN_2F:
      free_TARN_EVENT_2F((TARN_EVENT_2F *)event->ptrModel);
      break;
    default:
      free_TARN_EVENT_SV((TARN_EVENT_SV *)event->ptrModel);
      break;
    }
}

char *alloc_TARN_EVENT(TARN_EVENT **event, TARN_MC_AUX *mc_aux) {
  /* Variable declaration */
  TARN_EVENT_2F *event_2F = 0;
  TARN_EVENT_SV *event_sv = 0;

  /* Allocate memory for the general structure */
  if (!(*event = (TARN_EVENT *)malloc(sizeof(TARN_EVENT))))
    return "Memory allocation error in allocTargetNoteStruct";
  init_TARN_EVENT(*event);
  (*event)->iModelType = mc_aux->tnModel;

  /* Allocate memory for the model specific data */
  switch (mc_aux->tnModel) {
  case TN_2F:
    if (!(event_2F = (TARN_EVENT_2F *)malloc(sizeof(TARN_EVENT_2F))))
      return "Memory allocation error in allocTargetNoteStruct";
    init_TARN_EVENT_2F(event_2F);
    (*event)->ptrModel = (void *)event_2F;
    break;
  default:
    if (!(event_sv = (TARN_EVENT_SV *)malloc(sizeof(TARN_EVENT_SV))))
      return "Memory allocation error in allocTargetNoteStruct";
    init_TARN_EVENT_SV(event_sv);
    (*event)->ptrModel = (void *)event_sv;
    break;
  }

  /* return safely */
  return 0;
}
/* General Structure */
/* ----------------------------------------------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------------------------------------------
 */
/* Calculate the PV of the funding */
char *TargetNoteFunding(TARN_EVENT *event, double *res) {
  int l;
  double dFloatPV = 0.0;
  for (l = 0; l < event->num_fund_df; l++)
    dFloatPV +=
        event->fund_start[l] - event->fund_cvgMargin[l] * event->fund_end[l];
  res[0] = dFloatPV * TARGET_NOTE.m_dFundingNotional;

  return 0;
}
/* Calculate the PV of the funding */
/* -----------------------------------------------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------------------------------------------
 */
/* The coupon for the target note */
char *TargetNoteCoupon(TARN_EVENT *event,
                       double *res, // res[0] = fund  , res[1] = coupon
                       int *stop_path) {
  /* Variable declaration */
  double dCoupon, dNotionalRepay, dNextCumulCoupon, dFloater;
  char *err = 0;

  /* Calculate the coupon */
  dFloater = (event->fltr_start / event->fltr_end - 1.0) / event->fltr_cvg +
             event->fltr_spread;
  dCoupon = event->fixed_coupon + event->fltr_gearing * dFloater;

  /* Check to see if we need to floor or cap it */
  if (event->iIsFloored && dCoupon < event->dFloor)
    dCoupon = event->dFloor;
  if (event->iIsCapped && dCoupon > event->dCap)
    dCoupon = event->dFloor;

  /* Calculate the notional to be repaid */
  dNextCumulCoupon = TARGET_NOTE.m_dCumulCoupon + dCoupon * event->cumul_cvg;
  dNotionalRepay =
      1.0 - (TARGET_NOTE.m_dTarget - dNextCumulCoupon) / TARGET_NOTE.m_dSpread;
  dNotionalRepay = dNotionalRepay > 1.0 ? 1.0 : dNotionalRepay;
  dNotionalRepay = dNotionalRepay < TARGET_NOTE.m_dCouponNotionalRepaid
                       ? TARGET_NOTE.m_dCouponNotionalRepaid
                       : dNotionalRepay;

  /* Call the relevant coupon calculation function */
  switch (TARGET_NOTE.m_couponType) {
  case TN_ZC:
    err = TargetNoteCoupon_ZC(event, dCoupon, dNotionalRepay, res, stop_path);
    break;
  case TN_SWAP:
    err = TargetNoteCoupon_SWAP(event, dCoupon, dNotionalRepay, res, stop_path);
    break;
  case TN_FINAL:
    err = TargetNoteCoupon_FINAL(event, dCoupon, dFloater, dNotionalRepay, res,
                                 stop_path);
    break;
  case TN_SWAPKO:
    err =
        TargetNoteCoupon_SWAPKO(event, dCoupon, dNotionalRepay, res, stop_path);
    break;
  default:
    return "Unknown coupon type";
  }
  if (err)
    return err;

  if (*stop_path)
    return 0;

  /* Add the coupon to the cumulative value */
  TARGET_NOTE.m_dCumulCoupon += dCoupon * event->cumul_cvg;

  /* Return */
  return 0;
}
/* The coupon for the target note */
/* -----------------------------------------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------------------------------------
//

Payoff function for a ZC TargetNote

*/
char *TargetNoteCoupon_ZC(TARN_EVENT *event, double dCoupon,
                          double dNotionalRepay,
                          double *res, // res[0] = fund  , res[1] = coupon
                          int *stop_path) {
  /* Variable declaration */

  /* Initialize return values */
  res[0] = 0.0;
  res[1] = 0.0;
  *stop_path = 0;

  /* If it is the final period  , then pay the final coupon */
  if (event->iIsFinal && TARGET_NOTE.tarn->deal.bRepayCoupon) {
    res[1] = TARGET_NOTE.m_dCouponNotional * event->fixed_pay *
             TARGET_NOTE.m_dTarget;
    *stop_path = 1;
    res[0] += event->dKnockOutFee * event->fixed_pay;
    return 0;
  }

  /* if the target is exceeded pay the remaining coupon and cancel the deal */
  if (dNotionalRepay >= 1.0) {
    res[1] = event->fixed_pay * TARGET_NOTE.m_dTarget;
    res[1] *= TARGET_NOTE.m_dCouponNotional;
    *stop_path = 1;
    res[0] += event->dKnockOutFee * event->fixed_pay;
    return 0;
  }

  /* If the cumulative coupon exceeds the strike - spread * notional  ,
   * calculate the notional repayment */
  if (dNotionalRepay > TARGET_NOTE.m_dCouponNotionalRepaid) {
    double dIncPay = dNotionalRepay - TARGET_NOTE.m_dCouponNotionalRepaid;
    res[1] = dIncPay * TARGET_NOTE.m_dTarget * event->fixed_pay;
    TARGET_NOTE.m_dCouponNotionalRepaid = dNotionalRepay;
    TARGET_NOTE.m_dCouponNotional = 1.0 - dNotionalRepay;
    TARGET_NOTE.m_dFundingNotional = TARGET_NOTE.m_dCouponNotional;
    res[0] += dIncPay * event->dKnockOutFee * event->fixed_pay;
  }

  /* Return */
  return 0;
}

/* ------------------------------------------------------------------------------------------------------------------
//

Payoff function for the Swap

*/
char *TargetNoteCoupon_SWAP(TARN_EVENT *event, double dCoupon,
                            double dNotionalRepay,
                            double *res, // res[0] = fund  , res[1] = coupon
                            int *stop_path) {
  /* If it is the final period and the coupon is guaranteed  , then calculate
   * the coupon */
  if (event->iIsFinal && TARGET_NOTE.tarn->deal.bRepayCoupon) {
    res[1] = TARGET_NOTE.m_dCouponNotional * event->fixed_pay *
             /*event->fixed_cvg **/
             (TARGET_NOTE.m_dTarget - TARGET_NOTE.m_dCumulCoupon);
    *stop_path = 1;
    return 0;
  }

  /* if the target is exceeded pay the remaining coupon and cancel the deal */
  if (dNotionalRepay >= 1.0) {
    res[1] = event->fixed_pay /* event->fixed_cvg */
             * (TARGET_NOTE.m_dTarget - TARGET_NOTE.m_dCumulCoupon);
    res[1] *= TARGET_NOTE.m_dCouponNotional;
    res[0] += event->dKnockOutFee * event->fixed_pay;
    *stop_path = 1;
    return 0;
  }

  /* If the cumulative coupon exceeds the strike - spread * notional  ,
   * calculate the notional repayment */
  if (dNotionalRepay > TARGET_NOTE.m_dCouponNotionalRepaid) {
    double dIncPay = dNotionalRepay - TARGET_NOTE.m_dCouponNotionalRepaid;
    res[1] = dIncPay * /* event->fixed_cvg * */
             (TARGET_NOTE.m_dTarget - TARGET_NOTE.m_dCumulCoupon) *
             event->fixed_pay;
    TARGET_NOTE.m_dCouponNotionalRepaid = dNotionalRepay;
    TARGET_NOTE.m_dCouponNotional = 1.0 - dNotionalRepay;
    TARGET_NOTE.m_dFundingNotional = TARGET_NOTE.m_dCouponNotional;
    res[0] += dIncPay * event->dKnockOutFee * event->fixed_pay;
  }

  /* Calculate any coupon payments for the TN_SWAP case */
  res[1] += TARGET_NOTE.m_dCouponNotional * event->fixed_pay *
            event->fixed_cvg * dCoupon;

  /* Return */
  return 0;
}

/* ------------------------------------------------------------------------------------------------------------------
//

Payoff function for the	Knock-out Swap

*/
char *TargetNoteCoupon_SWAPKO(TARN_EVENT *event, double dCoupon,
                              double dNotionalRepay,
                              double *res, // res[0] = fund  , res[1] = coupon
                              int *stop_path) {
  /* if the target is exceeded pay the remaining coupon and cancel the deal */
  if (dNotionalRepay >= 1.0) {
    res[1] = event->fixed_pay * event->fixed_cvg * dCoupon *
             TARGET_NOTE.m_dCouponNotional;
    res[0] += event->dKnockOutFee * event->fixed_pay;
    *stop_path = 1;
    return 0;
  }

  /* Calculate any coupon payments */
  res[1] += TARGET_NOTE.m_dCouponNotional * event->fixed_pay *
            event->fixed_cvg * dCoupon;

  /* If the cumulative coupon exceeds the strike - spread * notional  ,
   * calculate the notional repayment */
  if (dNotionalRepay > TARGET_NOTE.m_dCouponNotionalRepaid) {
    double dIncPay = dNotionalRepay - TARGET_NOTE.m_dCouponNotionalRepaid;
    TARGET_NOTE.m_dCouponNotionalRepaid = dNotionalRepay;
    TARGET_NOTE.m_dCouponNotional = 1.0 - dNotionalRepay;
    TARGET_NOTE.m_dFundingNotional = TARGET_NOTE.m_dCouponNotional;
    res[0] += dIncPay * event->dKnockOutFee * event->fixed_pay;
  }

  /* Return */
  return 0;
}

/* ------------------------------------------------------------------------------------------------------------------
//

Payoff function for the accrued final payment coupon

*/
char *TargetNoteCoupon_FINAL(TARN_EVENT *event, double dCoupon, double dFloater,
                             double dNotionalRepay,
                             double *res, // res[0] = fund  , res[1] = coupon
                             int *stop_path) {
  /* Variable declaration */
  char *err = 0;

  /* Accrue any existing coupon value by the preivous fltr setting */
  TARGET_NOTE.m_dFinalCoupon *= TARGET_NOTE.m_dPrevAccrual;

  /* If it is the final period  , then calculate the coupon */
  if (event->iIsFinal) {
    double dIncPay;
    dIncPay = 1.0 - TARGET_NOTE.m_dCouponNotionalRepaid;
    TARGET_NOTE.m_dFinalCoupon += dIncPay * (1.0 + TARGET_NOTE.m_dTarget);
    res[1] = (TARGET_NOTE.m_dFinalCoupon - 1.0) * event->fixed_pay;
    *stop_path = 1;
    res[0] += event->dKnockOutFee * event->fixed_pay;
    return 0;
  }

  /* If the cumulative coupon exceeds the strike - spread * notional  ,
   * calculate the notional repayment */
  if (dNotionalRepay > TARGET_NOTE.m_dCouponNotionalRepaid) {
    double dIncPay = dNotionalRepay - TARGET_NOTE.m_dCouponNotionalRepaid;
    TARGET_NOTE.m_dFinalCoupon += dIncPay * (1.0 + TARGET_NOTE.m_dTarget);
    TARGET_NOTE.m_dCouponNotionalRepaid = dNotionalRepay;
    TARGET_NOTE.m_dCouponNotional = 1.0 - dNotionalRepay;
    res[0] += dIncPay * event->dKnockOutFee * event->fixed_pay;
  }

  /* Record the floater * coverage for possible use in next coupon */
  TARGET_NOTE.m_dPrevAccrual = 1.0 + dFloater * event->fixed_cvg;

  /* Return */
  return 0;
}
