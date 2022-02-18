
#include "CRIFCaller.h"

#include "CRIF.h"
#include "CRIFProdStruct.h"
#include "CRIFUtil.h"

// ----------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------
//                                 Calling routine for Callable Ratchet Inverse Floater (CRIF)
// ----------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------

Err CRIFAutocal(
    CMSCTS_MARKET       sMarket,
    CRIF_DEAL           sDeal,
    CMSCTS_CALIB        sCalibration,
    CRIF_PRICING_PARAMS sPricingParams,
    CRIF_OUTPUTS        sOutputs)
{
    Err           err = NULL;
    double        dFunding, dExotic;
    CRIF_SIMULARG sSimulArg = NULL;
    int           i;
    clock_t       time1, time2;

    time1 = clock();

    /* Fill the structures */
    err = crif_fill_check_all_struct(sDeal, sMarket, sPricingParams->sCouponPricingParams);

    if (err)
        goto FREE_RETURN;

    /* Calculate today's coupon values (all Funding + Historical Ratchet Coupons*/
    err = crif_calc_mkt_value(sDeal, sMarket, sPricingParams, sCalibration);
    if (err)
        goto FREE_RETURN;

    time2 = clock();
    smessage("CRIF IV Pricing: %.2f", (double)(time2 - time1) / CLOCKS_PER_SEC);

    /* fill IV's informations */
    dFunding = 0.0;

    for (i = sDeal->sFunding->sAux->iStartIndex; i <= sDeal->sFunding->sAux->iEndIndex; i++)
    {
        dFunding += sDeal->sFunding->sAux->dMarketValue[i];
    }

    sOutputs->dFunding = dFunding;

    dExotic = 0.0;

    for (i = sDeal->sExotic->sAux->iStartIndex; i <= sDeal->sExotic->sAux->iEndIndex; i++)
    {
        dExotic += sDeal->sExotic->sAux->dMarketValue[i];
    }

    sOutputs->dExotic = dExotic;

    /* Calculate the Option */

    if (sDeal->sExotic->sAux->iNbUsedNumCpn > 0)
    {
        /* Allocation */
        sSimulArg = calloc(1, sizeof(crif_simularg));

        if (!sSimulArg)
        {
            err = "Memory allocation faillure in crifAutocal";
            goto FREE_RETURN;
        }

        err = crif_allocate_simularg(sSimulArg);
        if (err)
            goto FREE_RETURN;

        /* First Calibrate the Model */
        err = crif_calibrate_model(sMarket, sDeal, sSimulArg, sCalibration, sOutputs);
        if (err)
            goto FREE_RETURN;

        time1 = clock();
        smessage("CRIF Calibration: %.2f", (double)(time1 - time2) / CLOCKS_PER_SEC);

        /* calculate the CV Exotic coupon value */
        err = crif_calc_cv_value(sMarket, sDeal, sSimulArg, sPricingParams, sOutputs);
        if (err)
            goto FREE_RETURN;

        /* Calculates the unadjusted option */
        err = crif_calculate_option(sMarket, sDeal, sSimulArg, sPricingParams, sOutputs);
        if (err)
            goto FREE_RETURN;

        time2 = clock();
        smessage("CRIF Option Pricing: %.2f", (double)(time2 - time1) / CLOCKS_PER_SEC);

        /* Calculates the switch adjsutment */

        // not yet implemented

        /* Final Call */
        sOutputs->dCall = sOutputs->dCallInit;
    }
    else
    {
        sOutputs->dCall = 0.0;
    }

    /* Fill the Outputs */
    err = crif_fill_outputs(sDeal, sSimulArg, sOutputs);

    if (err)
        goto FREE_RETURN;

FREE_RETURN:

    if (sSimulArg)
    {
        crif_free_simularg(sSimulArg);
        free(sSimulArg);
    }

    return err;
}

void crif_deal_set_default(CRIF_DEAL sDeal)
{
    sDeal->dPayRec             = -1.0;
    sDeal->sCall->dPayRec      = 1.0;
    sDeal->sCall->iIsCallable  = 1;
    sDeal->sExotic->iFixingLag = 2;
}

void crif_calibration_set_default(CMSCTS_CALIB sCalibration)
{
    sCalibration->iUseCalib = 1;
    //	sCalibration->cUndName = "CAL";

    sCalibration->iUseSV = 1;

    sCalibration->iNbFactor = 2;
    sCalibration->dLambda   = 0.1;
    sCalibration->dLGMAlpha = 1.2;
    sCalibration->dLGMGamma = 0.2;
    sCalibration->dLGMRho   = -0.85;

    sCalibration->iNbTimesSV     = 0;
    sCalibration->dTimesSV       = NULL;
    sCalibration->dAlphaSV       = NULL;
    sCalibration->dRho2SV        = NULL;
    sCalibration->dRhoSV         = NULL;
    sCalibration->dCalibSmileSTD = 1.0;
    sCalibration->dTStar         = 10.0;

    sCalibration->eCalibStrat = CMSCTS_DIAGCMS1;
    strcpy(sCalibration->cMaxCalibTenor, "30Y");
    sCalibration->iUseInArrearsFixing = 1;

    sCalibration->iFixLambda = 0;
    cpd_calib_set_default_param(sCalibration->sPrimCalibParams);
    cpd_calib_set_default_param(sCalibration->sSecCalibParams);
    diag_calib_lm_params_set_default_param(sCalibration->sLMParams);
    LGMSV_SetDefault_CalibParams(sCalibration->sLGMSVCalibParams);
    lgmsv_app_set_default_params_struct(sCalibration->sLGMSVNumerParams);
}

void crif_pricingparams_set_default(CRIF_PRICING_PARAMS sPricingParams)
{
    CRIF_CPN_NUMPARAMS    sCouponPricingParams;
    CMSCTS_SIMUL_PARAMS   sSimulParams;
    CMSCTS_RESERVE_PARAMS sReserveParams;
    Err                   err;

    sCouponPricingParams = sPricingParams->sCouponPricingParams;

    /* Coupon Params */
    sCouponPricingParams->dMCEBCallSpread = 0.0;

    /* Simul Params */
    sSimulParams = sPricingParams->sSimulParams;

    sSimulParams->dMCMaxTime     = 10.0 / 300.0;
    sSimulParams->iAdjustIV      = 0;
    sSimulParams->iNumeraireType = 0;
    sSimulParams->lNbPaths       = 20000;
    sSimulParams->lNbSteps       = 300;

    err = Fill_lgmSV_defaultParam(sSimulParams->sLGMSVParams);

    /* Reserve Params */
    sReserveParams = sPricingParams->sReserveParams;

    sReserveParams->iCalcSwitchAdjust = 1;
    sReserveParams->iGMANbFactor      = 1;
    sReserveParams->dGMAAlpha         = 0.01;
    sReserveParams->dGMAGamma         = 0.01;
    sReserveParams->dGMARho           = 0.0;

    sReserveParams->eCalibStrat = CMSCTS_DIAGCMS1;
    strcpy(sReserveParams->cLGMModelName, "CRIFAUTOCALRESUND");
    sReserveParams->iLGMNbFactor    = 2;
    sReserveParams->dLGMLambda      = 0.1;
    sReserveParams->iLGMFixLambda   = 0;
    sReserveParams->dLGMAlpha       = 1.2;
    sReserveParams->dLGMGamma       = 0.2;
    sReserveParams->dLGMRho         = -0.85;
    sReserveParams->dLGMLambdaShift = 0.0;
}

void crif_outputs_set_default(CRIF_OUTPUTS sOutputs)
{
    sOutputs->dCall    = 0.0;
    sOutputs->dExotic  = 0.0;
    sOutputs->dFunding = 0.0;
    sOutputs->dStd     = 0.0;

    sOutputs->iCalcExeProbas = 0;
}

Err crif_allocate_outputs(CRIF_OUTPUTS sOutputs)
{
    sOutputs->sInstDatas       = calloc(1, sizeof(cpd_calib_inst_data));
    sOutputs->sLGMSVModel      = calloc(1, sizeof(LGMSV_model));
    sOutputs->sFwdIVInfos      = calloc(1, sizeof(crif_fwdiv));
    sOutputs->sSwitchInstDatas = calloc(1, sizeof(cpd_calib_inst_data));

    if (!sOutputs->sInstDatas || !sOutputs->sLGMSVModel || !sOutputs->sFwdIVInfos ||
        !sOutputs->sSwitchInstDatas)
    {
        return "Memory allcation faillure in crif_allocate_outputs";
    }

    return NULL;
}

void crif_free_outputs(CRIF_OUTPUTS sOutputs)
{
    if (sOutputs->sInstDatas)
    {
        cpd_free_calib_inst_data(sOutputs->sInstDatas);
        free(sOutputs->sInstDatas);
    }

    if (sOutputs->sLGMSVModel)
    {
        free_LGMSV_model(sOutputs->sLGMSVModel);
        free(sOutputs->sLGMSVModel);
    }

    if (sOutputs->sFwdIVInfos)
    {
        crif_free_fwdiv(sOutputs->sFwdIVInfos);
    }

    if (sOutputs->sSwitchInstDatas)
    {
        cpd_free_calib_inst_data(sOutputs->sSwitchInstDatas);
        free(sOutputs->sSwitchInstDatas);
    }
}

Err crif_fill_outputs(CRIF_DEAL sDeal, CRIF_SIMULARG sSimulArg, CRIF_OUTPUTS sOutputs)
{
    Err           err = NULL;
    LGMSV_MODEL   sModel;
    CRIF_CALL     sCall;
    CRIF_CALL_AUX sCallAux;
    CRIF_MCARG    sMCArg;
    CRIF_PAYARG   sPayArg;
    CRIF_FWDIV    sFwdIVInfos;

    int i, iIndexStart, iIndexCall, iIndexEvt;

    /* Fill Model */
    if (sDeal->sExotic->sAux->iNbUsedNumCpn > 0)
    {
        sCall       = sDeal->sCall;
        sCallAux    = sCall->sAux;
        sModel      = sSimulArg->sModel;
        sMCArg      = sSimulArg->sMCArg;
        sFwdIVInfos = sOutputs->sFwdIVInfos;

        if (sOutputs->sLGMSVModel)
        {
            err = init_LGMSV_model(
                sOutputs->sLGMSVModel,
                sModel->lToday,
                sModel->iOne2F,
                sModel->iNbPWTime,
                sModel->dLambdaX,
                sModel->dPWTime,
                sModel->dSigma,
                sModel->dAlpha,
                sModel->dLambdaEps,
                sModel->dLvlEps,
                sModel->dRho,
                sModel->dTStar,
                sModel->dInitLGMAlpha,
                sModel->dLGMGamma,
                sModel->dInitLGMRho,
                sModel->dRho2);
            if (err)
                goto FREE_RETURN;

            Convert_Tstar_model(sOutputs->sLGMSVModel, sModel->dInitTStar);

            ConvertTS_LGMSV_to_LGM(
                sOutputs->sLGMSVModel->iNbPWTime,
                sOutputs->sLGMSVModel->dPWTime,
                sOutputs->sLGMSVModel->dSigma,
                sOutputs->sLGMSVModel->dLambdaX,
                sOutputs->sLGMSVModel->dTStar);

            for (i = 0; i < sOutputs->sLGMSVModel->iNbPWTime; i++)
            {
                sOutputs->sLGMSVModel->dAlpha[i] /= 2.0;
                sOutputs->sLGMSVModel->dLambdaEps[i] /= 2.0;
            }
        }

        /* Fwd IVs */
        if (sCallAux->iNbUsedCall > 0 && sFwdIVInfos)
        {
            iIndexStart = sCallAux->iStartIndex;

            err = crif_allocate_fwdiv(sFwdIVInfos, sCallAux->iNbUsedCall);
            if (err)
                goto FREE_RETURN;

            /* Go to the first call */
            iIndexEvt = 0;
            sPayArg   = (CRIF_PAYARG)(sMCArg->vPayoffParams[sMCArg->iIndexEvent[iIndexEvt]]);

            while (sPayArg && !sPayArg->iIsCall)
            {
                iIndexEvt++;
                if (iIndexEvt < sMCArg->iNbEvent)
                {
                    sPayArg = (CRIF_PAYARG)(sMCArg->vPayoffParams[sMCArg->iIndexEvent[iIndexEvt]]);
                }
                else
                {
                    sPayArg = NULL;
                }
            }

            for (i = 0; i < sFwdIVInfos->lNbCall; i++)
            {
                iIndexCall                = iIndexStart + i;
                sFwdIVInfos->lExeDates[i] = sCall->lExeDate[iIndexCall];

                /* Model IV */
                sFwdIVInfos->dMarketIV[i]   = sMCArg->sMCEBParams->dMarketFwdIV[iIndexEvt];
                sFwdIVInfos->dModelIV[i]    = sMCArg->sMCEBParams->dModelFwdIV[iIndexEvt];
                sFwdIVInfos->dFees[i]       = sMCArg->sMCEBParams->dFee[iIndexEvt];
                sFwdIVInfos->dNewModelIV[i] = sMCArg->dModelValue[iIndexEvt];

                if (sMCArg->sMCEBParams->dExeProba)
                    sFwdIVInfos->dExeProba[i] = sMCArg->sMCEBParams->dExeProba[iIndexEvt];
                if (sMCArg->sMCEBParams->dOneTimeCall)
                    sFwdIVInfos->dOTCall[i] = sMCArg->sMCEBParams->dOneTimeCall[iIndexEvt];
                if (sMCArg->sMCEBParams->dOneTimePartial)
                    sFwdIVInfos->dOTCallPartial[i] =
                        sMCArg->sMCEBParams->dOneTimePartial[iIndexEvt];

                iIndexEvt++;

                if (iIndexEvt < sMCArg->iNbEvent)
                {
                    sPayArg = (CRIF_PAYARG)(sMCArg->vPayoffParams[sMCArg->iIndexEvent[iIndexEvt]]);
                }
                else
                {
                    sPayArg = NULL;
                }

                while (sPayArg && !sPayArg->iIsCall)
                {
                    sFwdIVInfos->dMarketIV[i] += sMCArg->sMCEBParams->dMarketFwdIV[iIndexEvt];
                    sFwdIVInfos->dModelIV[i] += sMCArg->sMCEBParams->dModelFwdIV[iIndexEvt];
                    sFwdIVInfos->dFees[i] += sMCArg->sMCEBParams->dFee[iIndexEvt];
                    sFwdIVInfos->dNewModelIV[i] += sMCArg->dModelValue[iIndexEvt];

                    iIndexEvt++;
                    if (iIndexEvt < sMCArg->iNbEvent)
                    {
                        sPayArg =
                            (CRIF_PAYARG)(sMCArg->vPayoffParams[sMCArg->iIndexEvent[iIndexEvt]]);
                    }
                    else
                    {
                        sPayArg = NULL;
                    }
                }
            }
        }
    }

FREE_RETURN:

    return err;
}

Err crif_allocate_fwdiv(CRIF_FWDIV sFwdIVInfos, long lNbCall)
{
    sFwdIVInfos->lNbCall = lNbCall;

    sFwdIVInfos->lExeDates      = calloc(sFwdIVInfos->lNbCall, sizeof(long));
    sFwdIVInfos->dMarketIV      = calloc(sFwdIVInfos->lNbCall, sizeof(double));
    sFwdIVInfos->dModelIV       = calloc(sFwdIVInfos->lNbCall, sizeof(double));
    sFwdIVInfos->dNewModelIV    = calloc(sFwdIVInfos->lNbCall, sizeof(double));
    sFwdIVInfos->dExeProba      = calloc(sFwdIVInfos->lNbCall, sizeof(double));
    sFwdIVInfos->dFees          = calloc(sFwdIVInfos->lNbCall, sizeof(double));
    sFwdIVInfos->dOTCall        = calloc(sFwdIVInfos->lNbCall, sizeof(double));
    sFwdIVInfos->dOTCallPartial = calloc(sFwdIVInfos->lNbCall, sizeof(double));

    if (!sFwdIVInfos->lExeDates || !sFwdIVInfos->dMarketIV || !sFwdIVInfos->dModelIV ||
        !sFwdIVInfos->dExeProba || !sFwdIVInfos->dNewModelIV || !sFwdIVInfos->dFees ||
        !sFwdIVInfos->dOTCall || !sFwdIVInfos->dOTCallPartial)
    {
        return "Memory allocation faillure in crif_fill_outputs";
    }

    return NULL;
}

void crif_free_fwdiv(CRIF_FWDIV sFwdIVInfos)
{
    if (sFwdIVInfos->lExeDates)
        free(sFwdIVInfos->lExeDates);
    if (sFwdIVInfos->dMarketIV)
        free(sFwdIVInfos->dMarketIV);
    if (sFwdIVInfos->dModelIV)
        free(sFwdIVInfos->dModelIV);
    if (sFwdIVInfos->dNewModelIV)
        free(sFwdIVInfos->dNewModelIV);
    if (sFwdIVInfos->dExeProba)
        free(sFwdIVInfos->dExeProba);
    if (sFwdIVInfos->dFees)
        free(sFwdIVInfos->dFees);
    if (sFwdIVInfos->dOTCall)
        free(sFwdIVInfos->dOTCall);
    if (sFwdIVInfos->dOTCallPartial)
        free(sFwdIVInfos->dOTCallPartial);

    sFwdIVInfos->lExeDates      = NULL;
    sFwdIVInfos->dMarketIV      = NULL;
    sFwdIVInfos->dModelIV       = NULL;
    sFwdIVInfos->dNewModelIV    = NULL;
    sFwdIVInfos->dExeProba      = NULL;
    sFwdIVInfos->dFees          = NULL;
    sFwdIVInfos->dOTCall        = NULL;
    sFwdIVInfos->dOTCallPartial = NULL;
}
