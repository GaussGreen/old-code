#include "CMSCTSCaller.h"

#include "CMSCTSOption.h"
#include "CMSCTSPricing.h"

Err CMSCTSAutocal(
    CMSCTS_MARKET         sMarket,
    CMSCTS_DEAL           sDeal,
    CMSCTS_CALIB          sCalibration,
    CMSCTS_PRICING_PARAMS sPricingParams,
    CMSCTS_OUTPUTS        sOutputs)
{
    Err             err = NULL;
    double          dFunding, dExotic;
    CMSCTS_SIMULARG sSimulArg = NULL;
    int             i;
    clock_t         time1, time2;

    time1 = clock();

    /* Fill the structures */
    err = cmscts_fill_check_all_struct(sDeal, sMarket, sPricingParams->sCouponPricingParams);

    if (err)
        goto FREE_RETURN;

    /* Calculate today's coupon values */
    err = cmscts_calc_mkt_value(sDeal, sMarket, sPricingParams, sCalibration);

    if (err)
        goto FREE_RETURN;

    time2 = clock();
    smessage("CMSCTS IV Pricing: %.2f", (double)(time2 - time1) / CLOCKS_PER_SEC);

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
        sSimulArg = (cmscts_simularg*)calloc(1, sizeof(cmscts_simularg));

        if (!sSimulArg)
        {
            err = "Memory allocation faillure in CMSCTSAutocal";
            goto FREE_RETURN;
        }

        err = cmscts_allocate_simularg(sSimulArg);

        if (err)
            goto FREE_RETURN;

        /* First Calibrate the Model */
        err = cmscts_calibrate_model(sMarket, sDeal, sSimulArg, sCalibration, sOutputs);

        if (err)
            goto FREE_RETURN;

        time1 = clock();
        smessage("CMSCTS Calibration: %.2f", (double)(time1 - time2) / CLOCKS_PER_SEC);

        /* Calculates the unadjusted option */
        err = cmscts_calculate_option(sMarket, sDeal, sSimulArg, sPricingParams, sOutputs);

        if (err)
            goto FREE_RETURN;

        time2 = clock();
        smessage("CMSCTS Option Pricing: %.2f", (double)(time2 - time1) / CLOCKS_PER_SEC);

        /* Calculates the switch adjsutment */
        if (sDeal->sCall->sAux->iNbUsedCall > 1)
        {
            err = cmscts_calculate_switch_adjust(
                sMarket, sDeal, sSimulArg, sCalibration, sPricingParams, sOutputs);

            if (err)
                goto FREE_RETURN;

            time1 = clock();
            smessage("CMSCTS Switch Adjustment: %.2f", (double)(time1 - time2) / CLOCKS_PER_SEC);
        }
        else
        {
            sOutputs->dSwitchAdjust = 0.0;
        }

        /* Final Call */
        sOutputs->dCall = sOutputs->dCallInit + sOutputs->dSwitchAdjust;
    }
    else
    {
        sOutputs->dCall = 0.0;
    }

    /* Fill the Outputs */
    err = cmscts_fill_outputs(sDeal, sSimulArg, sOutputs);

    if (err)
        goto FREE_RETURN;

FREE_RETURN:

    if (sSimulArg)
    {
        cmscts_free_simularg(sSimulArg);
        free(sSimulArg);
    }

    return err;
}

void cmscts_market_set_default(CMSCTS_MARKET sMarket)
{
    sMarket->cFixedBasis[0]  = '\0';
    sMarket->cFixedFreq[0]   = '\0';
    sMarket->cFloatBasis[0]  = '\0';
    sMarket->cFundYcName[0]  = '\0';
    sMarket->cRefRateName[0] = '\0';
    sMarket->cVcName[0]      = '\0';
    sMarket->cYcName[0]      = '\0';

    sMarket->dFundSpotFx   = 0.0;
    sMarket->dStrikesInVol = NULL;
    sMarket->eVolType      = (SrtDiffusionType)0;
    sMarket->GetCashVol    = NULL;
    sMarket->iCashVol      = 0;
    sMarket->iEODExeFlag   = 0;
    sMarket->iEODFixFlag   = 0;
    sMarket->iEODPayFlag   = 0;
    sMarket->iIsSABR       = 0;
    //	sMarket->iIsSABRAF = 0;
    sMarket->lFundFxSpotDate  = 0;
    sMarket->lNumStrikesInVol = 0;
    sMarket->lToday           = 0;
}

void cmscts_deal_set_default(CMSCTS_DEAL sDeal)
{
    sDeal->sExotic->iAccuralFixLag = 2;
    strcpy(sDeal->sExotic->cAccrualFixFreq, "1D");
    sDeal->dPayRec             = -1.0;
    sDeal->sCall->dPayRec      = 1.0;
    sDeal->sExotic->iFixingLag = 2;
    sDeal->sKO->iFixingLag     = 2;
}

void cmscts_calibration_set_default(CMSCTS_CALIB sCalibration)
{
    sCalibration->iUseCalib = 1;
    sCalibration->cUndName  = NULL;

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

    sCalibration->eCalibStrat = CMSCTS_LONGSHORT;
    strcpy(sCalibration->cMaxCalibTenor, "30Y");
    sCalibration->iUseInArrearsFixing = 1;

    sCalibration->iFixLambda = 0;
    cpd_calib_set_default_param(sCalibration->sPrimCalibParams);
    cpd_calib_set_default_param(sCalibration->sSecCalibParams);
    diag_calib_lm_params_set_default_param(sCalibration->sLMParams);
    LGMSV_SetDefault_CalibParams(sCalibration->sLGMSVCalibParams);
    lgmsv_app_set_default_params_struct(sCalibration->sLGMSVNumerParams);
}

void cmscts_pricingparams_set_default(CMSCTS_PRICING_PARAMS sPricingParams)
{
    CMSCTS_CPN_NUMPARAMS  sCouponPricingParams;
    CMSCTS_SIMUL_PARAMS   sSimulParams;
    CMSCTS_RESERVE_PARAMS sReserveParams;
    Err                   err;

    sCouponPricingParams = sPricingParams->sCouponPricingParams;

    /* Coupon Params */
    sCouponPricingParams->iIVBuySell    = -1;
    sCouponPricingParams->dIVCallSpread = 10.0 / 10000;
    sCouponPricingParams->iIVTrimType   = 2;
    sCouponPricingParams->dIVMinFixTime = 1.0 / 16.0;
    sCouponPricingParams->iIVMaxFix     = 10000;

    sCouponPricingParams->iAccrueOnBarrier = 1;

    sCouponPricingParams->iCallBuySell    = -1;
    sCouponPricingParams->dCallCallSpread = 10.0 / 10000;
    sCouponPricingParams->iCallTrimType   = 1;
    sCouponPricingParams->dCallMinFixTime = 0.25;
    sCouponPricingParams->iCallMaxFix     = 2;

    sCouponPricingParams->iValueZero  = 1;
    sCouponPricingParams->dMinBarrier = -0.5;
    sCouponPricingParams->dMaxBarrier = 0.5;

    sCouponPricingParams->dMCEBCallSpread = 0.0;
    sCouponPricingParams->iKOBuySell      = -1;
    sCouponPricingParams->dKOCallSpread   = 10.0 / 10000;

    sCouponPricingParams->iUseLGMCorrel   = 0;
    sCouponPricingParams->dLGMCorrelAlpha = 0.0;
    sCouponPricingParams->dLGMCorrelGamma = 0.0;
    sCouponPricingParams->dLGMCorrelRho   = 0.0;

    cms_spread_set_default_model_params(sCouponPricingParams->sCMSSpreadModel);

    /* Simul Params */
    sSimulParams = sPricingParams->sSimulParams;

    sSimulParams->dMCMaxTime        = 10.0 / 300.0;
    sSimulParams->iAdjustIV         = 1;
    sSimulParams->iEstimateCallOnly = 0;
    sSimulParams->iNumeraireType    = 1;
    sSimulParams->lNbPaths          = 20000;
    sSimulParams->lNbSteps          = 300;
    sSimulParams->iTrimFreq         = 1;

    err = Fill_lgmSV_defaultParam(sSimulParams->sLGMSVParams);

    /* Reserve Params */
    sReserveParams = sPricingParams->sReserveParams;

    sReserveParams->iCalcSwitchAdjust = 1;
    sReserveParams->iGMANbFactor      = 1;
    sReserveParams->dGMAAlpha         = 0.01;
    sReserveParams->dGMAGamma         = 0.01;
    sReserveParams->dGMARho           = 0.0;

    sReserveParams->eCalibStrat = CMSCTS_DIAGCAP;
    strcpy(sReserveParams->cLGMModelName, "CMSCTSAUTOCALRESUND");
    sReserveParams->iLGMNbFactor    = 2;
    sReserveParams->dLGMLambda      = 0.1;
    sReserveParams->iLGMFixLambda   = 0;
    sReserveParams->dLGMAlpha       = 1.2;
    sReserveParams->dLGMGamma       = 0.2;
    sReserveParams->dLGMRho         = -0.85;
    sReserveParams->dLGMLambdaShift = 0.0;
}

void cmscts_outputs_set_default(CMSCTS_OUTPUTS sOutputs)
{
    sOutputs->dCall    = 0.0;
    sOutputs->dExotic  = 0.0;
    sOutputs->dFunding = 0.0;
    sOutputs->dStd     = 0.0;

    sOutputs->iCalcExeProbas = 0;
}

Err cmscts_allocate_outputs(CMSCTS_OUTPUTS sOutputs)
{
    sOutputs->sInstDatas       = (cpd_calib_inst_data*)calloc(1, sizeof(cpd_calib_inst_data));
    sOutputs->sLGMSVModel      = (LGMSV_model*)calloc(1, sizeof(LGMSV_model));
    sOutputs->sFwdIVInfos      = (cmscts_fwdiv*)calloc(1, sizeof(cmscts_fwdiv));
    sOutputs->sSwitchInstDatas = (cpd_calib_inst_data*)calloc(1, sizeof(cpd_calib_inst_data));

    if (!sOutputs->sInstDatas || !sOutputs->sLGMSVModel || !sOutputs->sFwdIVInfos ||
        !sOutputs->sSwitchInstDatas)
    {
        return "Memory allcation faillure in cmscts_allocate_outputs";
    }

    return NULL;
}

void cmscts_free_outputs(CMSCTS_OUTPUTS sOutputs)
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
        cmscts_free_fwdiv(sOutputs->sFwdIVInfos);
    }

    if (sOutputs->sSwitchInstDatas)
    {
        cpd_free_calib_inst_data(sOutputs->sSwitchInstDatas);
        free(sOutputs->sSwitchInstDatas);
    }
}

Err cmscts_fill_outputs(CMSCTS_DEAL sDeal, CMSCTS_SIMULARG sSimulArg, CMSCTS_OUTPUTS sOutputs)
{
    Err             err = NULL;
    LGMSV_MODEL     sModel;
    CMSCTS_CALL     sCall;
    CMSCTS_CALL_AUX sCallAux;
    CMSCTS_KO       sKO;
    CMSCTS_KO_AUX   sKOAux;
    CMSCTS_MCARG    sMCArg;
    CMSCTS_PAYARG   sPayArg;
    CMSCTS_FWDIV    sFwdIVInfos;
    int             i, iIndexStart, iIndexCall, iIndexKO, iIndexEvt;
    int             iIncludePartialPV;

    /* Fill Model */
    if (sDeal->sExotic->sAux->iNbUsedNumCpn > 0)
    {
        sCall       = sDeal->sCall;
        sCallAux    = sCall->sAux;
        sKO         = sDeal->sKO;
        sKOAux      = sKO->sAux;
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

            err = cnscts_allocate_fwdiv(sFwdIVInfos, sCallAux->iNbUsedCall);

            if (err)
                goto FREE_RETURN;

            /* Go to the first call */
            iIndexEvt = 0;
            sPayArg   = (CMSCTS_PAYARG)(sMCArg->vPayoffParams[sMCArg->iIndexEvent[iIndexEvt]]);

            while (sPayArg && !sPayArg->iIsCall)
            {
                iIndexEvt++;
                if (iIndexEvt < sMCArg->iNbEvent)
                {
                    sPayArg =
                        (CMSCTS_PAYARG)(sMCArg->vPayoffParams[sMCArg->iIndexEvent[iIndexEvt]]);
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

                if (sMCArg->sMCEBParams->dOneTimeCall)
                    sFwdIVInfos->dOTCall[i] = sMCArg->sMCEBParams->dOneTimeCall[iIndexEvt];
                if (sMCArg->sMCEBParams->dOneTimePartial)
                    sFwdIVInfos->dOTCallPartial[i] =
                        sMCArg->sMCEBParams->dOneTimePartial[iIndexEvt];

                iIndexEvt++;

                if (iIndexEvt < sMCArg->iNbEvent)
                {
                    sPayArg =
                        (CMSCTS_PAYARG)(sMCArg->vPayoffParams[sMCArg->iIndexEvent[iIndexEvt]]);
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
                            (CMSCTS_PAYARG)(sMCArg->vPayoffParams[sMCArg->iIndexEvent[iIndexEvt]]);
                    }
                    else
                    {
                        sPayArg = NULL;
                    }
                }
            }
        }
        else if (sKOAux->iNbUsedKO > 0 && sFwdIVInfos)
        {
            iIndexEvt = 0;
            sPayArg   = (CMSCTS_PAYARG)(sMCArg->vPayoffParams[sMCArg->iIndexEvent[iIndexEvt]]);

            iIndexStart       = sKOAux->iStartIndex;
            iIncludePartialPV = 0;

            if (!sPayArg->iIsKO && iIndexStart > 0)
            {
                /* we include first payment */
                iIndexStart--;

                err = cnscts_allocate_fwdiv(sFwdIVInfos, sKOAux->iNbUsedKO + 1);

                if (err)
                    goto FREE_RETURN;

                iIncludePartialPV = 1;
            }
            else
            {
                err = cnscts_allocate_fwdiv(sFwdIVInfos, sKOAux->iNbUsedKO);

                if (err)
                    goto FREE_RETURN;
            }

            for (i = 0; i < sFwdIVInfos->lNbCall; i++)
            {
                iIndexKO                  = iIndexStart + i;
                sFwdIVInfos->lExeDates[i] = sKO->lFixingDate[iIndexKO];

                /* Model IV */
                sFwdIVInfos->dMarketIV[i]   = sMCArg->sMCEBParams->dMarketFwdIV[iIndexEvt];
                sFwdIVInfos->dModelIV[i]    = sMCArg->sMCEBParams->dModelFwdIV[iIndexEvt];
                sFwdIVInfos->dFees[i]       = sMCArg->sMCEBParams->dFee[iIndexEvt];
                sFwdIVInfos->dNewModelIV[i] = sMCArg->dModelValue[iIndexEvt];

                iIndexEvt++;

                if (iIndexEvt < sMCArg->iNbEvent)
                {
                    sPayArg =
                        (CMSCTS_PAYARG)(sMCArg->vPayoffParams[sMCArg->iIndexEvent[iIndexEvt]]);
                }
                else
                {
                    sPayArg = NULL;
                }

                while (sPayArg && !sPayArg->iIsKO)
                {
                    sFwdIVInfos->dMarketIV[i] += sMCArg->sMCEBParams->dMarketFwdIV[iIndexEvt];
                    sFwdIVInfos->dModelIV[i] += sMCArg->sMCEBParams->dModelFwdIV[iIndexEvt];
                    sFwdIVInfos->dFees[i] += sMCArg->sMCEBParams->dFee[iIndexEvt];
                    sFwdIVInfos->dNewModelIV[i] += sMCArg->dModelValue[iIndexEvt];

                    iIndexEvt++;
                    if (iIndexEvt < sMCArg->iNbEvent)
                    {
                        sPayArg =
                            (CMSCTS_PAYARG)(sMCArg->vPayoffParams[sMCArg->iIndexEvent[iIndexEvt]]);
                    }
                    else
                    {
                        sPayArg = NULL;
                    }
                }
            }

            if (iIncludePartialPV)
            {
                for (i = sKOAux->iStartFundIdx[iIndexStart];
                     i < sKOAux->iStartFundIdx[iIndexStart + 1];
                     i++)
                {
                    sFwdIVInfos->dMarketIV[0] -=
                        sDeal->sFunding->sAux->dMarketValue[i] * sDeal->dPayRec;
                    sFwdIVInfos->dModelIV[0] -=
                        sDeal->sFunding->sAux->dMarketValue[i] * sDeal->dPayRec;
                    sFwdIVInfos->dNewModelIV[0] -=
                        sDeal->sFunding->sAux->dMarketValue[i] * sDeal->dPayRec;
                }

                if (sDeal->sExotic->sAux->iFloatingType[sKOAux->iStartCpnIdx[iIndexStart]] < 2)
                {
                    sFwdIVInfos->dMarketIV[0] += sDeal->sExotic->sAux->dFixedPV * sDeal->dPayRec;
                    sFwdIVInfos->dModelIV[0] += sDeal->sExotic->sAux->dFixedPV * sDeal->dPayRec;
                    sFwdIVInfos->dNewModelIV[0] += sDeal->sExotic->sAux->dFixedPV * sDeal->dPayRec;
                }
            }
        }
    }

FREE_RETURN:

    return err;
}

Err cnscts_allocate_fwdiv(CMSCTS_FWDIV sFwdIVInfos, long lNbCall)
{
    sFwdIVInfos->lNbCall = lNbCall;

    sFwdIVInfos->lExeDates      = (long*)calloc(sFwdIVInfos->lNbCall, sizeof(long));
    sFwdIVInfos->dMarketIV      = (double*)calloc(sFwdIVInfos->lNbCall, sizeof(double));
    sFwdIVInfos->dModelIV       = (double*)calloc(sFwdIVInfos->lNbCall, sizeof(double));
    sFwdIVInfos->dNewModelIV    = (double*)calloc(sFwdIVInfos->lNbCall, sizeof(double));
    sFwdIVInfos->dFees          = (double*)calloc(sFwdIVInfos->lNbCall, sizeof(double));
    sFwdIVInfos->dOTCall        = (double*)calloc(sFwdIVInfos->lNbCall, sizeof(double));
    sFwdIVInfos->dOTCallPartial = (double*)calloc(sFwdIVInfos->lNbCall, sizeof(double));

    if (!sFwdIVInfos->lExeDates || !sFwdIVInfos->dMarketIV || !sFwdIVInfos->dModelIV ||
        !sFwdIVInfos->dNewModelIV || !sFwdIVInfos->dFees || !sFwdIVInfos->dOTCall ||
        !sFwdIVInfos->dOTCallPartial)
    {
        return "Memory allocation faillure in cmscts_fill_outputs";
    }

    return NULL;
}

void cmscts_free_fwdiv(CMSCTS_FWDIV sFwdIVInfos)
{
    if (sFwdIVInfos->lExeDates)
        free(sFwdIVInfos->lExeDates);
    if (sFwdIVInfos->dMarketIV)
        free(sFwdIVInfos->dMarketIV);
    if (sFwdIVInfos->dModelIV)
        free(sFwdIVInfos->dModelIV);
    if (sFwdIVInfos->dNewModelIV)
        free(sFwdIVInfos->dNewModelIV);
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
    sFwdIVInfos->dFees          = NULL;
    sFwdIVInfos->dOTCall        = NULL;
    sFwdIVInfos->dOTCallPartial = NULL;
}