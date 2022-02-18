/*******************************************************************************
 *
 * FUNCTION      : EuropeanGRFN.c
 *
 * PURPOSE       :
 *
 * DESCRIPTION   :
 *
 * CALLS         :
 *
 *******************************************************************************/
#include "EuropeanGRFN.h"

#include "num_h_interp.h"
#include "srtaccess.h"

/* ---------------------------------------------------------------------------

        Free the GRFNEuropeanCorrMat

  ---------------------------------------------------------------------------- */
char* free_GRFNEuropeanCorrMat(GRFNEuropeanCorrMat* ptr)
{
    if (ptr)
    {
        if (ptr->dCorrMat)
            free_dmatrix(ptr->dCorrMat, 0, ptr->lNbPair - 1, 0, ptr->lNbDate - 1);

        if (ptr->cName1)
            free_svector_size(ptr->cName1, 0, ptr->lNbPair - 1, ptr->lStringSize);

        if (ptr->cName2)
            free_svector_size(ptr->cName2, 0, ptr->lNbPair - 1, ptr->lStringSize);

        if (ptr->dDate)
            free_dvector(ptr->dDate, 0, ptr->lNbDate - 1);
        free(ptr);
    }

    return NULL;
}

/* ---------------------------------------------------------------------------

        SrtGetGRFNEuropeanCorr

  ---------------------------------------------------------------------------- */
Err SrtGetGRFNEuropeanCorr(/* Inputs */
                           GRFNEuropeanCorrMat* ptrCorrMat,
                           double               dTime,
                           char*                cName1,
                           char*                cName2,

                           /* OutPuts */
                           double* dCorr)
{
    Err err    = NULL;
    int iIndex = 0;

    /* Find the pair */
    iIndex = 0;
    while ((iIndex < ptrCorrMat->lNbPair) &&
           ((strcmp(cName1, ptrCorrMat->cName1[iIndex]) != 0) ||
            (strcmp(cName2, ptrCorrMat->cName2[iIndex]) != 0)) &&
           ((strcmp(cName2, ptrCorrMat->cName1[iIndex]) != 0) ||
            (strcmp(cName1, ptrCorrMat->cName2[iIndex]) != 0)))
        iIndex++;

    if (iIndex == ptrCorrMat->lNbPair)
        return "GRFNEuropeanCorrMatrix : underlying pair not found";

    /* Linear interp + flat extrapolation */
    interp(/* Inputs */
           ptrCorrMat->dDate,
           ptrCorrMat->dCorrMat[iIndex],
           ptrCorrMat->lNbDate,
           dTime,
           0, /* Method */
           /* OutPuts */
           dCorr);

    return err;
}

/* ---------------------------------------------------------------------------

        Free the SABRUnd

  ---------------------------------------------------------------------------- */
char* SABR_und_free_struct(void* pUndDesc)
{
    SrtUndPtr   UndPtr;
    SrtEqDesc*  pUndEqDesc = NULL;
    SrtSABRUnd* ptrSABRUnd = NULL;

    /* Casting */
    UndPtr     = (SrtUndPtr)pUndDesc;
    pUndEqDesc = (SrtEqDesc*)UndPtr->spec_desc;

    ptrSABRUnd = (SrtSABRUnd*)pUndEqDesc->spec;

    if (ptrSABRUnd->dTime)
        free(ptrSABRUnd->dTime);
    if (ptrSABRUnd->dFwd)
        free(ptrSABRUnd->dFwd);
    if (ptrSABRUnd->dSigma)
        free(ptrSABRUnd->dSigma);
    if (ptrSABRUnd->dAlpha)
        free(ptrSABRUnd->dAlpha);
    if (ptrSABRUnd->dBeta)
        free(ptrSABRUnd->dBeta);
    if (ptrSABRUnd->dRho)
        free(ptrSABRUnd->dRho);

    if (pUndEqDesc->spec)
        free(pUndEqDesc->spec);
    if (UndPtr->spec_desc)
        free(UndPtr->spec_desc);

    /* */
    return NULL;
}

/* ---------------------------------------------------------------------------

        Definition of the SABRUnd

  ---------------------------------------------------------------------------- */
char* SrtInitSABRUnd(
    char*   cUndName,
    char*   cYieldCurveName,
    long    lNbDate,
    double* dTime,
    double* dFwd,
    double* dSigma,
    double* dAlpha,
    double* dBeta,
    double* dRho)
{
    Err           err      = NULL;
    SrtUndDesc*   pUndDesc = NULL;
    SrtUndListPtr und_list;
    SrtSABRUnd*   ptrSABRUnd = NULL;
    SrtEqDesc*    pUndEqDesc = NULL;
    SrtCurvePtr   crv;

    /* */
    if ((crv = lookup_curve(cYieldCurveName)) == NULL)
    {
        return "Fatal: (SrtInitSABRUnd) Cannot find yield curve";
    }

    /* Create the new underlying */
    und_list = get_underlying_list();
    pUndDesc = (SrtUndDesc*)calloc(1, sizeof(SrtUndDesc));

    strcpy(pUndDesc->underl_name, cUndName);
    strupper(pUndDesc->underl_name);
    strip_white_space(pUndDesc->underl_name);
    strcpy(pUndDesc->underl_lbl, "SABR_UND");
    pUndDesc->underl_ccy = get_curve_ccy(crv);
    ;
    pUndDesc->underl_type = EQUITY_UND;

    /* Create an Eq und */
    pUndEqDesc          = (SrtEqDesc*)calloc(1, sizeof(SrtEqDesc));
    pUndDesc->spec_desc = (void*)pUndEqDesc;
    strcpy(pUndEqDesc->disc_name, cYieldCurveName);

    /* Allocation of a SABRUnd */
    ptrSABRUnd          = (SrtSABRUnd*)calloc(1, sizeof(SrtSABRUnd));
    ptrSABRUnd->dTime   = dTime;
    ptrSABRUnd->dFwd    = dFwd;
    ptrSABRUnd->dSigma  = dSigma;
    ptrSABRUnd->dAlpha  = dAlpha;
    ptrSABRUnd->dBeta   = dBeta;
    ptrSABRUnd->dRho    = dRho;
    ptrSABRUnd->lNbDate = lNbDate;

    pUndEqDesc->spec = (void*)ptrSABRUnd;

    err = srt_f_lstins(
        und_list,
        pUndDesc->underl_name,
        0.0,
        OBJ_PTR_UND,
        pUndDesc,
        &SABR_und_free_struct,
        &(pUndDesc->underl_ticker));

    if (err)
        return err;

    return NULL;
}

/* ---------------------------------------------------------------------------

        SrtGetSABRUndParams

  ---------------------------------------------------------------------------- */
Err SrtGetSABRUndParams(/* Inputs */
                        SrtSABRUnd* ptrSABRUnd,
                        double      dTime,

                        /* OutPuts */
                        double* dFwd,
                        double* dSigma,
                        double* dAlpha,
                        double* dBeta,
                        double* dRho)
{
    Err err = NULL;

    /* Linear interp + flat extrapolation */
    interp(/* Inputs */
           ptrSABRUnd->dTime,
           ptrSABRUnd->dFwd,
           ptrSABRUnd->lNbDate,
           dTime,
           0, /* Method */
           /* OutPuts */
           dFwd);

    interp(/* Inputs */
           ptrSABRUnd->dTime,
           ptrSABRUnd->dSigma,
           ptrSABRUnd->lNbDate,
           dTime,
           0, /* Method */
           /* OutPuts */
           dSigma);

    interp(/* Inputs */
           ptrSABRUnd->dTime,
           ptrSABRUnd->dAlpha,
           ptrSABRUnd->lNbDate,
           dTime,
           0, /* Method */
           /* OutPuts */
           dAlpha);

    interp(/* Inputs */
           ptrSABRUnd->dTime,
           ptrSABRUnd->dBeta,
           ptrSABRUnd->lNbDate,
           dTime,
           0, /* Method */
           /* OutPuts */
           dBeta);

    interp(/* Inputs */
           ptrSABRUnd->dTime,
           ptrSABRUnd->dRho,
           ptrSABRUnd->lNbDate,
           dTime,
           0, /* Method */
           /* OutPuts */
           dRho);

    return err;
}

/* ----------------------------------------------------------------------------------------------------------

        PayOff_GRFNEuropean

  -----------------------------------------------------------------------------------------------------------
*/

Err PayOff_GRFNEuropean(/* Underlying value */
                        int     iNbUnderlying,
                        double* UnderlyingValue,

                        /* For Payoff description */
                        int   iNbProduct,
                        void* PayOffParam,

                        /* OutPut */
                        double* dPV)
{
    Err                  err               = NULL;
    GRFNEuropean_Payoff* ProductCallsParam = NULL;
    int                  iNum;
    double               dValue;

    /* Payoff parameter Cast to ProductofCalls */
    ProductCallsParam = (GRFNEuropean_Payoff*)PayOffParam;

    /* Put the value of the underlyings for GFN call */
    for (iNum = 0; iNum < iNbUnderlying; iNum++)
        ProductCallsParam->local->smp.und[ProductCallsParam->iIndex[iNum]].sv[SPOT] =
            UnderlyingValue[iNum];

    /* Call the GRFN evaluation */
    err = FIRSTEvalEvent(
        ProductCallsParam->global, /* GRFNCOMMSTRUCT */
        ProductCallsParam->local,  /* FIRSTMktAtT */
        ProductCallsParam->global->num_cols,
        2,        /* type_eval:Default is 2*/
        NULL,     /*	fwd Default is NULL	*/
        NULL,     /* cur:Default is NULL	*/
        dPV,      /* all the PVs  */
        &dValue); /* the last PV */

    /* Return result */
    return err;
}

/* ----------------------------------------------------------------------------------------------------------

        SrtGRFNEuropean
        Pricing of a european GRFN

  -----------------------------------------------------------------------------------------------------------
*/

Err SrtGRFNEuropean(/* Grfn inputs */
                    int      numeventdates,
                    long*    eventdates,
                    long     tableauRows,
                    long     tableauCols,
                    char***  tableauStrings,
                    int**    tableauMask,
                    long     auxWidth,
                    long*    auxLen,
                    double** aux,

                    /* Copula Pricing Parameters */
                    CopulaType   iCopulaType,
                    long         lNumPaths,
                    long         lNbPoints,
                    SrtMCSamType iMCType,
                    ModelType    iModelType,
                    double       dNStdBMMCalib,

                    /* Dom Market */
                    char* cDomYcName,

                    /* Correlation Matrix */
                    GRFNEuropeanCorrMat* ptrCorrMat,

                    /* OutPuts */
                    double*  dPv,
                    double** grfn_ss /* For GRFN_celll */)
{
    char     dummy_und[9];
    double** dummy_vol_curve = NULL;
    double** dummy_tau_curve = NULL;

    Err          err      = NULL;
    int          free_str = 0;
    FIRSTAllMkts xStr;
    SrtGrfnParam defParm;
    int          forback;
    int          flag = 0;
    FIRSTMktAtT* evt  = NULL;

    int num_col = 0, max_num_df = 0, num_evt = 0, num_und = 0;

    SrtUndPtr* und_ptr = NULL;
    long*      evt_dts = NULL;
    double*    evt_tms = NULL;
    long       lToday;
    int*       am         = NULL;
    int**      num_df_mat = NULL;
    long***    df_mat_dts = NULL;
    double***  df_mat_tms = NULL;

    SrtUndPtr und = NULL;

    int iNbUnderlying, iNum, iNumCol, iNumRow, iNumEvt;

    double *dForward = NULL, *dVol = NULL, *dAlpha = NULL, *dBeta = NULL, *dRho = NULL,
           **dCorrelation = NULL;

    int* iIndex = NULL;

    double dMaturity;
    long   lEvtDte;

    Generic_Model* GenericModel  = NULL;
    void*          SpecificModel = NULL;

    Generic_Copula      CopulaParam;
    GRFNEuropean_Payoff GRFNEuropeanPayoff;

    char**      tabUndName = NULL;
    SrtCurvePtr Crv;
    char*       cCcy        = NULL;
    double      dPv_of_past = 0.0;

    GrfnEvent* GrfEvt = NULL;
    int        iNumUnd;

    /* -------------------------------------------------------------------------------------------------------
            GRFN PART
            -------------------------------------------------------------------------------------------------------
     */

    /*	Initialise the GRFN tableau */

    /*	First, initialise the param struct */
    err              = srt_f_set_default_GrfnParams(&defParm);
    defParm.force_mc = 0;

    /* Create a dummy domestic underlying to satisfy FIRSTInitMktStruct */
    strcpy(dummy_und, "DummyUnd");
    dummy_vol_curve       = (double**)dmatrix(0, 1, 0, 1);
    dummy_tau_curve       = (double**)dmatrix(0, 1, 0, 1);
    dummy_vol_curve[0][0] = 1;
    dummy_vol_curve[1][0] = 0.01;
    dummy_tau_curve[0][0] = 1;
    dummy_tau_curve[1][0] = 0.01;
    err                   = SrtInitIRUnd(
        dummy_und,
        cDomYcName,
        "LGM1F",
        1,
        1,
        dummy_vol_curve,
        1,
        1,
        dummy_tau_curve,
        0,
        0,
        0,
        0,
        0,
        0,
        0.0,
        0,
        0,
        NULL);
    if (err)
    {
        goto FREE_RETURN;
    }

    /* FIRSTInitMktStruct */
    err = FIRSTInitMktStruct(
        numeventdates,
        eventdates,
        tableauRows,
        tableauCols,
        tableauStrings,
        tableauMask,
        auxWidth,
        auxLen,
        aux,
        dummy_und,
        &defParm,
        &forback,
        &xStr);
    if (err)
    {
        goto FREE_RETURN;
    }

    free_str = 1;

    /*	Now, lookup underlyings involved */
    err = FIRSTGetUndFromDeal(&xStr, &num_und, &und_ptr);

    if (err)
    {
        goto FREE_RETURN;
    }

    /* look for the today date */
    if ((Crv = lookup_curve(cDomYcName)) == NULL)
    {
        goto FREE_RETURN;
    }

    /* Extract today and currency from this curve */
    lToday = get_today_from_curve(Crv);
    cCcy   = get_curve_ccy(Crv);

    /* Get number of columns */
    err = FIRSTGetNumColFromDeal(&xStr, &num_col);
    if (err)
    {
        goto FREE_RETURN;
    }

    /*	Get the maximum number of dfs required	*/
    err = FIRSTGetMaxNumDfFromDeal(&xStr, &max_num_df);
    if (err)
    {
        goto FREE_RETURN;
    }

    /*	Next, get the time steps */
    err = FIRSTGetEvtDatesFromDeal(&xStr, &num_evt, &evt_dts, &evt_tms);

    if (err)
    {
        goto FREE_RETURN;
    }

    /*	Get events and maturities of discount factors required	*/
    err = FIRSTGetEventInfoFromDeal(
        &xStr, num_evt, num_und, max_num_df, &evt, &am, &num_df_mat, &df_mat_dts, &df_mat_tms);

    if (err)
    {
        goto FREE_RETURN;
    }

    /* Fill the df : No diffusion on it */
    for (iNumEvt = 0; iNumEvt < num_evt; iNumEvt++)
    {
        if (xStr.evt[iNumEvt].evt)
        {
            GrfEvt = xStr.evt[iNumEvt].evt;
            for (iNumUnd = 0; iNumUnd < xStr.und_info->no_of_underlyings; iNumUnd++)
            {
                for (iNum = 0; iNum < GrfEvt->dflen[iNumUnd]; iNum++)
                {
                    GrfEvt->df[iNumUnd][iNum] =
                        swp_f_df(GrfEvt->d, GrfEvt->dfd[iNumUnd][iNum], cDomYcName);
                }
            }
        }
    }

    /* -------------------------------------------------------------------------------------------------------
            Init for Copula pricing
            -------------------------------------------------------------------------------------------------------
     */

    /* Calculate the number of SABR_UND , the number of underlying */
    iNbUnderlying = 0;
    for (iNum = 0; iNum < num_und; iNum++)
    {
        if (!strcmp(und_ptr[iNum]->underl_lbl, "SABR_UND"))
        {
            /* Check that ccy is the same */
            if (strcmp(und_ptr[iNum]->underl_ccy, cCcy))
            {
                err = " Ccy of all underlying should be the same as the DomCcy";
                goto FREE_RETURN;
            }
            iNbUnderlying++;
        }
    }

    /* Allocation */
    if (iNbUnderlying > 0)
    {
        tabUndName = svector_size(0, iNbUnderlying - 1, 255);
        iIndex     = (int*)calloc(iNbUnderlying, sizeof(int));
    }

    /* Fill TabUndName and iIndex */
    iNbUnderlying = 0;
    for (iNum = 0; iNum < num_und; iNum++)
    {
        if (!strcmp(und_ptr[iNum]->underl_lbl, "SABR_UND"))
        {
            iIndex[iNbUnderlying] = iNum;
            strcpy(tabUndName[iNbUnderlying], und_ptr[iNum]->underl_name);
            iNbUnderlying++;
        }
    }

    /* Initialisation of the payoff param */
    GRFNEuropeanPayoff.global = &xStr;
    GRFNEuropeanPayoff.local  = evt;
    GRFNEuropeanPayoff.iIndex = iIndex;

    /* Initialisation of the Generic Model */
    err = initGenericModel(/* Input */
                           iNbUnderlying,
                           iModelType,

                           /* For BMM type */
                           dNStdBMMCalib,

                           /* OutPut */
                           &GenericModel,
                           &SpecificModel);

    if (err)
    {
        goto FREE_RETURN;
    }

    /* Memory allocation */
    if (iNbUnderlying > 0)
    {
        dForward = (double*)calloc(iNbUnderlying, sizeof(double));
        dVol     = (double*)calloc(iNbUnderlying, sizeof(double));
        dAlpha   = (double*)calloc(iNbUnderlying, sizeof(double));
        dBeta    = (double*)calloc(iNbUnderlying, sizeof(double));
        dRho     = (double*)calloc(iNbUnderlying, sizeof(double));

        /* Gestion of the correlation */
        dCorrelation = (double**)dmatrix(0, iNbUnderlying - 1, 0, iNbUnderlying - 1);
    }

    /* Initialise the generic copula parameter */
    CopulaParam.iCopulaType   = iCopulaType;
    CopulaParam.iMCType       = iMCType;
    CopulaParam.iNbUnderlying = iNbUnderlying;
    CopulaParam.lNbPoints     = lNbPoints;
    CopulaParam.lNumPaths     = lNumPaths;
    CopulaParam.iIsCorrMatrix = 1;
    CopulaParam.dCorrMatrix   = dCorrelation;
    CopulaParam.OtherParams   = NULL;

    /* -------------------------------------------------------------------------------------------------------
             Copula pricing of the GRFN
            -------------------------------------------------------------------------------------------------------
     */

    for (iNumEvt = 0; iNumEvt < num_evt; iNumEvt++)
    {
        lEvtDte   = evt_dts[iNumEvt];
        dMaturity = evt_tms[iNumEvt];

        GRFNEuropeanPayoff.local = xStr.evt + iNumEvt;

        /* Fill Arrays */
        for (iNum = 0; iNum < iNbUnderlying; iNum++)
        {
            /* Get the und */
            und = lookup_und(tabUndName[iNum]);

            /* Fill SABR parameter */
            err = SrtGetSABRUndParams(/* Inputs */
                                      ((SrtEqDesc*)und->spec_desc)->spec,
                                      (double)lEvtDte,

                                      /* OutPuts */
                                      &dForward[iNum],
                                      &dVol[iNum],
                                      &dAlpha[iNum],
                                      &dBeta[iNum],
                                      &dRho[iNum]);
        }

        /* Fill Corr Matrix */
        for (iNumRow = 0; iNumRow < iNbUnderlying; iNumRow++)
            for (iNumCol = 0; iNumCol < iNbUnderlying; iNumCol++)
            {
                if (iNumRow == iNumCol)
                    dCorrelation[iNumRow][iNumRow] = 1.0;
                else
                {
                    err = SrtGetGRFNEuropeanCorr(/* Inputs */
                                                 ptrCorrMat,
                                                 (double)(lEvtDte),
                                                 tabUndName[iNumRow],
                                                 tabUndName[iNumCol],

                                                 /* OutPuts */
                                                 &(dCorrelation[iNumRow][iNumCol]));

                    if (err)
                        goto FREE_RETURN;
                }
            }

        /* Copula Pricing */
        err = Copula_SABRCalib_GenericPayOff_Price(
            /* Copula Parameter */
            &CopulaParam,

            /* PayOff parameters */
            tableauCols,
            &GRFNEuropeanPayoff,

            /* For generating Marginales Distributions */
            dMaturity,
            dForward, /* array */
            GenericModel,

            /* Sabr Parameters */
            dVol,
            dAlpha,
            dBeta,
            dRho,

            SRT_LOGNORMAL,

            /* Copula function */
            GaussianCopulaGetSamples,

            /* Payoff Function */
            PayOff_GRFNEuropean,

            /* Results */
            dPv);

        if (err)
            goto FREE_RETURN;

        /* Copy the fwd for GRFN_CELL */
        for (iNum = 0; iNum < tableauCols; iNum++)
            grfn_ss[iNumEvt][iNum] = dPv[iNum];
    }

    /* Reset to 0 */
    for (iNum = 0; iNum < tableauCols; iNum++)
        dPv[iNum] = 0.0;

    /* Cumul the Pv of each evt dates and discount */
    for (iNumEvt = 0; iNumEvt < num_evt; iNumEvt++)
    {
        lEvtDte = evt_dts[iNumEvt];
        for (iNum = 0; iNum < tableauCols; iNum++)
            dPv[iNum] += grfn_ss[iNumEvt][iNum] * swp_f_df(lToday, lEvtDte, cDomYcName);
    }

FREE_RETURN:
    if (dummy_vol_curve)
        free_dmatrix(dummy_vol_curve, 0, 1, 0, 1);
    if (dummy_tau_curve)
        free_dmatrix(dummy_tau_curve, 0, 1, 0, 1);

    if (free_str)
    {
        FIRSTFreeUndFromDeal(num_und, &und_ptr);

        FIRSTFreeEvtDatesFromDeal(num_evt, &evt_dts, &evt_tms);

        FIRSTFreeMktStruct(&xStr);
    }

    if (GenericModel)
        free(GenericModel);
    if (SpecificModel)
        free(SpecificModel);

    if (dForward)
        free(dForward);
    if (dVol)
        free(dVol);
    if (dAlpha)
        free(dAlpha);
    if (dBeta)
        free(dBeta);
    if (dRho)
        free(dRho);
    if (dCorrelation)
        free_dmatrix(dCorrelation, 0, iNbUnderlying - 1, 0, iNbUnderlying - 1);
    if (tabUndName)
        free_svector_size(tabUndName, 0, iNbUnderlying - 1, 255);
    if (iIndex)
        free(iIndex);

    return err;
}
