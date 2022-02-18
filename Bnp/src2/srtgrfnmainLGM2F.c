/**********************************************************************
 *      Name: SrtGrfnMainLgm2F.c                                      *
 *  Function: Entry point to GRFN with raw data                       *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: L.C.				                                      *
 *      Date: 03/01/01                                                *
 *--------------------------------------------------------------------*
 *    Inputs: Raw data from anywhere (Excel or 2020)                  *
 *   Returns:                                                         *
 *   Globals: Expects mkt and request list structures to exist        *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 18/10/95 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/
#include "BGMEval.h"
#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "LGM2FMC.h"
#include "LGM2Fgrfn.h"
#include "LGM2Fpde.h"
#include "SrtAccess.h"
#include "math.h"
#include "srt_h_all.h"

//**********************************************************************
//*                                                                    *
//*             LGM2F PDE with constant Tau and Tau TS                 *
//*                                                                    *
//**********************************************************************

char* SrtGrfnLGM2Fpde(
    char*    underlying,
    int      numeventdates,
    long*    eventdates,
    long     tableauRows,
    long     tableauCols,
    char***  tableauStrings,
    int**    tableauMask,
    long     auxWidth,
    long*    auxLen,
    double** aux,
    int      is_end_of_day_fixing,
    int      is_end_of_day_payment,
    int      nstept,
    int      nstepx,
    int*     nb_prod,
    double** prod_val)
{
    int           free_str = 0;
    FIRSTAllMkts  xStr;
    SrtGrfnParam  defParm;
    GRFNPARMLGM2F grfn_prm;
    int           forback;
    int           flag = 0;
    long          nstp;

    double next_d;

    double *evt_tms = NULL, *time = NULL, *sig_time = NULL, *sigma = NULL, *ifr = NULL,
           *date = NULL;

    long* evt_dts = NULL;

    int*   is_event = NULL;
    void** void_prm = NULL;

    double *dff = NULL, *gam1 = NULL, *gam2 = NULL, *gam12 = NULL;

    long   nbSig;
    double alpha, gamma, rho, tau;

    long today, spot_date;

    int num_col = 0, max_num_df = 0, num_evt = 0, num_und = 0;

    char *domestic_name, *yc;

    SrtUndPtr *und_ptr = NULL, und = NULL;

    int i, j;

    Err err = NULL;

    /*	Initialise the GRFN tableau */

    /*	First, initialise the param struct */
    err              = srt_f_set_default_GrfnParams(&defParm);
    defParm.force_mc = 0;

    /* End of Day Fixing */
    if (is_end_of_day_fixing)
    {
        defParm.end_of_day_fixings = SRT_YES;
        defParm.end_of_day_flg     = SRT_YES;
    }
    else
    {
        defParm.end_of_day_fixings = SRT_NO;
        defParm.end_of_day_flg     = SRT_NO;
    }

    /* End of Day Payment */
    if (is_end_of_day_payment)
    {
        defParm.end_of_day_payment = SRT_YES;
    }
    else
    {
        defParm.end_of_day_payment = SRT_NO;
    }

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
        underlying,
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

    if (num_und != 1)
    {
        err = "Product should involve only one underlying";
        goto FREE_RETURN;
    }

    und = und_ptr[0];

    /* look for the underlying name */
    und = lookup_und(underlying);
    if (!und)
    {
        err = "cannot find the underlying";
        goto FREE_RETURN;
    }

    if (get_mdltype_from_irund(und) == LGM)
    {
        domestic_name = und->underl_name;
    }
    else
    {
        err = "Model must be LGM2F";
        goto FREE_RETURN;
    }

    if (strcmp(domestic_name, und->underl_name))
    {
        err = "Tableau uses different underlying";
        goto FREE_RETURN;
    }

    /* look for the today date */
    today     = get_today_from_underlying(und);
    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    yc = (char*)get_ycname_from_irund(und);

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

    /* Get the term structure */
    err = Get_LGM2F_TermStructure(
        domestic_name, &sig_time, &sigma, &nbSig, &tau, &alpha, &gamma, &rho);
    if (err)
    {
        goto FREE_RETURN;
    }

    /* discretise in time			*/

    nstp = num_evt;

    time = (double*)calloc(nstp, sizeof(double));
    if (!time)
    {
        err = "Memory allocation error (1) in SrtGrfnLgm2FPde";
        goto FREE_RETURN;
    }
    memcpy(time, xStr.tms, nstp * sizeof(double));

    /*	Fill the time vector */
    err = fill_time_vector(&time, &nstp, 0, NULL, 0, NULL, nstept);
    if (err)
    {
        goto FREE_RETURN;
    }

    void_prm = (void**)calloc(nstp, sizeof(void*));
    is_event = (int*)calloc(nstp, sizeof(int));
    date     = (double*)calloc(nstp, sizeof(double));
    ifr      = (double*)calloc(nstp, sizeof(double));

    if (max_num_df > 0)
    {
        dff   = dvector(0, max_num_df - 1);
        gam1  = dvector(0, max_num_df - 1);
        gam2  = dvector(0, max_num_df - 1);
        gam12 = dvector(0, max_num_df - 1);
    }

    if (!void_prm || !is_event || !ifr || !date ||
        ((!dff || !gam1 || !gam2 || !gam12) && (max_num_df > 0)))
    {
        err = "Memory allocation failure";
        goto FREE_RETURN;
    }

    j = xStr.num_evt - 1;

    next_d = evt_dts[j] + 1;

    for (i = nstp - 1; i >= 0; i--)
    {
        date[i] = ((long)(today + time[i] * 365.0 + 1.0E-08)) * 1.0;
        // time[i] = (date[i] - today) / 365.0;

        if (next_d > date[i])
        {
            ifr[i] = swp_f_zr(date[i], next_d, yc);
        }
        else
        {
            ifr[i] = swp_f_zr(date[i], date[i] + 1, yc);
        }

        if (j >= 0 && fabs(time[i] - evt_tms[j]) < 1.0E-08)
        {
            grfn_prm = malloc(sizeof(grfn_parm_lgm));

            grfn_prm->global = &xStr;
            grfn_prm->local  = xStr.evt + j;

            grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
            grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
            grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];

            grfn_prm->dff   = dff;
            grfn_prm->gam1  = gam1;
            grfn_prm->gam2  = gam2;
            grfn_prm->gam12 = gam12;

            is_event[i] = 1;
            void_prm[i] = (void*)grfn_prm;

            j--;
            while (j >= 0 && xStr.evt[j].evt == NULL)
            {
                j--;
            }
        }
        else if (j >= 0 && xStr.am[j])
        {
            grfn_prm         = malloc(sizeof(grfn_parm_lgm));
            grfn_prm->global = &xStr;
            grfn_prm->local  = xStr.evt + j;

            grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
            grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
            grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];

            grfn_prm->dff   = dff;
            grfn_prm->gam1  = gam1;
            grfn_prm->gam2  = gam2;
            grfn_prm->gam12 = gam12;

            is_event[i] = 1;
            void_prm[i] = (void*)grfn_prm;
        }
        else
        {
            is_event[i] = 0;
            void_prm[i] = NULL;
        }
        next_d = date[i];
    }

    /* put the middle point */
    /*
    for (i=0; i<nstp-1; i++)
    {
            ifr[i] = 0.5 * (ifr[i] + ifr[i+1]);
    }
    */

    /*	Eventually! call to function */

    *prod_val = dvector(0, num_col - 1);

    if (!*prod_val)
    {
        err = "Memory allocation failure";
        goto FREE_RETURN;
    }

    /* Check to see if there are any event dates */
    if (eventdates[numeventdates - 1] >= today)
    {
        err = lgm2f_adi(
            nstp,
            time,
            date,
            nstepx,
            1.0 / tau,
            sig_time,
            sigma,
            nbSig,
            alpha,
            gamma,
            rho,
            void_prm,
            is_event,
            ifr,
            yc,
            payoff_lgm2fTau_pde,
            num_col,
            *prod_val);
    }
    if (err)
    {
        goto FREE_RETURN;
    }

    *nb_prod = num_col;

    /*	Add PV of Past */
    (*prod_val)[num_col - 1] += xStr.gd->pv_of_past;

FREE_RETURN:

    if (free_str)
    {
        FIRSTFreeUndFromDeal(num_und, &und_ptr);

        FIRSTFreeEvtDatesFromDeal(nstp, &evt_dts, &evt_tms);

        FIRSTFreeMktStruct(&xStr);
    }

    if (void_prm)
    {
        for (i = 0; i < nstp; i++)
        {
            if (void_prm[i])
            {
                grfn_prm = (GRFNPARMLGM2F)void_prm[i];
                free(grfn_prm);
            }
        }

        free(void_prm);
    }

    if (max_num_df > 0)
    {
        if (dff)
            free_dvector(dff, 0, max_num_df - 1);
        if (gam1)
            free_dvector(gam1, 0, max_num_df - 1);
        if (gam2)
            free_dvector(gam2, 0, max_num_df - 1);
        if (gam12)
            free_dvector(gam12, 0, max_num_df - 1);
    }

    if (is_event)
        free(is_event);
    if (date)
        free(date);
    if (ifr)
        free(ifr);
    if (time)
        free(time);

    if (sig_time)
        free(sig_time);
    if (sigma)
        free(sigma);

    return err;
}

char* SrtGrfnLGM2FTaupde(
    char*    underlying,
    int      numeventdates,
    long*    eventdates,
    long     tableauRows,
    long     tableauCols,
    char***  tableauStrings,
    int**    tableauMask,
    long     auxWidth,
    long*    auxLen,
    double** aux,
    int      is_end_of_day_fixing,
    int      is_end_of_day_payment,
    int      nstept,
    int      nstepx,
    int*     nb_prod,
    double** prod_val)
{
    int           free_str = 0;
    FIRSTAllMkts  xStr;
    SrtGrfnParam  defParm;
    GRFNPARMLGM2F grfn_prm;
    int           forback;
    int           flag = 0;
    long          nstp;

    double next_d;

    double *evt_tms = NULL, *time = NULL, *sigma_time = NULL, *lambda_time = NULL, *sigma = NULL,
           *lambda = NULL, *ifr = NULL, *date = NULL;

    long* evt_dts = NULL;

    int*   is_event = NULL;
    void** void_prm = NULL;

    double *dff = NULL, *gam1 = NULL, *gam2 = NULL, *gam12 = NULL;

    long   nb_sigma, nb_lambda;
    double alpha, gamma, rho;

    long today, spot_date;

    int num_col = 0, max_num_df = 0, num_evt = 0, num_und = 0;

    char *domestic_name, *yc;

    SrtUndPtr *und_ptr = NULL, und = NULL;

    int i, j;

    Err err = NULL;

    /*	Initialise the GRFN tableau */

    /*	First, initialise the param struct */
    err              = srt_f_set_default_GrfnParams(&defParm);
    defParm.force_mc = 0;

    /* End of Day Fixing */
    if (is_end_of_day_fixing)
    {
        defParm.end_of_day_fixings = SRT_YES;
        defParm.end_of_day_flg     = SRT_YES;
    }
    else
    {
        defParm.end_of_day_fixings = SRT_NO;
        defParm.end_of_day_flg     = SRT_NO;
    }

    /* End of Day Payment */
    if (is_end_of_day_payment)
    {
        defParm.end_of_day_payment = SRT_YES;
    }
    else
    {
        defParm.end_of_day_payment = SRT_NO;
    }

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
        underlying,
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

    if (num_und != 1)
    {
        err = "Product should involve only one underlying";
        goto FREE_RETURN;
    }

    und = und_ptr[0];

    /* look for the underlying name */
    und = lookup_und(underlying);
    if (!und)
    {
        err = "cannot find the underlying";
        goto FREE_RETURN;
    }

    if (get_mdltype_from_irund(und) == LGM)
    {
        domestic_name = und->underl_name;
    }
    else
    {
        err = "Model must be LGM2F";
        goto FREE_RETURN;
    }

    if (strcmp(domestic_name, und->underl_name))
    {
        err = "Tableau uses different underlying";
        goto FREE_RETURN;
    }

    /* look for the today date */
    today     = get_today_from_underlying(und);
    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    yc = (char*)get_ycname_from_irund(und);

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

    /* Get the term structure */
    err = Get_LGM2F_TermStructure2(
        domestic_name,
        &sigma,
        &sigma_time,
        &nb_sigma,
        &lambda,
        &lambda_time,
        &nb_lambda,
        &alpha,
        &gamma,
        &rho);
    if (err)
    {
        goto FREE_RETURN;
    }

    /* discretise in time			*/

    nstp = num_evt;

    time = (double*)calloc(nstp, sizeof(double));
    if (!time)
    {
        err = "Memory allocation error (1) in SrtGrfnLgm2FPde";
        goto FREE_RETURN;
    }
    memcpy(time, xStr.tms, nstp * sizeof(double));

    /*	Fill the time vector */
    err = fill_time_vector(&time, &nstp, 0, NULL, 0, NULL, nstept);
    if (err)
    {
        goto FREE_RETURN;
    }

    void_prm = (void**)calloc(nstp, sizeof(void*));
    is_event = (int*)calloc(nstp, sizeof(int));
    date     = (double*)calloc(nstp, sizeof(double));
    ifr      = (double*)calloc(nstp, sizeof(double));

    if (max_num_df > 0)
    {
        dff   = dvector(0, max_num_df - 1);
        gam1  = dvector(0, max_num_df - 1);
        gam2  = dvector(0, max_num_df - 1);
        gam12 = dvector(0, max_num_df - 1);
    }

    if (!void_prm || !is_event || !ifr || !date ||
        ((!dff || !gam1 || !gam2 || !gam12) && (max_num_df > 0)))
    {
        err = "Memory allocation failure";
        goto FREE_RETURN;
    }

    j = xStr.num_evt - 1;

    next_d = evt_dts[j] + 1;

    for (i = nstp - 1; i >= 0; i--)
    {
        date[i] = ((long)(today + time[i] * 365.0 + 1.0E-08)) * 1.0;
        // time[i] = (date[i] - today) / 365.0;

        if (next_d > date[i])
        {
            ifr[i] = swp_f_zr(date[i], next_d, yc);
        }
        else
        {
            ifr[i] = swp_f_zr(date[i], date[i] + 1, yc);
        }

        if (j >= 0 && fabs(time[i] - evt_tms[j]) < 1.0E-08)
        {
            grfn_prm = malloc(sizeof(grfn_parm_lgm));

            grfn_prm->global = &xStr;
            grfn_prm->local  = xStr.evt + j;

            grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
            grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
            grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];

            grfn_prm->dff   = dff;
            grfn_prm->gam1  = gam1;
            grfn_prm->gam2  = gam2;
            grfn_prm->gam12 = gam12;

            is_event[i] = 1;
            void_prm[i] = (void*)grfn_prm;

            j--;
            while (j >= 0 && xStr.evt[j].evt == NULL)
            {
                j--;
            }
        }
        else if (j >= 0 && xStr.am[j])
        {
            grfn_prm         = malloc(sizeof(grfn_parm_lgm));
            grfn_prm->global = &xStr;
            grfn_prm->local  = xStr.evt + j;

            grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
            grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
            grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];

            grfn_prm->dff   = dff;
            grfn_prm->gam1  = gam1;
            grfn_prm->gam2  = gam2;
            grfn_prm->gam12 = gam12;

            is_event[i] = 1;
            void_prm[i] = (void*)grfn_prm;
        }
        else
        {
            is_event[i] = 0;
            void_prm[i] = NULL;
        }
        next_d = date[i];
    }

    /*	Eventually! call to function */

    *prod_val = dvector(0, num_col - 1);

    if (!*prod_val)
    {
        err = "Memory allocation failure";
        goto FREE_RETURN;
    }

    /* Slower version but much more accurate
            in case of negative tau	*/
    /*
    err = lgm2fTau2_adi(	nstp,
                                            time,
                                            date,
                                            nstepx,
                                            0,
                                            sigma,
                                            sigma_time,
                                            nb_sigma,
                                            lambda,
                                            lambda_time,
                                            nb_lambda,
                                            alpha,
                                            gamma,
                                            rho,
                                            void_prm,
                                            is_event,
                                            ifr,
                                            yc,
                                            payoff_lgm2fTau_pde,
                                            num_col,
                                            *prod_val);	*/

    err = lgm2fTau_adi(
        nstp,
        time,
        date,
        nstepx,
        0,
        sigma,
        sigma_time,
        nb_sigma,
        lambda,
        lambda_time,
        nb_lambda,
        alpha,
        gamma,
        rho,
        void_prm,
        is_event,
        ifr,
        yc,
        payoff_lgm2fTau_pde,
        num_col,
        *prod_val);

    if (err)
    {
        goto FREE_RETURN;
    }

    *nb_prod = num_col;

    /*	Add PV of Past */
    (*prod_val)[num_col - 1] += xStr.gd->pv_of_past;

FREE_RETURN:

    if (free_str)
    {
        FIRSTFreeUndFromDeal(num_und, &und_ptr);

        FIRSTFreeEvtDatesFromDeal(nstp, &evt_dts, &evt_tms);

        FIRSTFreeMktStruct(&xStr);
    }

    if (void_prm)
    {
        for (i = 0; i < nstp; i++)
        {
            if (void_prm[i])
            {
                grfn_prm = (GRFNPARMLGM2F)void_prm[i];
                free(grfn_prm);
            }
        }

        free(void_prm);
    }

    if (max_num_df > 0)
    {
        if (dff)
            free_dvector(dff, 0, max_num_df - 1);
        if (gam1)
            free_dvector(gam1, 0, max_num_df - 1);
        if (gam2)
            free_dvector(gam2, 0, max_num_df - 1);
        if (gam12)
            free_dvector(gam12, 0, max_num_df - 1);
    }

    if (is_event)
        free(is_event);
    if (date)
        free(date);
    if (ifr)
        free(ifr);
    if (time)
        free(time);

    if (sigma)
        free(sigma);
    if (sigma_time)
        free(sigma_time);
    if (lambda)
        free(lambda);
    if (lambda_time)
        free(lambda_time);

    return err;
}

char* SrtGrfnLGM2FTaupde2(
    char*    underlying,
    int      numeventdates,
    long*    eventdates,
    long     tableauRows,
    long     tableauCols,
    char***  tableauStrings,
    int**    tableauMask,
    long     auxWidth,
    long*    auxLen,
    double** aux,
    int      is_end_of_day_fixing,
    int      is_end_of_day_payment,
    int      nstept,
    int      nstepx,
    int*     nb_prod,
    double** prod_val)
{
    int           free_str = 0;
    FIRSTAllMkts  xStr;
    SrtGrfnParam  defParm;
    GRFNPARMLGM2F grfn_prm;
    int           forback;
    int           flag = 0;
    long          nstp;

    double next_d;

    double *evt_tms = NULL, *time = NULL, *sigma_time = NULL, *lambda_time = NULL, *sigma = NULL,
           *lambda = NULL, *ifr = NULL, *date = NULL;

    long* evt_dts = NULL;

    int*   is_event = NULL;
    void** void_prm = NULL;

    double *dff = NULL, *gam1 = NULL, *gam2 = NULL, *gam12 = NULL;

    long   nb_sigma, nb_lambda;
    double alpha, gamma, rho;

    long today, spot_date;

    int num_col = 0, max_num_df = 0, num_evt = 0, num_und = 0;

    char *domestic_name, *yc;

    SrtUndPtr *und_ptr = NULL, und = NULL;

    int i, j;

    Err err = NULL;

    /*	Initialise the GRFN tableau */

    /*	First, initialise the param struct */
    err              = srt_f_set_default_GrfnParams(&defParm);
    defParm.force_mc = 0;

    /* End of Day Fixing */
    if (is_end_of_day_fixing)
    {
        defParm.end_of_day_fixings = SRT_YES;
        defParm.end_of_day_flg     = SRT_YES;
    }
    else
    {
        defParm.end_of_day_fixings = SRT_NO;
        defParm.end_of_day_flg     = SRT_NO;
    }

    /* End of Day Payment */
    if (is_end_of_day_payment)
    {
        defParm.end_of_day_payment = SRT_YES;
    }
    else
    {
        defParm.end_of_day_payment = SRT_NO;
    }

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
        underlying,
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

    if (num_und != 1)
    {
        err = "Product should involve only one underlying";
        goto FREE_RETURN;
    }

    und = und_ptr[0];

    /* look for the underlying name */
    und = lookup_und(underlying);
    if (!und)
    {
        err = "cannot find the underlying";
        goto FREE_RETURN;
    }

    if (get_mdltype_from_irund(und) == LGM)
    {
        domestic_name = und->underl_name;
    }
    else
    {
        err = "Model must be LGM2F";
        goto FREE_RETURN;
    }

    if (strcmp(domestic_name, und->underl_name))
    {
        err = "Tableau uses different underlying";
        goto FREE_RETURN;
    }

    /* look for the today date */
    today     = get_today_from_underlying(und);
    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    yc = (char*)get_ycname_from_irund(und);

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

    /* Get the term structure */
    err = Get_LGM2F_TermStructure2(
        domestic_name,
        &sigma,
        &sigma_time,
        &nb_sigma,
        &lambda,
        &lambda_time,
        &nb_lambda,
        &alpha,
        &gamma,
        &rho);
    if (err)
    {
        goto FREE_RETURN;
    }

    /* discretise in time			*/

    nstp = num_evt;

    time = (double*)calloc(nstp, sizeof(double));
    if (!time)
    {
        err = "Memory allocation error (1) in SrtGrfnLgm2FPde";
        goto FREE_RETURN;
    }
    memcpy(time, xStr.tms, nstp * sizeof(double));

    /*	Fill the time vector */
    err = fill_time_vector(&time, &nstp, 0, NULL, 0, NULL, nstept);
    if (err)
    {
        goto FREE_RETURN;
    }

    void_prm = (void**)calloc(nstp, sizeof(void*));
    is_event = (int*)calloc(nstp, sizeof(int));
    date     = (double*)calloc(nstp, sizeof(double));
    ifr      = (double*)calloc(nstp, sizeof(double));

    if (max_num_df > 0)
    {
        dff   = dvector(0, max_num_df - 1);
        gam1  = dvector(0, max_num_df - 1);
        gam2  = dvector(0, max_num_df - 1);
        gam12 = dvector(0, max_num_df - 1);
    }

    if (!void_prm || !is_event || !ifr || !date ||
        ((!dff || !gam1 || !gam2 || !gam12) && (max_num_df > 0)))
    {
        err = "Memory allocation failure";
        goto FREE_RETURN;
    }

    j = xStr.num_evt - 1;

    next_d = evt_dts[j] + 1;

    for (i = nstp - 1; i >= 0; i--)
    {
        date[i] = ((long)(today + time[i] * 365.0 + 1.0E-08)) * 1.0;
        // time[i] = (date[i] - today) / 365.0;

        if (next_d > date[i])
        {
            ifr[i] = swp_f_zr(date[i], next_d, yc);
        }
        else
        {
            ifr[i] = swp_f_zr(date[i], date[i] + 1, yc);
        }

        if (j >= 0 && fabs(time[i] - evt_tms[j]) < 1.0E-08)
        {
            grfn_prm = malloc(sizeof(grfn_parm_lgm));

            grfn_prm->global = &xStr;
            grfn_prm->local  = xStr.evt + j;

            grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
            grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
            grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];

            grfn_prm->dff   = dff;
            grfn_prm->gam1  = gam1;
            grfn_prm->gam2  = gam2;
            grfn_prm->gam12 = gam12;

            is_event[i] = 1;
            void_prm[i] = (void*)grfn_prm;

            j--;
            while (j >= 0 && xStr.evt[j].evt == NULL)
            {
                j--;
            }
        }
        else if (j >= 0 && xStr.am[j])
        {
            grfn_prm         = malloc(sizeof(grfn_parm_lgm));
            grfn_prm->global = &xStr;
            grfn_prm->local  = xStr.evt + j;

            grfn_prm->num_df = xStr.evt[j].evt->dflen[0];
            grfn_prm->df_tms = xStr.evt[j].evt->dft[0];
            grfn_prm->df_dts = xStr.evt[j].evt->dfd[0];

            grfn_prm->dff   = dff;
            grfn_prm->gam1  = gam1;
            grfn_prm->gam2  = gam2;
            grfn_prm->gam12 = gam12;

            is_event[i] = 1;
            void_prm[i] = (void*)grfn_prm;
        }
        else
        {
            is_event[i] = 0;
            void_prm[i] = NULL;
        }
        next_d = date[i];
    }

    /*	Eventually! call to function */

    *prod_val = dvector(0, num_col - 1);

    if (!*prod_val)
    {
        err = "Memory allocation failure";
        goto FREE_RETURN;
    }

    err = lgm2fTau2_adi(
        nstp,
        time,
        date,
        nstepx,
        0,
        sigma,
        sigma_time,
        nb_sigma,
        lambda,
        lambda_time,
        nb_lambda,
        alpha,
        gamma,
        rho,
        void_prm,
        is_event,
        ifr,
        yc,
        payoff_lgm2fTau_pde,
        num_col,
        *prod_val);

    if (err)
    {
        goto FREE_RETURN;
    }

    *nb_prod = num_col;

    /*	Add PV of Past */
    (*prod_val)[num_col - 1] += xStr.gd->pv_of_past;

FREE_RETURN:

    if (free_str)
    {
        FIRSTFreeUndFromDeal(num_und, &und_ptr);

        FIRSTFreeEvtDatesFromDeal(nstp, &evt_dts, &evt_tms);

        FIRSTFreeMktStruct(&xStr);
    }

    if (void_prm)
    {
        for (i = 0; i < nstp; i++)
        {
            if (void_prm[i])
            {
                grfn_prm = (GRFNPARMLGM2F)void_prm[i];
                free(grfn_prm);
            }
        }

        free(void_prm);
    }

    if (max_num_df > 0)
    {
        if (dff)
            free_dvector(dff, 0, max_num_df - 1);
        if (gam1)
            free_dvector(gam1, 0, max_num_df - 1);
        if (gam2)
            free_dvector(gam2, 0, max_num_df - 1);
        if (gam12)
            free_dvector(gam12, 0, max_num_df - 1);
    }

    if (is_event)
        free(is_event);
    if (date)
        free(date);
    if (ifr)
        free(ifr);
    if (time)
        free(time);

    if (sigma)
        free(sigma);
    if (sigma_time)
        free(sigma_time);
    if (lambda)
        free(lambda);
    if (lambda_time)
        free(lambda_time);

    return err;
}

//**********************************************************************
//*                                                                    *
//*                   LGM2F MC with constant Tau                       *
//*                                                                    *
//**********************************************************************

static double Phi_Func(double x, double T, double s, double t)
{
    double result;

    result = (exp(-x * (T - t)) - exp(-x * (T - s))) / x;

    return result;
}

static double Etha_Func(double x, double T, double s, double t)
{
    double result;

    result = (t - s - Phi_Func(x, T, s, t)) / x;

    return result;
}

static double B_Func(double x, double T, double s, double t)
{
    double result;

    result = -(t * t - s * s) / 2.0 +
             1.0 / x * ((t - 1.0 / x) * exp(-x * (T - t)) - (s - 1.0 / x) * exp(-x * (T - s)));

    return result;
}

static double Psi_Func(double x, double y, double T, double s, double t)
{
    double result;

    result = 1.0 / (x * y) *
             (t - s - Phi_Func(x, T, s, t) - Phi_Func(y, T, s, t) + Phi_Func(x + y, T, s, t));

    return result;
}

static double Psi2_Func(double x, double y, double Tx, double Ty, double s, double t)
{
    double result;

    result = 1.0 / (x * y) *
             (t - s - Phi_Func(x, Tx, s, t) - Phi_Func(y, Ty, s, t) +
              exp(-x * (Tx - Ty)) * Phi_Func(x + y, Ty, s, t));

    return result;
}

Err fill_mc_init_lgm2f(
    int       do_jump,
    long      pay_date,
    double    pay_time,
    double*   date,
    double*   time,
    long      nb_dates,
    double*   sig_dates,
    long      nb_sig_dates,
    double*   sig_curve_dom,
    double    lda_dom,
    double    alpha_dom,
    double    gamma_dom,
    double    rho_dom,
    char*     dom_yc,
    double*   dom_fwd1,
    double*   dom_fwd2,
    double*   dom_exp1,
    double*   dom_exp2,
    double*   dom_phi1,
    double*   dom_phi2,
    double*   dom_phi12,
    double*   dom_gam1_fwd,
    double*   dom_gam2_fwd,
    double*   dom_bond_pay,
    double*   dom_gam1_pay,
    double*   dom_gam2_pay,
    double*** covar)
{
    double alpha_dom2;
    double lda_dom2, sig_dom, sig_dom2;
    double T1, T2, start_date, end_date, start_mat, end_mat;

    double phif_dom1, phif_dom2, phif_dom11, phif_dom12, phif_dom21, phif_dom22;
    double phi_dom1, phi_dom2, phi_dom12;
    double etha_dom1_Tend, etha_dom1_Tpay, etha_dom2_Tend, etha_dom2_Tpay;
    double expect_dom1, expect_dom2, expect_dom12, expect_dom21;
    double adj_pay_dom1, adj_pay_dom12, adj_pay_dom2, adj_pay_dom21;

    double mat_pay, mat, pay_mat, pay_mat2;
    double zc_dom, zc_pay;
    double exp_dom_mat, exp_dom_mat2, exp_dom_pay, exp_dom_pay2;
    int    i, k;
    long   StartIndex, EndIndex;

    double** cov;

    Err err = NULL;

    alpha_dom2 = alpha_dom * alpha_dom;
    lda_dom2   = lda_dom + gamma_dom;

    dom_fwd1[0] = 0.0;
    dom_fwd2[0] = 0.0;

    dom_phi1[0]  = 0.0;
    dom_phi2[0]  = 0.0;
    dom_phi12[0] = 0.0;

    mat_pay = (pay_date - date[0]) / 365.0;

    start_date = date[0];
    start_mat  = time[0];
    StartIndex = Get_Index(start_mat, sig_dates, nb_sig_dates);

    for (k = 0; k < nb_dates - 1; k++)
    {
        end_date = date[k + 1];

        if (do_jump)
        {
            pay_date = (long)(date[k + 1] + 0.5);
            pay_time = time[k + 1];
        }

        end_mat  = time[k + 1];
        EndIndex = Get_Index(end_mat, sig_dates, nb_sig_dates);
        mat      = (end_mat - start_mat);
        pay_mat  = (pay_date - start_date) / 365.0;
        pay_mat2 = (pay_date - end_date) / 365.0;

        zc_dom = swp_f_zr(start_date, end_date, dom_yc);

        // QTexpect of the log Fx

        dom_gam1_fwd[k] = 1.0 / lda_dom * (1 - exp(-lda_dom * mat));
        dom_gam2_fwd[k] = 1.0 / lda_dom2 * (1 - exp(-lda_dom2 * mat));

        // For Df to pay date
        zc_pay = swp_f_zr(start_date, pay_date, dom_yc);

        dom_gam1_pay[k] = 1.0 / lda_dom * (1 - exp(-lda_dom * pay_mat));
        dom_gam2_pay[k] = 1.0 / lda_dom2 * (1 - exp(-lda_dom2 * pay_mat));
        dom_bond_pay[k] =
            exp(-zc_pay * pay_mat -
                0.5 * (dom_gam1_pay[k] * dom_gam1_pay[k] * dom_phi1[k] +
                       dom_gam2_pay[k] * dom_gam2_pay[k] * dom_phi2[k]) -
                dom_gam1_pay[k] * dom_gam2_pay[k] * dom_phi12[k]);

        //	Variables initialisation
        exp_dom_mat  = exp(-lda_dom * mat);
        exp_dom_mat2 = exp(-lda_dom2 * mat);
        exp_dom_pay  = exp(-lda_dom * pay_mat2);
        exp_dom_pay2 = exp(-lda_dom2 * pay_mat2);

        phi_dom1 = phi_dom2 = phi_dom12 = 0.0;
        expect_dom1 = expect_dom12 = expect_dom2 = expect_dom21 = 0.0;
        adj_pay_dom1 = adj_pay_dom12 = adj_pay_dom2 = adj_pay_dom21 = 0.0;

        for (i = StartIndex; i < EndIndex + 1; i++)
        {
            if (i > StartIndex)
            {
                T1 = sig_dates[i - 1];
            }
            else
            {
                // First part
                T1 = start_mat;
            }

            if (i == EndIndex || StartIndex == EndIndex)
            {
                // Last part
                T2 = end_mat;
            }
            else
            {
                T2 = sig_dates[i];
            }

            sig_dom = sig_dom2 = sig_curve_dom[i];
            sig_dom2 *= sig_dom2;

            phif_dom1  = Phi_Func(lda_dom, end_mat, T1, T2);
            phif_dom2  = Phi_Func(lda_dom2, end_mat, T1, T2);
            phif_dom11 = Phi_Func(2.0 * lda_dom, end_mat, T1, T2);
            phif_dom22 = Phi_Func(2.0 * lda_dom2, end_mat, T1, T2);
            phif_dom12 = Phi_Func(lda_dom + lda_dom2, end_mat, T1, T2);
            phif_dom21 = phif_dom12;

            etha_dom1_Tend = Etha_Func(lda_dom, end_mat, T1, T2);
            etha_dom1_Tpay = Etha_Func(lda_dom, mat_pay, T1, T2);
            etha_dom2_Tend = Etha_Func(lda_dom2, end_mat, T1, T2);
            etha_dom2_Tpay = Etha_Func(lda_dom2, mat_pay, T1, T2);

            // Domestic phi and expectations
            phi_dom1 += sig_dom2 * phif_dom11;
            phi_dom2 += sig_dom2 * phif_dom22;
            phi_dom12 += sig_dom2 * phif_dom12;

            expect_dom1 += sig_dom2 * (phif_dom1 - phif_dom11);
            expect_dom12 += sig_dom2 * (phif_dom1 - phif_dom12);

            adj_pay_dom1 += sig_dom2 * (phif_dom1 - exp_dom_pay * phif_dom11);
            adj_pay_dom12 += sig_dom2 * (phif_dom1 - exp_dom_pay2 * phif_dom12);

            expect_dom2 += sig_dom2 * (phif_dom2 - phif_dom22);
            expect_dom21 += sig_dom2 * (phif_dom2 - phif_dom21);

            adj_pay_dom2 += sig_dom2 * (phif_dom2 - exp_dom_pay2 * phif_dom22);
            adj_pay_dom21 += sig_dom2 * (phif_dom2 - exp_dom_pay * phif_dom21);
        }

        dom_exp1[k + 1] = exp_dom_mat;
        dom_exp2[k + 1] = exp_dom_mat2;

        dom_phi1[k + 1] = dom_phi1[k] * exp_dom_mat * exp_dom_mat + phi_dom1;
        dom_phi2[k + 1] = dom_phi2[k] * exp_dom_mat2 * exp_dom_mat2 + alpha_dom2 * phi_dom2;
        dom_phi12[k + 1] =
            dom_phi12[k] * exp_dom_mat * exp_dom_mat2 + rho_dom * alpha_dom * phi_dom12;

        dom_fwd1[k + 1] =
            dom_phi1[k] * exp_dom_mat * Phi_Func(-lda_dom, start_mat, start_mat, end_mat) +
            expect_dom1 / lda_dom +
            dom_phi12[k] * exp_dom_mat * Phi_Func(-lda_dom2, start_mat, start_mat, end_mat) +
            rho_dom * alpha_dom * expect_dom12 / lda_dom2 - adj_pay_dom1 / lda_dom -
            rho_dom * alpha_dom * adj_pay_dom12 / lda_dom2;

        dom_fwd2[k + 1] =
            dom_phi2[k] * exp_dom_mat2 * Phi_Func(-lda_dom2, start_mat, start_mat, end_mat) +
            alpha_dom2 * expect_dom2 / lda_dom2 +
            dom_phi12[k] * exp_dom_mat2 * Phi_Func(-lda_dom, start_mat, start_mat, end_mat) +
            rho_dom * alpha_dom * expect_dom21 / lda_dom - adj_pay_dom2 * alpha_dom2 / lda_dom2 -
            rho_dom * alpha_dom * adj_pay_dom21 / lda_dom;

        cov = covar[k + 1];

        cov[0][0] = phi_dom1;
        cov[0][1] = cov[1][0] = rho_dom * alpha_dom * phi_dom12;
        cov[1][1]             = alpha_dom2 * phi_dom2;
        cov[1][0]             = cov[0][1];

        start_date = end_date;
        start_mat  = end_mat;
        StartIndex = EndIndex;
    }

    dom_gam1_pay[nb_dates - 1] = 0.0;
    dom_gam2_pay[nb_dates - 1] = 0.0;
    dom_bond_pay[nb_dates - 1] = 1.0;

    return err;
}

// LGM2F MC with constant lambda
char* SrtGrfnLGM2FMC(
    char*     lgmund,
    int       numeventdates,
    long*     eventdates,
    long      tableauRows,
    long*     tableauCols,
    char***   tableauStrings,
    int**     tableauMask,
    long      auxWidth,
    long*     auxLen,
    double**  aux,
    int       is_end_of_day_fixing,
    int       is_end_of_day_payment,
    long      num_paths,
    int       do_pecs,
    int       do_jump,
    double*** prod_val)
{
    int          free_str = 0;
    FIRSTAllMkts xStr;
    SrtGrfnParam defParm;
    int          forback;
    int          flag = 0;
    long         nstp;

    double *time = NULL, *date = NULL;

    double* sig_dom = NULL;

    double *dom_fwd1 = NULL, *dom_fwd2 = NULL, *dom_exp1 = NULL, *dom_exp2 = NULL, *dom_phi1 = NULL,
           *dom_phi2 = NULL, *dom_phi12 = NULL, *dom_gam1_fwd = NULL, *dom_gam2_fwd = NULL,
           *dom_bond_pay = NULL, *dom_gam1_pay = NULL, *dom_gam2_pay = NULL, ***covar = NULL;

    void**       void_prm = NULL;
    GRFNPARMMC2F grfn_prm;

    long today, spot_date;
    int  i, k, num_col, num_und;

    SrtUndPtr *und_ptr = NULL, dom_und;

    double lam_dom, lam_dom2, tau_dom, alpha_dom, gamma_dom, rho_dom;

    char* dom_yc;
    char* domestic_name;

    double pay_time, df;

    double *sigma_date_dom = NULL, *sigma_dom = NULL;

    long pay_date;

    long sigma_n_dom;
    int  dom_idx;

    clock_t t1, t2;

    Err err = NULL;

    t1 = clock();

    /*	Initialise the GRFN tableau */

    /*	First, initialise the param struct */

    err                        = srt_f_set_default_GrfnParams(&defParm);
    defParm.num_MCarlo_paths   = num_paths;
    defParm.max_time_per_slice = 1000;
    defParm.min_nodes_per_path = 1;
    defParm.force_mc           = 1;
    defParm.jumping            = 1;

    /* End of Day Fixing */
    if (is_end_of_day_fixing)
    {
        defParm.end_of_day_fixings = SRT_YES;
        defParm.end_of_day_flg     = SRT_YES;
    }
    else
    {
        defParm.end_of_day_fixings = SRT_NO;
        defParm.end_of_day_flg     = SRT_NO;
    }

    /* End of Day Payment */
    if (is_end_of_day_payment)
    {
        defParm.end_of_day_payment = SRT_YES;
    }
    else
    {
        defParm.end_of_day_payment = SRT_NO;
    }

    err = FIRSTInitMktStruct(
        numeventdates,
        eventdates,
        tableauRows,
        *tableauCols,
        tableauStrings,
        tableauMask,
        auxWidth,
        auxLen,
        aux,
        lgmund,
        &defParm,
        &forback,
        &xStr);

    if (err)
        goto FREE_RETURN;

    free_str = 1;

    /*	Now, lookup underlyings involved */
    err = FIRSTGetUndFromDeal(&xStr, &num_und, &und_ptr);

    if (err)
    {
        goto FREE_RETURN;
    }

    if (num_und != 1)
    {
        err = "Product should involve only one underlying";
        goto FREE_RETURN;
    }

    dom_und = und_ptr[0];

    /* look for the underlying name */
    dom_und = lookup_und(lgmund);
    if (!dom_und)
    {
        err = "cannot find the underlying";
        goto FREE_RETURN;
    }

    if (get_mdltype_from_irund(dom_und) == LGM)
    {
        domestic_name = dom_und->underl_name;
    }
    else
    {
        err = "Model must be LGM2F";
        goto FREE_RETURN;
    }

    if (strcmp(domestic_name, dom_und->underl_name))
    {
        err = "Tableau uses different underlying";
        goto FREE_RETURN;
    }

    /* look for the today date */
    today     = get_today_from_underlying(dom_und);
    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    num_col = xStr.num_cols;

    /*	Next, get the time steps */

    /*	Copy event dates */
    nstp = xStr.num_evt;
    while (nstp >= 1 && xStr.evt[nstp - 1].evt == NULL)
    {
        nstp--;
    }
    if (nstp < 1)
    {
        err = "No event in Tableau";
        goto FREE_RETURN;
    }

    time = (double*)calloc(nstp, sizeof(double));
    date = (double*)calloc(nstp, sizeof(double));

    if (!time || !date)
    {
        err = "Memory allocation error (1) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    memcpy(time, xStr.tms, nstp * sizeof(double));

    for (i = 0; i < nstp; i++)
    {
        date[i] = today + DAYS_IN_YEAR * time[i];

        if (i > 0 && date[i] - date[i - 1] >= 1)
        {
            date[i] = (long)(date[i] + 1.0e-08);
            time[i] = YEARS_IN_DAY * (date[i] - today);
        }
    }

    if (time[0] > 0)
    {
        /* add the zero time */
        num_f_add_number(&nstp, &time, 0);
        num_f_sort_vector(nstp, time);
        nstp -= 1;
        num_f_add_number(&nstp, &date, today);
        num_f_sort_vector(nstp, date);
        flag = 1;
    }

    if (do_jump)
    {
        pay_date = (long)(date[1] + 1.0E-8);
        pay_time = time[1];
    }
    else
    {
        pay_date = (long)(date[nstp - 1] + 1.0E-8);
        pay_time = time[nstp - 1];
    }

    /* Get all the term structures */
    /* Get the term structure */
    err = Get_LGM2F_TermStructure(
        domestic_name,
        &sigma_date_dom,
        &sigma_dom,
        &sigma_n_dom,
        &tau_dom,
        &alpha_dom,
        &gamma_dom,
        &rho_dom);

    if (err)
    {
        goto FREE_RETURN;
    }

    lam_dom  = 1.0 / tau_dom;
    lam_dom2 = lam_dom + gamma_dom;

    /*	Get Fx spot and yield curves */
    dom_yc = (char*)get_ycname_from_irund(dom_und);

    /*	Get distributions */

    dom_fwd1     = (double*)calloc(nstp, sizeof(double));
    dom_fwd2     = (double*)calloc(nstp, sizeof(double));
    dom_exp1     = (double*)calloc(nstp, sizeof(double));
    dom_exp2     = (double*)calloc(nstp, sizeof(double));
    dom_phi1     = (double*)calloc(nstp, sizeof(double));
    dom_phi2     = (double*)calloc(nstp, sizeof(double));
    dom_phi12    = (double*)calloc(nstp, sizeof(double));
    dom_gam1_fwd = (double*)calloc(nstp, sizeof(double));
    dom_gam2_fwd = (double*)calloc(nstp, sizeof(double));
    dom_bond_pay = (double*)calloc(nstp, sizeof(double));
    dom_gam1_pay = (double*)calloc(nstp, sizeof(double));
    dom_gam2_pay = (double*)calloc(nstp, sizeof(double));

    covar = f3tensor(0, nstp - 1, 0, 1, 0, 1);

    if (!dom_fwd1 || !dom_fwd2 || !dom_phi1 || !dom_phi2 || !dom_phi12 || !dom_exp1 || !dom_exp2 ||
        !dom_gam1_fwd || !dom_gam2_fwd || !dom_bond_pay || !dom_gam1_pay || !dom_gam2_pay || !covar)
    {
        err = "Memory allocation error (3) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    fill_mc_init_lgm2f(
        do_jump,
        pay_date,
        pay_time,
        date,
        time,
        nstp,
        sigma_date_dom,
        sigma_n_dom,
        sigma_dom,
        lam_dom,
        alpha_dom,
        gamma_dom,
        rho_dom,
        dom_yc,
        dom_fwd1,
        dom_fwd2,
        dom_exp1,
        dom_exp2,
        dom_phi1,
        dom_phi2,
        dom_phi12,
        dom_gam1_fwd,
        dom_gam2_fwd,
        dom_bond_pay,
        dom_gam1_pay,
        dom_gam2_pay,
        covar);

    /*	Fill product structure */

    void_prm = (void**)calloc(nstp, sizeof(void*));

    if (!void_prm)
    {
        err = "Memory allocation error (4) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    dom_idx = 0;

    for (i = xStr.num_evt - 1; i >= 0; i--)
    {
        if (xStr.evt[i].evt)
        {
            grfn_prm         = malloc(sizeof(grfn_parm_mc2F));
            grfn_prm->global = &xStr;
            grfn_prm->local  = xStr.evt + i;

            grfn_prm->num_dom_df = xStr.evt[i].evt->dflen[dom_idx];
            grfn_prm->dom_df_tms = xStr.evt[i].evt->dft[dom_idx];
            grfn_prm->dom_df_dts = xStr.evt[i].evt->dfd[dom_idx];

            if (grfn_prm->num_dom_df > 0)
            {
                grfn_prm->dom_dff   = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam1  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam2  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam12 = dvector(0, grfn_prm->num_dom_df - 1);

                if (!grfn_prm->dom_dff || !grfn_prm->dom_gam1 || !grfn_prm->dom_gam2 ||
                    !grfn_prm->dom_gam12)
                {
                    err = "Memory allocation error (8) in SrtGrfn5DFXMc";
                    goto FREE_RETURN;
                }

                for (k = 0; k < grfn_prm->num_dom_df; k++)
                {
                    grfn_prm->dom_dff[k] =
                        swp_f_df(xStr.dts[i], grfn_prm->dom_df_dts[k], (char*)dom_yc);
                    grfn_prm->dom_gam1[k] =
                        (1.0 - exp(-lam_dom * grfn_prm->dom_df_tms[k])) / lam_dom;
                    grfn_prm->dom_gam2[k] =
                        (1.0 - exp(-lam_dom2 * grfn_prm->dom_df_tms[k])) / lam_dom2;
                    grfn_prm->dom_gam12[k] =
                        -0.5 *
                            (grfn_prm->dom_gam1[k] * grfn_prm->dom_gam1[k] * dom_phi1[i + flag] +
                             grfn_prm->dom_gam2[k] * grfn_prm->dom_gam2[k] * dom_phi2[i + flag]) -
                        grfn_prm->dom_gam1[k] * grfn_prm->dom_gam2[k] * dom_phi12[i + flag];
                }

                grfn_prm->do_dom = 1;
            }
            else
            {
                grfn_prm->do_dom = 0;
            }

            void_prm[i + flag] = (void*)grfn_prm;
        }
        else
        {
            void_prm[i + flag] = NULL;
        }
    }

    /*	Eventually! call to function */

    *prod_val = dmatrix(0, num_col - 1, 0, 1);

    if (!(*prod_val))
    {
        err = "Memory allocation error";
        goto FREE_RETURN;
    }

    t2 = clock();

    smessage("Phase 1 -preprocessing, time in sec: %.2f", (double)(t2 - t1) / CLOCKS_PER_SEC);

    err = mc_main_lgm2f(
        /*	Time data */
        num_paths,
        num_col,
        time,
        date,
        nstp,
        do_jump,
        dom_fwd1,
        dom_fwd2,
        dom_exp1,
        dom_exp2,
        dom_phi1,
        dom_phi2,
        dom_phi12,
        dom_gam1_fwd,
        dom_gam2_fwd,
        dom_bond_pay,
        dom_gam1_pay,
        dom_gam2_pay,
        covar,
        void_prm,
        do_pecs,
        0,
        NULL,
        NULL,
        NULL,
        /*	Payoff function */
        grfn_payoff_lgm2f_mc, /*	Result */
        *prod_val);

    if (err)
        goto FREE_RETURN;

    df = swp_f_zr(today, pay_date, dom_yc);
    df = exp(-df * pay_time);

    *tableauCols = num_col;
    for (i = 0; i < num_col; i++)
    {
        (*prod_val)[i][0] *= df;
        (*prod_val)[i][1] *= df;
    }

    /*	Add PV of Past */
    (*prod_val)[num_col - 1][0] += xStr.gd->pv_of_past;

FREE_RETURN:

    if (time)
        free(time);
    if (date)
        free(date);
    if (dom_fwd1)
        free(dom_fwd1);
    if (dom_fwd2)
        free(dom_fwd2);
    if (dom_exp1)
        free(dom_exp1);
    if (dom_exp2)
        free(dom_exp2);
    if (dom_phi1)
        free(dom_phi1);
    if (dom_phi2)
        free(dom_phi2);
    if (dom_phi12)
        free(dom_phi12);
    if (dom_gam1_fwd)
        free(dom_gam1_fwd);
    if (dom_gam2_fwd)
        free(dom_gam2_fwd);
    if (dom_bond_pay)
        free(dom_bond_pay);
    if (dom_gam1_pay)
        free(dom_gam1_pay);
    if (dom_gam2_pay)
        free(dom_gam2_pay);

    if (covar)
        free_f3tensor(covar, 0, nstp - 1, 0, 1, 0, 1);

    if (sigma_date_dom)
        free(sigma_date_dom);
    if (sigma_dom)
        free(sigma_dom);

    if (void_prm)
    {
        for (i = 0; i < nstp; i++)
        {
            if (void_prm[i])
            {
                grfn_prm = (GRFNPARMMC2F)void_prm[i];

                if (grfn_prm->do_dom && grfn_prm->num_dom_df > 0)
                {
                    if (grfn_prm->dom_dff)
                        free_dvector(grfn_prm->dom_dff, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam1)
                        free_dvector(grfn_prm->dom_gam1, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam2)
                        free_dvector(grfn_prm->dom_gam2, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam12)
                        free_dvector(grfn_prm->dom_gam12, 0, grfn_prm->num_dom_df - 1);
                }

                free(grfn_prm);
            }
        }

        free(void_prm);
    }

    if (free_str)
    {
        FIRSTFreeMktStruct(&xStr);
    }

    return err;
}

// LGM2F MCEB with constant lambda
char* SrtGrfnLGM2FMCEB(
    char*      lgmund,
    int        numeventdates,
    long*      eventdates,
    int*       nUsedEventDates,
    int*       optimise,
    double*    fwd_iv,
    MCEBPARAMS params,
    long*      resRows,
    long       tableauRows,
    long*      tableauCols,
    char***    tableauStrings,
    int**      tableauMask,
    long       auxWidth,
    long*      auxLen,
    double**   aux,
    int        is_end_of_day_fixing,
    int        is_end_of_day_payment,
    long       num_paths,
    int        do_pecs,
    int        do_jump,
    double***  prod_val)
{
    int          free_str = 0;
    FIRSTAllMkts xStr;
    SrtGrfnParam defParm;
    int          forback;
    int          flag = 0;
    long         nstp;
    int*         ivOptimise = 0;

    double *time = NULL, *date = NULL;

    double* sig_dom = NULL;

    double *dom_fwd1 = NULL, *dom_fwd2 = NULL, *dom_exp1 = NULL, *dom_exp2 = NULL, *dom_phi1 = NULL,
           *dom_phi2 = NULL, *dom_phi12 = NULL, *dom_gam1_fwd = NULL, *dom_gam2_fwd = NULL,
           *dom_bond_pay = NULL, *dom_gam1_pay = NULL, *dom_gam2_pay = NULL, ***covar = NULL;

    void**       void_prm = NULL;
    GRFNPARMMC2F grfn_prm;

    long today, spot_date;
    int  i, k, num_col, num_und;

    SrtUndPtr *und_ptr = NULL, dom_und;

    double lam_dom, lam_dom2, tau_dom, alpha_dom, gamma_dom, rho_dom;

    char* dom_yc;
    char* domestic_name;

    double pay_time, df;

    double *sigma_date_dom = NULL, *sigma_dom = NULL;

    long pay_date;

    long sigma_n_dom;
    int  dom_idx;

    clock_t t1, t2;

    Err err = NULL;

    t1 = clock();

    /*	Initialise the GRFN tableau */

    /*	First, initialise the param struct */

    err                        = srt_f_set_default_GrfnParams(&defParm);
    defParm.num_MCarlo_paths   = num_paths;
    defParm.max_time_per_slice = 1000;
    defParm.min_nodes_per_path = 1;
    defParm.force_mc           = 1;
    defParm.jumping            = 1;

    /* End of Day Fixing */
    if (is_end_of_day_fixing)
    {
        defParm.end_of_day_fixings = SRT_YES;
        defParm.end_of_day_flg     = SRT_YES;
    }
    else
    {
        defParm.end_of_day_fixings = SRT_NO;
        defParm.end_of_day_flg     = SRT_NO;
    }

    /* End of Day Payment */
    if (is_end_of_day_payment)
    {
        defParm.end_of_day_payment = SRT_YES;
    }
    else
    {
        defParm.end_of_day_payment = SRT_NO;
    }

    err = FIRSTInitMktStruct(
        numeventdates,
        eventdates,
        tableauRows,
        *tableauCols,
        tableauStrings,
        tableauMask,
        auxWidth,
        auxLen,
        aux,
        lgmund,
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

    if (num_und != 1)
    {
        err = "Product should involve only one underlying";
        goto FREE_RETURN;
    }

    dom_und = und_ptr[0];

    /* look for the underlying name */
    dom_und = lookup_und(lgmund);
    if (!dom_und)
    {
        err = "cannot find the underlying";
        goto FREE_RETURN;
    }

    if (get_mdltype_from_irund(dom_und) == LGM)
    {
        domestic_name = dom_und->underl_name;
    }
    else
    {
        err = "Model must be LGM2F";
        goto FREE_RETURN;
    }

    if (strcmp(domestic_name, dom_und->underl_name))
    {
        err = "Tableau uses different underlying";
        goto FREE_RETURN;
    }

    /* look for the today date */
    today     = get_today_from_underlying(dom_und);
    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    num_col = xStr.num_cols;

    /*	Next, get the time steps */

    /*	Copy event dates */
    nstp = xStr.num_evt;
    while (nstp >= 1 && xStr.evt[nstp - 1].evt == NULL)
    {
        nstp--;
    }
    if (nstp < 1)
    {
        err = "No event in Tableau";
        goto FREE_RETURN;
    }

    time = (double*)calloc(nstp, sizeof(double));
    date = (double*)calloc(nstp, sizeof(double));

    if (!time || !date)
    {
        err = "Memory allocation error (1) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    memcpy(time, xStr.tms, nstp * sizeof(double));

    for (i = 0; i < nstp; i++)
    {
        date[i] = today + DAYS_IN_YEAR * time[i];

        if (i > 0 && date[i] - date[i - 1] >= 1)
        {
            date[i] = (long)(date[i] + 1.0e-08);
            time[i] = YEARS_IN_DAY * (date[i] - today);
        }
    }

    if (time[0] > 0)
    {
        /* add the zero time */
        num_f_add_number(&nstp, &time, 0);
        num_f_sort_vector(nstp, time);
        nstp -= 1;
        num_f_add_number(&nstp, &date, today);
        num_f_sort_vector(nstp, date);
        flag = 1;
    }

    if (do_jump)
    {
        pay_date = (long)(date[1] + 1.0E-8);
        pay_time = time[1];
    }
    else
    {
        pay_date = (long)(date[nstp - 1] + 1.0E-8);
        pay_time = time[nstp - 1];
    }

    /* Get all the term structures */
    /* Get the term structure */
    err = Get_LGM2F_TermStructure(
        domestic_name,
        &sigma_date_dom,
        &sigma_dom,
        &sigma_n_dom,
        &tau_dom,
        &alpha_dom,
        &gamma_dom,
        &rho_dom);

    if (err)
    {
        goto FREE_RETURN;
    }

    lam_dom  = 1.0 / tau_dom;
    lam_dom2 = lam_dom + gamma_dom;

    /*	Get Fx spot and yield curves */
    dom_yc = (char*)get_ycname_from_irund(dom_und);

    /*	Get distributions */

    dom_fwd1     = (double*)calloc(nstp, sizeof(double));
    dom_fwd2     = (double*)calloc(nstp, sizeof(double));
    dom_exp1     = (double*)calloc(nstp, sizeof(double));
    dom_exp2     = (double*)calloc(nstp, sizeof(double));
    dom_phi1     = (double*)calloc(nstp, sizeof(double));
    dom_phi2     = (double*)calloc(nstp, sizeof(double));
    dom_phi12    = (double*)calloc(nstp, sizeof(double));
    dom_gam1_fwd = (double*)calloc(nstp, sizeof(double));
    dom_gam2_fwd = (double*)calloc(nstp, sizeof(double));
    dom_bond_pay = (double*)calloc(nstp, sizeof(double));
    dom_gam1_pay = (double*)calloc(nstp, sizeof(double));
    dom_gam2_pay = (double*)calloc(nstp, sizeof(double));

    covar = f3tensor(0, nstp - 1, 0, 1, 0, 1);

    if (!dom_fwd1 || !dom_fwd2 || !dom_phi1 || !dom_phi2 || !dom_phi12 || !dom_exp1 || !dom_exp2 ||
        !dom_gam1_fwd || !dom_gam2_fwd || !dom_bond_pay || !dom_gam1_pay || !dom_gam2_pay || !covar)
    {
        err = "Memory allocation error (3) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    fill_mc_init_lgm2f(
        do_jump,
        pay_date,
        pay_time,
        date,
        time,
        nstp,
        sigma_date_dom,
        sigma_n_dom,
        sigma_dom,
        lam_dom,
        alpha_dom,
        gamma_dom,
        rho_dom,
        dom_yc,
        dom_fwd1,
        dom_fwd2,
        dom_exp1,
        dom_exp2,
        dom_phi1,
        dom_phi2,
        dom_phi12,
        dom_gam1_fwd,
        dom_gam2_fwd,
        dom_bond_pay,
        dom_gam1_pay,
        dom_gam2_pay,
        covar);

    /*	Fill product structure */

    void_prm = (void**)calloc(nstp, sizeof(void*));

    if (!void_prm)
    {
        err = "Memory allocation error (4) in SrtGrfn5DFXMc";
        goto FREE_RETURN;
    }

    dom_idx = 0;

    for (i = xStr.num_evt - 1; i >= 0; i--)
    {
        if (xStr.evt[i].evt)
        {
            grfn_prm         = malloc(sizeof(grfn_parm_mc2F));
            grfn_prm->global = &xStr;
            grfn_prm->local  = xStr.evt + i;

            grfn_prm->num_dom_df = xStr.evt[i].evt->dflen[dom_idx];
            grfn_prm->dom_df_tms = xStr.evt[i].evt->dft[dom_idx];
            grfn_prm->dom_df_dts = xStr.evt[i].evt->dfd[dom_idx];

            if (grfn_prm->num_dom_df > 0)
            {
                grfn_prm->dom_dff   = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam1  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam2  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam12 = dvector(0, grfn_prm->num_dom_df - 1);

                if (!grfn_prm->dom_dff || !grfn_prm->dom_gam1 || !grfn_prm->dom_gam2 ||
                    !grfn_prm->dom_gam12)
                {
                    err = "Memory allocation error (8) in SrtGrfn5DFXMc";
                    goto FREE_RETURN;
                }

                for (k = 0; k < grfn_prm->num_dom_df; k++)
                {
                    grfn_prm->dom_dff[k] =
                        swp_f_df(xStr.dts[i], grfn_prm->dom_df_dts[k], (char*)dom_yc);
                    grfn_prm->dom_gam1[k] =
                        (1.0 - exp(-lam_dom * grfn_prm->dom_df_tms[k])) / lam_dom;
                    grfn_prm->dom_gam2[k] =
                        (1.0 - exp(-lam_dom2 * grfn_prm->dom_df_tms[k])) / lam_dom2;
                    grfn_prm->dom_gam12[k] =
                        -0.5 *
                            (grfn_prm->dom_gam1[k] * grfn_prm->dom_gam1[k] * dom_phi1[i + flag] +
                             grfn_prm->dom_gam2[k] * grfn_prm->dom_gam2[k] * dom_phi2[i + flag]) -
                        grfn_prm->dom_gam1[k] * grfn_prm->dom_gam2[k] * dom_phi12[i + flag];
                }

                grfn_prm->do_dom = 1;
            }
            else
            {
                grfn_prm->do_dom = 0;
            }

            void_prm[i + flag] = (void*)grfn_prm;
        }
        else
        {
            void_prm[i + flag] = NULL;
        }
    }

    /*	Eventually! call to function */

    *tableauCols     = num_col;
    *resRows         = max(num_col + 1, nstp);
    *nUsedEventDates = xStr.num_evt;
    df               = swp_f_df(today, pay_date, dom_yc);

    /* create an optimisation vector from the input */
    ivOptimise = ivector(0, nstp - 1);
    if (flag)
        ivOptimise[0] = 0;
    for (i = flag; i < nstp; i++)
        ivOptimise[i] = optimise[i + numeventdates - xStr.num_evt - flag];

    if (params->iMultiIndex)
    {
        params->iNbIndex = params->iMultiIndex;
    }
    else
    {
        params->iNbIndex = 1;
    }

    mceb_allocate_params(params, nstp);

    if (params->iAdjustIV)
    {
        if (flag)
        {
            params->dMarketFwdIV[0] = 0.0;
        }

        for (i = flag; i < nstp; i++)
        {
            params->dMarketFwdIV[i] = fwd_iv[i + numeventdates - xStr.num_evt - flag] / df;
        }
    }

    *prod_val = dmatrix(0, *resRows - 1, 0, 2 + params->iNbIndex);

    if (!(*prod_val))
    {
        err = "Memory allocation error";
        goto FREE_RETURN;
    }

    t2 = clock();

    smessage("Phase 1 -preprocessing, time in sec: %.2f", (double)(t2 - t1) / CLOCKS_PER_SEC);

    err = mc_main_lgm2f(
        /*	Time data */
        num_paths,
        num_col,
        time,
        date,
        nstp,
        do_jump,
        dom_fwd1,
        dom_fwd2,
        dom_exp1,
        dom_exp2,
        dom_phi1,
        dom_phi2,
        dom_phi12,
        dom_gam1_fwd,
        dom_gam2_fwd,
        dom_bond_pay,
        dom_gam1_pay,
        dom_gam2_pay,
        covar,
        void_prm,
        do_pecs,
        1,
        ivOptimise,
        params,
        NULL,
        /*	Payoff function */
        grfn_payoff_lgm2f_mc, /*	Result */
        *prod_val);

    if (err)
        goto FREE_RETURN;

    /* Recopy Barrier / CoefLin for the moment */
    for (i = 0; i < nstp; i++)
    {
        (*prod_val)[i][2] = params->dBarrier[i];

        for (k = 0; k < params->iNbIndex; k++)
        {
            (*prod_val)[i][3 + k] = params->dCoefLin[i][k + 1];
        }
    }

    for (i = 0; i < num_col + 1; i++)
    {
        (*prod_val)[i][0] *= df;
        (*prod_val)[i][1] *= df;
    }

    for (i = 0; i < nstp; i++)
    {
        if (params->iCalcOneTime)
            params->dOneTimeCall[i] *= df;
        if (params->iCalcOneTimePartial)
            params->dOneTimePartial[i] *= df;
        if (params->iCalcIV || params->iAdjustIV)
            params->dModelFwdIV[i] *= df;
    }

    if (flag)
    {
        for (i = 0; i < nstp - 1; i++)
        {
            (*prod_val)[i][2] = (*prod_val)[i + 1][2];

            for (k = 0; k < params->iNbIndex; k++)
            {
                (*prod_val)[i][3 + k] = (*prod_val)[i + 1][3 + k];
            }
        }

        mceb_shift_extrainfos(params);
    }

    /*	Add PV of Past */
    (*prod_val)[num_col - 1][0] += xStr.gd->pv_of_past;
    (*prod_val)[num_col][0] += xStr.gd->pv_of_past;

FREE_RETURN:

    if (time)
        free(time);
    if (date)
        free(date);
    if (dom_fwd1)
        free(dom_fwd1);
    if (dom_fwd2)
        free(dom_fwd2);
    if (dom_exp1)
        free(dom_exp1);
    if (dom_exp2)
        free(dom_exp2);
    if (dom_phi1)
        free(dom_phi1);
    if (dom_phi2)
        free(dom_phi2);
    if (dom_phi12)
        free(dom_phi12);
    if (dom_gam1_fwd)
        free(dom_gam1_fwd);
    if (dom_gam2_fwd)
        free(dom_gam2_fwd);
    if (dom_bond_pay)
        free(dom_bond_pay);
    if (dom_gam1_pay)
        free(dom_gam1_pay);
    if (dom_gam2_pay)
        free(dom_gam2_pay);

    if (covar)
        free_f3tensor(covar, 0, nstp - 1, 0, 1, 0, 1);

    if (sigma_date_dom)
        free(sigma_date_dom);
    if (sigma_dom)
        free(sigma_dom);
    if (ivOptimise)
        free_ivector(ivOptimise, 0, nstp - 1);

    if (void_prm)
    {
        for (i = 0; i < nstp; i++)
        {
            if (void_prm[i])
            {
                grfn_prm = (GRFNPARMMC2F)void_prm[i];

                if (grfn_prm->do_dom && grfn_prm->num_dom_df > 0)
                {
                    if (grfn_prm->dom_dff)
                        free_dvector(grfn_prm->dom_dff, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam1)
                        free_dvector(grfn_prm->dom_gam1, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam2)
                        free_dvector(grfn_prm->dom_gam2, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam12)
                        free_dvector(grfn_prm->dom_gam12, 0, grfn_prm->num_dom_df - 1);
                }

                free(grfn_prm);
            }
        }

        free(void_prm);
    }

    if (free_str)
    {
        FIRSTFreeMktStruct(&xStr);
    }

    return err;
}

//**********************************************************************
//*                                                                    *
//*               LGM2F MC with Tau TS (J.M.L Nov 2003)                *
//*                                                                    *
//**********************************************************************

// Merges the term structures of lambda and sigma
Err merge_lambda_sigma_ts(
    double*  sigma,
    double*  sigma_time,
    int      nb_sigma,
    double*  lambda,
    double*  lambda_time,
    int      nb_lambda,
    double** ts_time,
    double** lam,
    double** sig,
    int*     nb_new_time)
{
    int i;
    int nb_ts;

    Err err = NULL;

    *ts_time = NULL;
    *lam     = NULL;
    *sig     = NULL;

    /* First create a target TS */

    (*ts_time) = (double*)calloc(nb_sigma, sizeof(double));

    if (!(*ts_time))
    {
        err = "Memory allocation failure in merge_lambda_sigma_ts";
        goto FREE_RETURN;
    }

    memcpy((*ts_time), sigma_time, nb_sigma * sizeof(double));
    nb_ts = nb_sigma;

    for (i = 0; i < nb_lambda; i++)
    {
        num_f_add_number(&nb_ts, ts_time, lambda_time[i]);
    }

    num_f_sort_vector(nb_ts, *ts_time);
    num_f_unique_vector(&nb_ts, *ts_time);

    (*sig) = (double*)calloc(nb_ts, sizeof(double));
    (*lam) = (double*)calloc(nb_ts, sizeof(double));

    if (!(*sig) || !(*lam))
    {
        err = "Memory allocation failure in merge_lambda_sigma_ts";
        goto FREE_RETURN;
    }

    for (i = 0; i < nb_ts; i++)
    {
        (*sig)[i] = sigma[Get_Index((*ts_time)[i], sigma_time, nb_sigma)];
        (*lam)[i] = lambda[Get_Index((*ts_time)[i], lambda_time, nb_lambda)];
    }

    (*nb_new_time) = nb_ts;

FREE_RETURN:
    if (err)
    {
        if (*ts_time)
            free(*ts_time);
        if (*sig)
            free(*sig);
        if (*lam)
            free(*lam);
    }

    return err;

    return NULL;
}

// Calculates exp_lambda = exp(-integral(s,T,lambda(w),w))
// when lambda(w) is piecewise constant
Err exp_lambda(
    double s, double T, double* lambda_time, double* lambda, long n_lambda, double* result)
{
    long   i, StartIndex, EndIndex;
    double T1, T2;
    double temp = 0.0, sum = 0.0;

    StartIndex = Get_Index(s, lambda_time, n_lambda);
    EndIndex   = Get_Index(T, lambda_time, n_lambda);

    for (i = StartIndex; i < EndIndex + 1; i++)
    {
        if (i > StartIndex)
        {
            T1 = lambda_time[i - 1];
        }
        else
        {
            // First part
            T1 = s;
        }

        if (i == EndIndex || StartIndex == EndIndex)
        {
            // Last part
            T2 = T;
        }
        else
        {
            T2 = lambda_time[i];
        }

        temp = lambda[i] * (T2 - T1);

        sum += temp;
    }
    *result = exp(-sum);
    return NULL;
}

// Calculates gamma_lambda = integral(s,T,exp(-integral(s,u,lambda(w),w),du
// when lambda(w) is piecewise constant
Err gamma_lambda(
    double s, double T, double* lambda_time, double* lambda, long n_lambda, double* result)
{
    long   i, StartIndex, EndIndex;
    double T1, T2;
    double aux = 1.0, temp = 0.0, sum = 0.0;

    StartIndex = Get_Index(s, lambda_time, n_lambda);
    EndIndex   = Get_Index(T, lambda_time, n_lambda);

    for (i = StartIndex; i < EndIndex + 1; i++)
    {
        if (i > StartIndex)
        {
            T1 = lambda_time[i - 1];
        }
        else
        {
            // First part
            T1 = s;
        }

        if (i == EndIndex || StartIndex == EndIndex)
        {
            // Last part
            T2 = T;
        }
        else
        {
            T2 = lambda_time[i];
        }

        temp = (1 - exp(-lambda[i] * (T2 - T1))) / lambda[i];

        sum += aux * temp;

        aux = aux * exp(-lambda[i] * (T2 - T1));
    }
    *result = sum;
    return NULL;
}

// Calculates intermediary results for the calculation of the forwards
Err calculate_Bk_rhoCk(
    double  Tk,
    double  Tk1,
    long    k,
    double  Tn1,
    double* sigma_dates,
    long    nb_sigma_dates,
    double* sigma_dom,
    double* lda_dom,
    double* lda_dom2,
    double  alpha_dom,
    double  gamma_dom,
    double  rho_dom,
    double* Bk_rhoCk_fwd1,
    double* Bk_rhoCk_fwd2)
{
    Err    err;
    long   i;
    double T1i, T2i, a1, a2, alpha_dom2, sigma_dom2, sigma_dom_k;
    double sumBk1, sumCk1, Bk1, Ck1;
    double sumBk2, sumCk2, Bk2, Ck2;
    double phi_minus_lda_dom, phi_lda_dom, phi_lda_dom2, phi_minus_lda_dom2;

    // Initialization
    phi_minus_lda_dom  = Phi_Func(-lda_dom[k], Tk, Tk, Tk1);
    phi_lda_dom        = Phi_Func(lda_dom[k], Tk, Tk, Tk1);
    phi_minus_lda_dom2 = Phi_Func(-lda_dom2[k], Tk, Tk, Tk1);
    phi_lda_dom2       = Phi_Func(lda_dom2[k], Tk, Tk, Tk1);

    Bk1 = Ck1 = sumBk1 = sumCk1 = 0.0;
    Bk2 = Ck2 = sumBk2 = sumCk2 = 0.0;

    alpha_dom2 = alpha_dom * alpha_dom;

    // Between [0, k] inclusive
    for (i = 0; i < k + 1; i++)
    {
        if (i == 0)
        {
            T1i = 0.0;
        }
        else
        {
            // First part
            T1i = sigma_dates[i - 1];
        }

        if (i == k)
        {
            // Last part
            T2i = Tk;
        }
        else
        {
            T2i = sigma_dates[i];
        }

        err = exp_lambda(T2i, Tk, sigma_dates, lda_dom, nb_sigma_dates, &a1);
        err = exp_lambda(T2i, Tk, sigma_dates, lda_dom2, nb_sigma_dates, &a2);

        sigma_dom2 = sigma_dom[i] * sigma_dom[i];

        // To be used for Bk1 and Ck1
        sumBk1 +=
            sigma_dom2 / (2 * lda_dom[i]) * a1 * a1 * (1 - exp(-2 * lda_dom[i] * (T2i - T1i)));
        sumCk1 += alpha_dom * sigma_dom2 / (lda_dom[i] + lda_dom2[i]) * a1 * a2 *
                  (1 - exp(-(lda_dom[i] + lda_dom2[i]) * (T2i - T1i)));

        // To be used for Bk2 (the sum for Ck2 is identical to Ck1)
        sumBk2 += alpha_dom2 * sigma_dom2 / (2 * lda_dom2[i]) * a2 * a2 *
                  (1 - exp(-2 * lda_dom2[i] * (T2i - T1i)));
    }

    // On [k, k+1]
    sigma_dom_k = sigma_dom[k];
    err         = exp_lambda(Tk, Tn1, sigma_dates, lda_dom, nb_sigma_dates, &a1);
    err         = exp_lambda(Tk, Tn1, sigma_dates, lda_dom2, nb_sigma_dates, &a2);

    // For fwd1
    Bk1 = sigma_dom_k * sigma_dom_k / (2 * lda_dom[k]) * (phi_lda_dom - phi_minus_lda_dom);
    Bk1 += phi_minus_lda_dom * sumBk1;
    Bk1 *= a1;

    Ck1 = alpha_dom * sigma_dom_k * sigma_dom_k / (lda_dom[k] + lda_dom2[k]) *
          (phi_lda_dom - phi_minus_lda_dom2);
    Ck1 += phi_minus_lda_dom2 * sumCk1;
    Ck1 *= a1;

    // For fwd2
    Bk2 = alpha_dom2 * sigma_dom_k * sigma_dom_k / (2 * lda_dom2[k]) *
          (phi_lda_dom2 - phi_minus_lda_dom2);
    Bk2 += phi_minus_lda_dom2 * sumBk2;
    Bk2 *= a2;

    Ck2 = alpha_dom * sigma_dom_k * sigma_dom_k / (lda_dom[k] + lda_dom2[k]) *
          (phi_lda_dom2 - phi_minus_lda_dom);
    Ck2 += phi_minus_lda_dom * sumCk1;
    Ck2 *= a2;

    // Final results
    *Bk_rhoCk_fwd1 = Bk1 + rho_dom * Ck1;
    *Bk_rhoCk_fwd2 = Bk2 + rho_dom * Ck2;

    return NULL;
}

// Calculates intermediary results for the calculation of the forwards
Err calculate_Dk_rhoEk(
    double  Tk,
    double  Tk1,
    long    k,
    double  Tpay,
    long    iTpay,
    double  Tn1,
    double* sigma_dates,
    long    nb_sigma_dates,
    double* sigma_dom,
    double* lda_dom,
    double* lda_dom2,
    double  alpha_dom,
    double  gamma_dom,
    double  rho_dom,
    double* Dk_rhoEk_fwd1,
    double* Dk_rhoEk_fwd2)
{
    Err    err;
    long   i;
    double Dk1, Ek1, sumDk1, sumEk1, Dk2, Ek2, T1i, T2i, a1, a2, sigma_dom_k;
    double phi_lda_dom, phi_2_lda_dom, phi_lda_dom_lda_dom2, phi_lda_dom2, phi_2_lda_dom2;

    // Initialization
    phi_lda_dom          = Phi_Func(lda_dom[k], Tk1, Tk, Tk1);
    phi_2_lda_dom        = Phi_Func(2 * lda_dom[k], Tk1, Tk, Tk1);
    phi_lda_dom_lda_dom2 = Phi_Func(lda_dom[k] + lda_dom2[k], Tk1, Tk, Tk1);

    phi_lda_dom2   = Phi_Func(lda_dom2[k], Tk1, Tk, Tk1);
    phi_2_lda_dom2 = Phi_Func(2 * lda_dom2[k], Tk1, Tk, Tk1);

    Dk1 = Ek1 = sumDk1 = sumEk1 = Dk2 = Ek2 = 0.0;

    // On [k+1, Tpay]
    for (i = k; i < iTpay + 1; i++)
    {
        if (i == k)
        {
            // Last part
            T1i = Tk1;
        }
        else
        {
            T1i = sigma_dates[i - 1];
        }

        if (i == iTpay)
        {
            // Last part
            T2i = Tpay;
        }
        else
        {
            T2i = sigma_dates[i];
        }

        err = exp_lambda(Tk1, T1i, sigma_dates, lda_dom, nb_sigma_dates, &a1);
        err = exp_lambda(Tk1, T1i, sigma_dates, lda_dom2, nb_sigma_dates, &a2);

        // The sums are common to both forwards
        sumDk1 += a1 * (1 - exp(-lda_dom[i] * (T2i - T1i))) / (lda_dom[i]);
        sumEk1 += a2 * (1 - exp(-lda_dom2[i] * (T2i - T1i))) / (lda_dom2[i]);
    }

    // On [k, k+1]
    sigma_dom_k = sigma_dom[k];
    err         = exp_lambda(Tk1, Tn1, sigma_dates, lda_dom, nb_sigma_dates, &a1);
    err         = exp_lambda(Tk1, Tn1, sigma_dates, lda_dom2, nb_sigma_dates, &a2);

    // For fwd1
    Dk1 = (phi_lda_dom - phi_2_lda_dom) / (lda_dom[k]);
    Dk1 += phi_2_lda_dom * sumDk1;
    Dk1 *= sigma_dom_k * sigma_dom_k * a1;

    Ek1 = (phi_lda_dom - phi_lda_dom_lda_dom2) / (lda_dom2[k]);
    Ek1 += phi_lda_dom_lda_dom2 * sumEk1;
    Ek1 *= sigma_dom_k * sigma_dom_k * alpha_dom * a1;

    // For fwd2
    Dk2 = (phi_lda_dom2 - phi_2_lda_dom2) / (lda_dom2[k]);
    Dk2 += phi_2_lda_dom2 * sumEk1;
    Dk2 *= alpha_dom * alpha_dom * sigma_dom_k * sigma_dom_k * a2;

    Ek2 = (phi_lda_dom2 - phi_lda_dom_lda_dom2) / (lda_dom[k]);
    Ek2 += phi_lda_dom_lda_dom2 * sumDk1;
    Ek2 *= sigma_dom_k * sigma_dom_k * alpha_dom * a2;

    *Dk_rhoEk_fwd1 = Dk1 + rho_dom * Ek1;
    *Dk_rhoEk_fwd2 = Dk2 + rho_dom * Ek2;

    return NULL;
}

Err fill_mc_init_lgm2f_lambda(
    int       do_jump,
    long      pay_date,
    double    pay_time,
    double*   date,
    double*   time,
    long      nb_dates,
    double*   sig_dates,
    long      nb_sig_dates,
    double*   sig_curve_dom,
    double*   lda_dom,
    double*   lda_dom2,
    double    alpha_dom,
    double    gamma_dom,
    double    rho_dom,
    char*     dom_yc,
    double*   dom_fwd1,
    double*   dom_fwd2,
    double*   dom_exp1,
    double*   dom_exp2,
    double*   dom_phi1,
    double*   dom_phi2,
    double*   dom_phi12,
    double*   dom_gam1_fwd,
    double*   dom_gam2_fwd,
    double*   dom_bond_pay,
    double*   dom_gam1_pay,
    double*   dom_gam2_pay,
    double*** covar)
{
    double alpha_dom2;
    double sig_dom, sig_dom2;
    double T1, T2, start_date, end_date, start_mat, end_mat;

    double phi_dom1, phi_dom2, phi_dom12;
    double mat_pay, mat, pay_mat, pay_mat2;
    double zc_dom, zc_pay;
    double exp_dom_mat, exp_dom_mat2, exp_dom_pay, exp_dom_pay2;
    int    k, n;
    long   StartIndex, EndIndex, pay_time_Index;

    double aux_exp_lambda1, aux_exp_lambda2, var1_tn1_tn, var2_tn1_tn, var12_tn1_tn;

    double Bk_rhoCk_fwd1, Dk_rhoEk_fwd1, Bn_rhoCn_fwd1, Dn_rhoEn_fwd1, Bk_rhoCk_fwd2, Dk_rhoEk_fwd2,
        Bn_rhoCn_fwd2, Dn_rhoEn_fwd2;

    double** cov;

    Err err = NULL;

    alpha_dom2 = alpha_dom * alpha_dom;

    dom_fwd1[0] = 0.0;
    dom_fwd2[0] = 0.0;

    dom_phi1[0]  = 0.0;
    dom_phi2[0]  = 0.0;
    dom_phi12[0] = 0.0;

    mat_pay = (pay_date - date[0]) / 365.0;

    start_date = date[0];
    start_mat  = time[0];
    StartIndex = Get_Index(start_mat, sig_dates, nb_sig_dates);

    // For all GRFN event dates
    for (n = 0; n < nb_dates - 1; n++)
    {
        end_date = date[n + 1];

        if (do_jump)
        {
            pay_date = (long)(date[n + 1] + 0.5);
            pay_time = time[n + 1];
        }

        end_mat        = time[n + 1];
        EndIndex       = Get_Index(end_mat, sig_dates, nb_sig_dates);
        pay_time_Index = Get_Index(pay_time, sig_dates, nb_sig_dates);
        mat            = (end_mat - start_mat);
        pay_mat        = (pay_date - start_date) / 365.0;
        pay_mat2       = (pay_date - end_date) / 365.0;

        zc_dom = swp_f_zr(start_date, end_date, dom_yc);

        // Gamma function for reconstruction formula
        err = gamma_lambda(start_mat, end_mat, sig_dates, lda_dom, nb_sig_dates, &dom_gam1_fwd[n]);
        err = gamma_lambda(start_mat, end_mat, sig_dates, lda_dom2, nb_sig_dates, &dom_gam2_fwd[n]);

        // For Df to pay date
        zc_pay = swp_f_zr(start_date, pay_date, dom_yc);

        err = gamma_lambda(start_mat, pay_time, sig_dates, lda_dom, nb_sig_dates, &dom_gam1_pay[n]);
        err =
            gamma_lambda(start_mat, pay_time, sig_dates, lda_dom2, nb_sig_dates, &dom_gam2_pay[n]);

        dom_bond_pay[n] =
            exp(-zc_pay * pay_mat -
                0.5 * (dom_gam1_pay[n] * dom_gam1_pay[n] * dom_phi1[n] +
                       dom_gam2_pay[n] * dom_gam2_pay[n] * dom_phi2[n]) -
                dom_gam1_pay[n] * dom_gam2_pay[n] * dom_phi12[n]);

        //	Variables initialisation
        err = exp_lambda(start_mat, end_mat, sig_dates, lda_dom, nb_sig_dates, &exp_dom_mat);
        err = exp_lambda(start_mat, end_mat, sig_dates, lda_dom2, nb_sig_dates, &exp_dom_mat2);
        err = exp_lambda(end_mat, pay_time, sig_dates, lda_dom, nb_sig_dates, &exp_dom_pay);
        err = exp_lambda(end_mat, pay_time, sig_dates, lda_dom2, nb_sig_dates, &exp_dom_pay2);

        phi_dom1 = phi_dom2 = phi_dom12 = 0.0;
        var1_tn1_tn = var2_tn1_tn = var12_tn1_tn = 0.0;
        Bn_rhoCn_fwd1 = Dn_rhoEn_fwd1 = Bn_rhoCn_fwd2 = Dn_rhoEn_fwd2 = 0.0;

        // For all [Tk,Tk+1] between both event dates Tn and Tn+1 where all params are constants
        for (k = StartIndex; k < EndIndex + 1; k++)
        {
            if (k > StartIndex)
            {
                T1 = sig_dates[k - 1];
            }
            else
            {
                // First part
                T1 = start_mat;
            }

            if (k == EndIndex || StartIndex == EndIndex)
            {
                // Last part
                T2 = end_mat;
            }
            else
            {
                T2 = sig_dates[k];
            }

            sig_dom = sig_dom2 = sig_curve_dom[k];
            sig_dom2 *= sig_dom2;

            // Conditional variance var1_tn1_tn, var2_tn1_tn, var12_tn1_tn for phi1, phi2, phi12
            err = exp_lambda(T2, end_mat, sig_dates, lda_dom, nb_sig_dates, &aux_exp_lambda1);
            var1_tn1_tn +=
                sig_dom2 * aux_exp_lambda1 * aux_exp_lambda1 * Phi_Func(2 * lda_dom[k], T2, T1, T2);

            err = exp_lambda(T2, end_mat, sig_dates, lda_dom2, nb_sig_dates, &aux_exp_lambda2);
            var2_tn1_tn += sig_dom2 * aux_exp_lambda2 * aux_exp_lambda2 *
                           Phi_Func(2 * lda_dom2[k], T2, T1, T2);

            var12_tn1_tn += sig_dom2 * aux_exp_lambda1 * aux_exp_lambda2 *
                            Phi_Func(lda_dom[k] + lda_dom2[k], T2, T1, T2);

            // Terms for the forward reconstruction
            err = calculate_Bk_rhoCk(
                T1,
                T2,
                k,
                end_mat,
                sig_dates,
                nb_sig_dates,
                sig_curve_dom,
                lda_dom,
                lda_dom2,
                alpha_dom,
                gamma_dom,
                rho_dom,
                &Bk_rhoCk_fwd1,
                &Bk_rhoCk_fwd2);

            Bn_rhoCn_fwd1 += Bk_rhoCk_fwd1;
            Bn_rhoCn_fwd2 += Bk_rhoCk_fwd2;

            err = calculate_Dk_rhoEk(
                T1,
                T2,
                k,
                pay_time,
                pay_time_Index,
                end_mat,
                sig_dates,
                nb_sig_dates,
                sig_curve_dom,
                lda_dom,
                lda_dom2,
                alpha_dom,
                gamma_dom,
                rho_dom,
                &Dk_rhoEk_fwd1,
                &Dk_rhoEk_fwd2);

            Dn_rhoEn_fwd1 += Dk_rhoEk_fwd1;
            Dn_rhoEn_fwd2 += Dk_rhoEk_fwd2;
        }

        dom_exp1[n + 1] = exp_dom_mat;
        dom_exp2[n + 1] = exp_dom_mat2;

        dom_phi1[n + 1] = dom_phi1[n] * exp_dom_mat * exp_dom_mat + var1_tn1_tn;
        dom_phi2[n + 1] = dom_phi2[n] * exp_dom_mat2 * exp_dom_mat2 + alpha_dom2 * var2_tn1_tn;
        dom_phi12[n + 1] =
            dom_phi12[n] * exp_dom_mat * exp_dom_mat2 + rho_dom * alpha_dom * var12_tn1_tn;

        // Expectations of the state variables under the QTpay_time probability measure
        dom_fwd1[n + 1] = Bn_rhoCn_fwd1 - Dn_rhoEn_fwd1;
        dom_fwd2[n + 1] = Bn_rhoCn_fwd2 - Dn_rhoEn_fwd2;

        cov = covar[n + 1];

        cov[0][0] = var1_tn1_tn;
        cov[0][1] = cov[1][0] = rho_dom * alpha_dom * var12_tn1_tn;
        cov[1][1]             = alpha_dom2 * var2_tn1_tn;
        cov[1][0]             = cov[0][1];

        start_date = end_date;
        start_mat  = end_mat;
        StartIndex = EndIndex;
    }

    dom_gam1_pay[nb_dates - 1] = 0.0;
    dom_gam2_pay[nb_dates - 1] = 0.0;
    dom_bond_pay[nb_dates - 1] = 1.0;

    return err;
}

// LGM2F MC with lambda term structure
char* SrtGrfnLGM2FMClambda(
    char*     lgmund,
    int       numeventdates,
    long*     eventdates,
    long      tableauRows,
    long*     tableauCols,
    char***   tableauStrings,
    int**     tableauMask,
    long      auxWidth,
    long*     auxLen,
    double**  aux,
    int       is_end_of_day_fixing,
    int       is_end_of_day_payment,
    long      num_paths,
    int       do_pecs,
    int       do_jump,
    double*** prod_val)
{
    int          free_str = 0;
    FIRSTAllMkts xStr;
    SrtGrfnParam defParm;
    int          forback;
    int          flag = 0;
    long         nstp;

    double *time = NULL, *date = NULL;

    double* sig_dom = NULL;

    double *dom_fwd1 = NULL, *dom_fwd2 = NULL, *dom_exp1 = NULL, *dom_exp2 = NULL, *dom_phi1 = NULL,
           *dom_phi2 = NULL, *dom_phi12 = NULL, *dom_gam1_fwd = NULL, *dom_gam2_fwd = NULL,
           *dom_bond_pay = NULL, *dom_gam1_pay = NULL, *dom_gam2_pay = NULL, ***covar = NULL;

    // New variables for lambda TS
    double *new_time = NULL, *new_lambda = NULL, *new_sigma = NULL, *new_lambda2 = NULL;

    long n_new_time;

    double *lam_dom = NULL, *lam_dom_time = NULL;

    void**       void_prm = NULL;
    GRFNPARMMC2F grfn_prm;

    long today, spot_date;
    int  i, k, num_col, num_und;

    SrtUndPtr *und_ptr = NULL, dom_und;

    double alpha_dom, gamma_dom, rho_dom;

    char* dom_yc;
    char* domestic_name;

    double pay_time, df;

    double *sigma_date_dom = NULL, *sigma_dom = NULL;

    long pay_date;

    long n_sigma_dom, n_lam_dom;
    int  dom_idx;

    clock_t t1, t2;

    Err err = NULL;

    t1 = clock();

    /*	Initialise the GRFN tableau */

    /*	First, initialise the param struct */

    err                        = srt_f_set_default_GrfnParams(&defParm);
    defParm.num_MCarlo_paths   = num_paths;
    defParm.max_time_per_slice = 1000;
    defParm.min_nodes_per_path = 1;
    defParm.force_mc           = 1;
    defParm.jumping            = 1;

    /* End of Day Fixing */
    if (is_end_of_day_fixing)
    {
        defParm.end_of_day_fixings = SRT_YES;
        defParm.end_of_day_flg     = SRT_YES;
    }
    else
    {
        defParm.end_of_day_fixings = SRT_NO;
        defParm.end_of_day_flg     = SRT_NO;
    }

    /* End of Day Payment */
    if (is_end_of_day_payment)
    {
        defParm.end_of_day_payment = SRT_YES;
    }
    else
    {
        defParm.end_of_day_payment = SRT_NO;
    }

    err = FIRSTInitMktStruct(
        numeventdates,
        eventdates,
        tableauRows,
        *tableauCols,
        tableauStrings,
        tableauMask,
        auxWidth,
        auxLen,
        aux,
        lgmund,
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

    if (num_und != 1)
    {
        err = "Product should involve only one underlying";
        goto FREE_RETURN;
    }

    dom_und = und_ptr[0];

    /* look for the underlying name */
    dom_und = lookup_und(lgmund);
    if (!dom_und)
    {
        err = "cannot find the underlying";
        goto FREE_RETURN;
    }

    if (get_mdltype_from_irund(dom_und) == LGM)
    {
        domestic_name = dom_und->underl_name;
    }
    else
    {
        err = "Model must be LGM2F";
        goto FREE_RETURN;
    }

    if (strcmp(domestic_name, dom_und->underl_name))
    {
        err = "Tableau uses different underlying";
        goto FREE_RETURN;
    }

    /* look for the today date */
    today     = get_today_from_underlying(dom_und);
    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    num_col = xStr.num_cols;

    /*	Next, get the time steps */

    /*	Copy event dates */
    nstp = xStr.num_evt;
    while (nstp >= 1 && xStr.evt[nstp - 1].evt == NULL)
    {
        nstp--;
    }
    if (nstp < 1)
    {
        err = "No event in Tableau";
        goto FREE_RETURN;
    }

    time = (double*)calloc(nstp, sizeof(double));
    date = (double*)calloc(nstp, sizeof(double));

    if (!time || !date)
    {
        err = "Memory allocation error (1) in SrtGrfnLGM2FMClambda";
        goto FREE_RETURN;
    }

    memcpy(time, xStr.tms, nstp * sizeof(double));

    for (i = 0; i < nstp; i++)
    {
        date[i] = today + DAYS_IN_YEAR * time[i];

        if (i > 0 && date[i] - date[i - 1] >= 1)
        {
            date[i] = (long)(date[i] + 1.0e-08);
            time[i] = YEARS_IN_DAY * (date[i] - today);
        }
    }

    if (time[0] > 0)
    {
        /* add the zero time */
        num_f_add_number(&nstp, &time, 0);
        num_f_sort_vector(nstp, time);
        nstp -= 1;
        num_f_add_number(&nstp, &date, today);
        num_f_sort_vector(nstp, date);
        flag = 1;
    }

    if (do_jump)
    {
        pay_date = (long)(date[1] + 1.0E-8);
        pay_time = time[1];
    }
    else
    {
        pay_date = (long)(date[nstp - 1] + 1.0E-8);
        pay_time = time[nstp - 1];
    }

    // Get the term structure with lambda term structure
    err = Get_LGM2F_TermStructure2(
        domestic_name,
        &sigma_dom,
        &sigma_date_dom,
        &n_sigma_dom,
        &lam_dom,
        &lam_dom_time,
        &n_lam_dom,
        &alpha_dom,
        &gamma_dom,
        &rho_dom);

    if (err)
    {
        goto FREE_RETURN;
    }

    // Merge the term structures of lambda and sigmas
    err = merge_lambda_sigma_ts(
        sigma_dom,
        sigma_date_dom,
        n_sigma_dom,
        lam_dom,
        lam_dom_time,
        n_lam_dom,
        &new_time,
        &new_lambda,
        &new_sigma,
        &n_new_time);
    if (err)
    {
        goto FREE_RETURN;
    }

    // Create the array of new_lambda2[i] = new_lambda[i] + gamma
    new_lambda2 = (double*)calloc(n_new_time, sizeof(double));

    if (!new_lambda2)
    {
        err = "Memory allocation failure in SrtGrfnLGM2FMClambda";
        goto FREE_RETURN;
    }

    for (i = 0; i < n_new_time; i++)
    {
        new_lambda2[i] = new_lambda[i] + gamma_dom;
    }

    // Checks that ALL new_lambda[i], new_lambda2[i], new_lambda[i] + new_lambda2[i] <>0
    for (i = 0; i < n_new_time; i++)
    {
        if (new_lambda[i] == 0.0)
        {
            err = "One of the lambda1 = 0.0";
            goto FREE_RETURN;
        }
        if (new_lambda2[i] == 0.0)
        {
            err = "One of the lambda2 = lambda1 + gamma = 0.0";
            goto FREE_RETURN;
        }
        if (new_lambda[i] + new_lambda2[i] == 0.0)
        {
            err = "One of the lambda1 + lambda2 = 2*lambda1 + gamma = 0.0";
            goto FREE_RETURN;
        }
    }

    /*	Get yield curve */
    dom_yc = (char*)get_ycname_from_irund(dom_und);

    /*	Get distributions */

    dom_fwd1     = (double*)calloc(nstp, sizeof(double));
    dom_fwd2     = (double*)calloc(nstp, sizeof(double));
    dom_exp1     = (double*)calloc(nstp, sizeof(double));
    dom_exp2     = (double*)calloc(nstp, sizeof(double));
    dom_phi1     = (double*)calloc(nstp, sizeof(double));
    dom_phi2     = (double*)calloc(nstp, sizeof(double));
    dom_phi12    = (double*)calloc(nstp, sizeof(double));
    dom_gam1_fwd = (double*)calloc(nstp, sizeof(double));
    dom_gam2_fwd = (double*)calloc(nstp, sizeof(double));
    dom_bond_pay = (double*)calloc(nstp, sizeof(double));
    dom_gam1_pay = (double*)calloc(nstp, sizeof(double));
    dom_gam2_pay = (double*)calloc(nstp, sizeof(double));

    covar = f3tensor(0, nstp - 1, 0, 1, 0, 1);

    if (!dom_fwd1 || !dom_fwd2 || !dom_phi1 || !dom_phi2 || !dom_phi12 || !dom_exp1 || !dom_exp2 ||
        !dom_gam1_fwd || !dom_gam2_fwd || !dom_bond_pay || !dom_gam1_pay || !dom_gam2_pay || !covar)
    {
        err = "Memory allocation error (3) in SrtGrfnLGM2FMClambda";
        goto FREE_RETURN;
    }

    // Calculates all parameters needed for Monte-Carlo
    fill_mc_init_lgm2f_lambda(
        do_jump,
        pay_date,
        pay_time,
        date,
        time,
        nstp,
        new_time,
        n_new_time,
        new_sigma,
        new_lambda,
        new_lambda2,
        alpha_dom,
        gamma_dom,
        rho_dom,
        dom_yc,
        dom_fwd1,
        dom_fwd2,
        dom_exp1,
        dom_exp2,
        dom_phi1,
        dom_phi2,
        dom_phi12,
        dom_gam1_fwd,
        dom_gam2_fwd,
        dom_bond_pay,
        dom_gam1_pay,
        dom_gam2_pay,
        covar);

    /*	Fill product structure */

    void_prm = (void**)calloc(nstp, sizeof(void*));

    if (!void_prm)
    {
        err = "Memory allocation error (4) in SrtGrfnLGM2FMClambda";
        goto FREE_RETURN;
    }

    dom_idx = 0;

    for (i = xStr.num_evt - 1; i >= 0; i--)
    {
        if (xStr.evt[i].evt)
        {
            grfn_prm         = malloc(sizeof(grfn_parm_mc2F));
            grfn_prm->global = &xStr;
            grfn_prm->local  = xStr.evt + i;

            grfn_prm->num_dom_df = xStr.evt[i].evt->dflen[dom_idx];
            grfn_prm->dom_df_tms = xStr.evt[i].evt->dft[dom_idx];
            grfn_prm->dom_df_dts = xStr.evt[i].evt->dfd[dom_idx];

            if (grfn_prm->num_dom_df > 0)
            {
                grfn_prm->dom_dff   = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam1  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam2  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam12 = dvector(0, grfn_prm->num_dom_df - 1);

                if (!grfn_prm->dom_dff || !grfn_prm->dom_gam1 || !grfn_prm->dom_gam2 ||
                    !grfn_prm->dom_gam12)
                {
                    err = "Memory allocation error (8) in SrtGrfnLGM2FMClambda";
                    goto FREE_RETURN;
                }

                for (k = 0; k < grfn_prm->num_dom_df; k++)
                {
                    // For this event date i, reconstruct all future zero coupons needed at time k
                    // with k>i
                    grfn_prm->dom_dff[k] =
                        swp_f_df(xStr.dts[i], grfn_prm->dom_df_dts[k], (char*)dom_yc);
                    err = gamma_lambda(
                        xStr.tms[i],
                        xStr.tms[i] + grfn_prm->dom_df_tms[k],
                        new_time,
                        new_lambda,
                        n_new_time,
                        &grfn_prm->dom_gam1[k]);
                    err = gamma_lambda(
                        xStr.tms[i],
                        xStr.tms[i] + grfn_prm->dom_df_tms[k],
                        new_time,
                        new_lambda2,
                        n_new_time,
                        &grfn_prm->dom_gam2[k]);
                    grfn_prm->dom_gam12[k] =
                        -0.5 *
                            (grfn_prm->dom_gam1[k] * grfn_prm->dom_gam1[k] * dom_phi1[i + flag] +
                             grfn_prm->dom_gam2[k] * grfn_prm->dom_gam2[k] * dom_phi2[i + flag]) -
                        grfn_prm->dom_gam1[k] * grfn_prm->dom_gam2[k] * dom_phi12[i + flag];
                }

                grfn_prm->do_dom = 1;
            }
            else
            {
                grfn_prm->do_dom = 0;
            }

            void_prm[i + flag] = (void*)grfn_prm;
        }
        else
        {
            void_prm[i + flag] = NULL;
        }
    }

    /*	Eventually! call to function */

    *prod_val = dmatrix(0, num_col - 1, 0, 1);

    if (!(*prod_val))
    {
        err = "Memory allocation error in SrtGrfnLGM2FMClambda";
        goto FREE_RETURN;
    }

    t2 = clock();

    smessage("Phase 1 -preprocessing, time in sec: %.2f", (double)(t2 - t1) / CLOCKS_PER_SEC);

    err = mc_main_lgm2f(
        /*	Time data */
        num_paths,
        num_col,
        time,
        date,
        nstp,
        do_jump,
        dom_fwd1,
        dom_fwd2,
        dom_exp1,
        dom_exp2,
        dom_phi1,
        dom_phi2,
        dom_phi12,
        dom_gam1_fwd,
        dom_gam2_fwd,
        dom_bond_pay,
        dom_gam1_pay,
        dom_gam2_pay,
        covar,
        void_prm,
        do_pecs,
        0,
        NULL,
        NULL,
        NULL,
        /*	Payoff function */
        grfn_payoff_lgm2f_mc, /*	Result */
        *prod_val);

    if (err)
        goto FREE_RETURN;

    df = swp_f_zr(today, pay_date, dom_yc);
    df = exp(-df * pay_time);

    *tableauCols = num_col;
    for (i = 0; i < num_col; i++)
    {
        (*prod_val)[i][0] *= df;
        (*prod_val)[i][1] *= df;
    }

    /*	Add PV of Past */
    (*prod_val)[num_col - 1][0] += xStr.gd->pv_of_past;

FREE_RETURN:

    if (time)
        free(time);
    if (date)
        free(date);
    if (dom_fwd1)
        free(dom_fwd1);
    if (dom_fwd2)
        free(dom_fwd2);
    if (dom_exp1)
        free(dom_exp1);
    if (dom_exp2)
        free(dom_exp2);
    if (dom_phi1)
        free(dom_phi1);
    if (dom_phi2)
        free(dom_phi2);
    if (dom_phi12)
        free(dom_phi12);
    if (dom_gam1_fwd)
        free(dom_gam1_fwd);
    if (dom_gam2_fwd)
        free(dom_gam2_fwd);
    if (dom_bond_pay)
        free(dom_bond_pay);
    if (dom_gam1_pay)
        free(dom_gam1_pay);
    if (dom_gam2_pay)
        free(dom_gam2_pay);

    // Free new structures
    if (new_time)
        free(new_time);
    if (new_lambda)
        free(new_lambda);
    if (new_lambda2)
        free(new_lambda2);
    if (new_sigma)
        free(new_sigma);

    if (covar)
        free_f3tensor(covar, 0, nstp - 1, 0, 1, 0, 1);

    if (sigma_date_dom)
        free(sigma_date_dom);
    if (sigma_dom)
        free(sigma_dom);

    if (void_prm)
    {
        for (i = 0; i < nstp; i++)
        {
            if (void_prm[i])
            {
                grfn_prm = (GRFNPARMMC2F)void_prm[i];

                if (grfn_prm->do_dom && grfn_prm->num_dom_df > 0)
                {
                    if (grfn_prm->dom_dff)
                        free_dvector(grfn_prm->dom_dff, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam1)
                        free_dvector(grfn_prm->dom_gam1, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam2)
                        free_dvector(grfn_prm->dom_gam2, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam12)
                        free_dvector(grfn_prm->dom_gam12, 0, grfn_prm->num_dom_df - 1);
                }

                free(grfn_prm);
            }
        }

        free(void_prm);
    }

    if (free_str)
    {
        FIRSTFreeMktStruct(&xStr);
    }

    return err;
}

// LGM2F MCEB with lambda term structure
char* SrtGrfnLGM2FMCEBlambda(
    char*      lgmund,
    int        numeventdates,
    long*      eventdates,
    int*       nUsedEventDates,
    int*       optimise,
    double*    fwd_iv,
    MCEBPARAMS params,
    long*      resRows,
    long       tableauRows,
    long*      tableauCols,
    char***    tableauStrings,
    int**      tableauMask,
    long       auxWidth,
    long*      auxLen,
    double**   aux,
    int        is_end_of_day_fixing,
    int        is_end_of_day_payment,
    long       num_paths,
    int        do_pecs,
    int        do_jump,
    double***  prod_val)
{
    int          free_str = 0;
    FIRSTAllMkts xStr;
    SrtGrfnParam defParm;
    int          forback;
    int          flag = 0;
    long         nstp;

    double *time = NULL, *date = NULL;

    double* sig_dom = NULL;

    double *dom_fwd1 = NULL, *dom_fwd2 = NULL, *dom_exp1 = NULL, *dom_exp2 = NULL, *dom_phi1 = NULL,
           *dom_phi2 = NULL, *dom_phi12 = NULL, *dom_gam1_fwd = NULL, *dom_gam2_fwd = NULL,
           *dom_bond_pay = NULL, *dom_gam1_pay = NULL, *dom_gam2_pay = NULL, ***covar = NULL;

    // New variables for lambda TS
    double *new_time = NULL, *new_lambda = NULL, *new_sigma = NULL, *new_lambda2 = NULL;

    long n_new_time;

    double *lam_dom = NULL, *lam_dom_time = NULL;

    int* ivOptimise = NULL;

    void**       void_prm = NULL;
    GRFNPARMMC2F grfn_prm;

    long today, spot_date;
    int  i, k, num_col, num_und;

    SrtUndPtr *und_ptr = NULL, dom_und;

    double alpha_dom, gamma_dom, rho_dom;

    char* dom_yc;
    char* domestic_name;

    double pay_time, df;

    double *sigma_date_dom = NULL, *sigma_dom = NULL;

    long pay_date;

    long n_sigma_dom, n_lam_dom;
    int  dom_idx;

    clock_t t1, t2;

    Err err = NULL;

    t1 = clock();

    /*	Initialise the GRFN tableau */

    /*	First, initialise the param struct */

    err                        = srt_f_set_default_GrfnParams(&defParm);
    defParm.num_MCarlo_paths   = num_paths;
    defParm.max_time_per_slice = 1000;
    defParm.min_nodes_per_path = 1;
    defParm.force_mc           = 1;
    defParm.jumping            = 1;

    /* End of Day Fixing */
    if (is_end_of_day_fixing)
    {
        defParm.end_of_day_fixings = SRT_YES;
        defParm.end_of_day_flg     = SRT_YES;
    }
    else
    {
        defParm.end_of_day_fixings = SRT_NO;
        defParm.end_of_day_flg     = SRT_NO;
    }

    /* End of Day Payment */
    if (is_end_of_day_payment)
    {
        defParm.end_of_day_payment = SRT_YES;
    }
    else
    {
        defParm.end_of_day_payment = SRT_NO;
    }

    err = FIRSTInitMktStruct(
        numeventdates,
        eventdates,
        tableauRows,
        *tableauCols,
        tableauStrings,
        tableauMask,
        auxWidth,
        auxLen,
        aux,
        lgmund,
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

    if (num_und != 1)
    {
        err = "Product should involve only one underlying";
        goto FREE_RETURN;
    }

    dom_und = und_ptr[0];

    /* look for the underlying name */
    dom_und = lookup_und(lgmund);
    if (!dom_und)
    {
        err = "cannot find the underlying";
        goto FREE_RETURN;
    }

    if (get_mdltype_from_irund(dom_und) == LGM)
    {
        domestic_name = dom_und->underl_name;
    }
    else
    {
        err = "Model must be LGM2F";
        goto FREE_RETURN;
    }

    if (strcmp(domestic_name, dom_und->underl_name))
    {
        err = "Tableau uses different underlying";
        goto FREE_RETURN;
    }

    /* look for the today date */
    today     = get_today_from_underlying(dom_und);
    spot_date = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

    num_col = xStr.num_cols;

    /*	Next, get the time steps */

    /*	Copy event dates */
    nstp = xStr.num_evt;
    while (nstp >= 1 && xStr.evt[nstp - 1].evt == NULL)
    {
        nstp--;
    }
    if (nstp < 1)
    {
        err = "No event in Tableau";
        goto FREE_RETURN;
    }

    time = (double*)calloc(nstp, sizeof(double));
    date = (double*)calloc(nstp, sizeof(double));

    if (!time || !date)
    {
        err = "Memory allocation error (1) in SrtGrfnLGM2FMCEBlambda";
        goto FREE_RETURN;
    }

    memcpy(time, xStr.tms, nstp * sizeof(double));

    for (i = 0; i < nstp; i++)
    {
        date[i] = today + DAYS_IN_YEAR * time[i];

        if (i > 0 && date[i] - date[i - 1] >= 1)
        {
            date[i] = (long)(date[i] + 1.0e-08);
            time[i] = YEARS_IN_DAY * (date[i] - today);
        }
    }

    if (time[0] > 0)
    {
        /* add the zero time */
        num_f_add_number(&nstp, &time, 0);
        num_f_sort_vector(nstp, time);
        nstp -= 1;
        num_f_add_number(&nstp, &date, today);
        num_f_sort_vector(nstp, date);
        flag = 1;
    }

    if (do_jump)
    {
        pay_date = (long)(date[1] + 1.0E-8);
        pay_time = time[1];
    }
    else
    {
        pay_date = (long)(date[nstp - 1] + 1.0E-8);
        pay_time = time[nstp - 1];
    }

    // Get the term structure with lambda term structure
    err = Get_LGM2F_TermStructure2(
        domestic_name,
        &sigma_dom,
        &sigma_date_dom,
        &n_sigma_dom,
        &lam_dom,
        &lam_dom_time,
        &n_lam_dom,
        &alpha_dom,
        &gamma_dom,
        &rho_dom);

    if (err)
    {
        goto FREE_RETURN;
    }

    // Merge the term structures of lambda and sigmas
    err = merge_lambda_sigma_ts(
        sigma_dom,
        sigma_date_dom,
        n_sigma_dom,
        lam_dom,
        lam_dom_time,
        n_lam_dom,
        &new_time,
        &new_lambda,
        &new_sigma,
        &n_new_time);
    if (err)
    {
        goto FREE_RETURN;
    }

    // Create the array of new_lambda2[i] = new_lambda[i] + gamma
    new_lambda2 = (double*)calloc(n_new_time, sizeof(double));

    if (!new_lambda2)
    {
        err = "Memory allocation failure in SrtGrfnLGM2FMCEBlambda";
        goto FREE_RETURN;
    }

    for (i = 0; i < n_new_time; i++)
    {
        new_lambda2[i] = new_lambda[i] + gamma_dom;
    }

    // Checks that ALL new_lambda[i], new_lambda2[i], new_lambda[i] + new_lambda2[i] <>0
    for (i = 0; i < n_new_time; i++)
    {
        if (new_lambda[i] == 0.0)
        {
            err = "One of the lambda1 = 0.0";
            goto FREE_RETURN;
        }
        if (new_lambda2[i] == 0.0)
        {
            err = "One of the lambda2 = lambda1 + gamma = 0.0";
            goto FREE_RETURN;
        }
        if (new_lambda[i] + new_lambda2[i] == 0.0)
        {
            err = "One of the lambda1 + lambda2 = 2*lambda1 + gamma = 0.0";
            goto FREE_RETURN;
        }
    }

    /*	Get Fx spot and yield curves */
    dom_yc = (char*)get_ycname_from_irund(dom_und);

    /*	Get distributions */

    dom_fwd1     = (double*)calloc(nstp, sizeof(double));
    dom_fwd2     = (double*)calloc(nstp, sizeof(double));
    dom_exp1     = (double*)calloc(nstp, sizeof(double));
    dom_exp2     = (double*)calloc(nstp, sizeof(double));
    dom_phi1     = (double*)calloc(nstp, sizeof(double));
    dom_phi2     = (double*)calloc(nstp, sizeof(double));
    dom_phi12    = (double*)calloc(nstp, sizeof(double));
    dom_gam1_fwd = (double*)calloc(nstp, sizeof(double));
    dom_gam2_fwd = (double*)calloc(nstp, sizeof(double));
    dom_bond_pay = (double*)calloc(nstp, sizeof(double));
    dom_gam1_pay = (double*)calloc(nstp, sizeof(double));
    dom_gam2_pay = (double*)calloc(nstp, sizeof(double));

    covar = f3tensor(0, nstp - 1, 0, 1, 0, 1);

    if (!dom_fwd1 || !dom_fwd2 || !dom_phi1 || !dom_phi2 || !dom_phi12 || !dom_exp1 || !dom_exp2 ||
        !dom_gam1_fwd || !dom_gam2_fwd || !dom_bond_pay || !dom_gam1_pay || !dom_gam2_pay || !covar)
    {
        err = "Memory allocation error (3) in SrtGrfnLGM2FMCEBlambda";
        goto FREE_RETURN;
    }

    // Calculates all parameters needed for Monte-Carlo
    fill_mc_init_lgm2f_lambda(
        do_jump,
        pay_date,
        pay_time,
        date,
        time,
        nstp,
        new_time,
        n_new_time,
        new_sigma,
        new_lambda,
        new_lambda2,
        alpha_dom,
        gamma_dom,
        rho_dom,
        dom_yc,
        dom_fwd1,
        dom_fwd2,
        dom_exp1,
        dom_exp2,
        dom_phi1,
        dom_phi2,
        dom_phi12,
        dom_gam1_fwd,
        dom_gam2_fwd,
        dom_bond_pay,
        dom_gam1_pay,
        dom_gam2_pay,
        covar);

    /*	Fill product structure */

    void_prm = (void**)calloc(nstp, sizeof(void*));

    if (!void_prm)
    {
        err = "Memory allocation error (4) in SrtGrfnLGM2FMCEBlambda";
        goto FREE_RETURN;
    }

    dom_idx = 0;

    for (i = xStr.num_evt - 1; i >= 0; i--)
    {
        if (xStr.evt[i].evt)
        {
            grfn_prm         = malloc(sizeof(grfn_parm_mc2F));
            grfn_prm->global = &xStr;
            grfn_prm->local  = xStr.evt + i;

            grfn_prm->num_dom_df = xStr.evt[i].evt->dflen[dom_idx];
            grfn_prm->dom_df_tms = xStr.evt[i].evt->dft[dom_idx];
            grfn_prm->dom_df_dts = xStr.evt[i].evt->dfd[dom_idx];

            if (grfn_prm->num_dom_df > 0)
            {
                grfn_prm->dom_dff   = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam1  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam2  = dvector(0, grfn_prm->num_dom_df - 1);
                grfn_prm->dom_gam12 = dvector(0, grfn_prm->num_dom_df - 1);

                if (!grfn_prm->dom_dff || !grfn_prm->dom_gam1 || !grfn_prm->dom_gam2 ||
                    !grfn_prm->dom_gam12)
                {
                    err = "Memory allocation error (8) in SrtGrfnLGM2FMCEBlambda";
                    goto FREE_RETURN;
                }

                for (k = 0; k < grfn_prm->num_dom_df; k++)
                {
                    // For this event date i, reconstruct all future zero coupons needed at time k
                    // with k>i
                    grfn_prm->dom_dff[k] =
                        swp_f_df(xStr.dts[i], grfn_prm->dom_df_dts[k], (char*)dom_yc);
                    err = gamma_lambda(
                        xStr.tms[i],
                        xStr.tms[i] + grfn_prm->dom_df_tms[k],
                        new_time,
                        new_lambda,
                        n_new_time,
                        &grfn_prm->dom_gam1[k]);
                    err = gamma_lambda(
                        xStr.tms[i],
                        xStr.tms[i] + grfn_prm->dom_df_tms[k],
                        new_time,
                        new_lambda2,
                        n_new_time,
                        &grfn_prm->dom_gam2[k]);
                    grfn_prm->dom_gam12[k] =
                        -0.5 *
                            (grfn_prm->dom_gam1[k] * grfn_prm->dom_gam1[k] * dom_phi1[i + flag] +
                             grfn_prm->dom_gam2[k] * grfn_prm->dom_gam2[k] * dom_phi2[i + flag]) -
                        grfn_prm->dom_gam1[k] * grfn_prm->dom_gam2[k] * dom_phi12[i + flag];
                }
                grfn_prm->do_dom = 1;
            }
            else
            {
                grfn_prm->do_dom = 0;
            }

            void_prm[i + flag] = (void*)grfn_prm;
        }
        else
        {
            void_prm[i + flag] = NULL;
        }
    }

    /*	Eventually! call to function */

    *tableauCols     = num_col;
    *resRows         = max(num_col + 1, nstp);
    *nUsedEventDates = xStr.num_evt;
    df               = swp_f_df(today, pay_date, dom_yc);

    /* create an optimisation vector from the input */
    ivOptimise = ivector(0, nstp - 1);
    if (flag)
        ivOptimise[0] = 0;
    for (i = flag; i < nstp; i++)
        ivOptimise[i] = optimise[i + numeventdates - xStr.num_evt - flag];

    if (params->iMultiIndex)
    {
        params->iNbIndex = params->iMultiIndex;
    }
    else
    {
        params->iNbIndex = 1;
    }

    mceb_allocate_params(params, nstp);

    if (params->iAdjustIV)
    {
        if (flag)
        {
            params->dMarketFwdIV[0] = 0.0;
        }

        for (i = flag; i < nstp; i++)
        {
            params->dMarketFwdIV[i] = fwd_iv[i + numeventdates - xStr.num_evt - flag] / df;
        }
    }

    *prod_val = dmatrix(0, *resRows - 1, 0, 2 + params->iNbIndex);

    if (!(*prod_val))
    {
        err = "Memory allocation error";
        goto FREE_RETURN;
    }

    t2 = clock();

    smessage("Phase 1 -preprocessing, time in sec: %.2f", (double)(t2 - t1) / CLOCKS_PER_SEC);

    err = mc_main_lgm2f(
        /*	Time data */
        num_paths,
        num_col,
        time,
        date,
        nstp,
        do_jump,
        dom_fwd1,
        dom_fwd2,
        dom_exp1,
        dom_exp2,
        dom_phi1,
        dom_phi2,
        dom_phi12,
        dom_gam1_fwd,
        dom_gam2_fwd,
        dom_bond_pay,
        dom_gam1_pay,
        dom_gam2_pay,
        covar,
        void_prm,
        do_pecs,
        1,
        ivOptimise,
        params,
        NULL,
        /*	Payoff function */
        grfn_payoff_lgm2f_mc, /*	Result */
        *prod_val);

    if (err)
        goto FREE_RETURN;

    /* Recopy Barrier / CoefLin for the moment */
    for (i = 0; i < nstp; i++)
    {
        (*prod_val)[i][2] = params->dBarrier[i];

        for (k = 0; k < params->iNbIndex; k++)
        {
            (*prod_val)[i][3 + k] = params->dCoefLin[i][k + 1];
        }
    }

    for (i = 0; i < num_col + 1; i++)
    {
        (*prod_val)[i][0] *= df;
        (*prod_val)[i][1] *= df;
    }

    for (i = 0; i < nstp; i++)
    {
        if (params->iCalcOneTime)
            params->dOneTimeCall[i] *= df;
        if (params->iCalcOneTimePartial)
            params->dOneTimePartial[i] *= df;
        if (params->iCalcIV || params->iAdjustIV)
            params->dModelFwdIV[i] *= df;
    }

    if (flag)
    {
        for (i = 0; i < nstp - 1; i++)
        {
            (*prod_val)[i][2] = (*prod_val)[i + 1][2];

            for (k = 0; k < params->iNbIndex; k++)
            {
                (*prod_val)[i][3 + k] = (*prod_val)[i + 1][3 + k];
            }
        }

        mceb_shift_extrainfos(params);
    }

    /*	Add PV of Past */
    (*prod_val)[num_col - 1][0] += xStr.gd->pv_of_past;
    (*prod_val)[num_col][0] += xStr.gd->pv_of_past;

FREE_RETURN:

    if (time)
        free(time);
    if (date)
        free(date);
    if (dom_fwd1)
        free(dom_fwd1);
    if (dom_fwd2)
        free(dom_fwd2);
    if (dom_exp1)
        free(dom_exp1);
    if (dom_exp2)
        free(dom_exp2);
    if (dom_phi1)
        free(dom_phi1);
    if (dom_phi2)
        free(dom_phi2);
    if (dom_phi12)
        free(dom_phi12);
    if (dom_gam1_fwd)
        free(dom_gam1_fwd);
    if (dom_gam2_fwd)
        free(dom_gam2_fwd);
    if (dom_bond_pay)
        free(dom_bond_pay);
    if (dom_gam1_pay)
        free(dom_gam1_pay);
    if (dom_gam2_pay)
        free(dom_gam2_pay);

    // Free new structures
    if (new_time)
        free(new_time);
    if (new_lambda)
        free(new_lambda);
    if (new_lambda2)
        free(new_lambda2);
    if (new_sigma)
        free(new_sigma);

    if (covar)
        free_f3tensor(covar, 0, nstp - 1, 0, 1, 0, 1);

    if (sigma_date_dom)
        free(sigma_date_dom);
    if (sigma_dom)
        free(sigma_dom);
    if (ivOptimise)
        free_ivector(ivOptimise, 0, nstp - 1);

    if (void_prm)
    {
        for (i = 0; i < nstp; i++)
        {
            if (void_prm[i])
            {
                grfn_prm = (GRFNPARMMC2F)void_prm[i];

                if (grfn_prm->do_dom && grfn_prm->num_dom_df > 0)
                {
                    if (grfn_prm->dom_dff)
                        free_dvector(grfn_prm->dom_dff, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam1)
                        free_dvector(grfn_prm->dom_gam1, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam2)
                        free_dvector(grfn_prm->dom_gam2, 0, grfn_prm->num_dom_df - 1);
                    if (grfn_prm->dom_gam12)
                        free_dvector(grfn_prm->dom_gam12, 0, grfn_prm->num_dom_df - 1);
                }

                free(grfn_prm);
            }
        }

        free(void_prm);
    }

    if (free_str)
    {
        FIRSTFreeMktStruct(&xStr);
    }

    return err;
}
