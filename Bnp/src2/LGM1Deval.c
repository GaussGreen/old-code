/****************************************************************************/
/* Purpose */
/****************************************************************************/

/****************************************************************************/
/* EXTERNAL ROUTINES USED */
/****************************************************************************/
/*

*/
/****************************************************************************/
/* Headers */
/****************************************************************************/
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_lgmUSprotos.h"
#include "srt_h_lgmeval.h"
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmtypes.h"

/****************************************************************************/
/* Definition */
/****************************************************************************/

#define SRT_MIN_REDUC_FACT 1.10

/****************************************************************************/
/* Prototypes of private functions */
/****************************************************************************/

/* static counter_getvolfudgetimeswap_new; */ /* for debug */

/* calculates the integer part and remainder of x */
static long IntPart(double x, double* remainder);

/* Interpolates value at k+dk from an array of values for convolver */
static double InterpfromArray(long k, double dk, long nk, double* TheArray);

/* same as before, but extrapolates flat */
static double InterpfromArrayTS(long k, double dk, long nk, double* TheArray);

/* Computes integration weights for convolver */
static void GenConvWeights(double* w, double h, long m);

/* Calculates the integration weights for the kink correction */
static void GetCorrection(double* w, double x);

/* Calculates the payoff at each exercise date for general European option */
static LGMErr GenEurpayoff(
    double**    payoff,
    long        nx,
    double*     x,
    double*     reduction,
    LGMDealType dealtype,
    void*       dealPtr,
    LGMCalParm* CalReq,
    Date        tNow,
    String      ycName,
    LGM_TS*     tsPtr,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double*));                    /* swaption/cap exponents (beta) */

/* Calculates the payoff at each exercise date for simple midatlantic deal */
static LGMErr MidAtpayoff(
    double**    payoff,
    long        nx,
    double*     x,
    double*     reduction,
    LGMDealType dealtype,
    void*       dealPtr,
    LGMCalParm* CalReq,
    Date        tNow,
    String      ycName,
    LGM_TS*     tsPtr,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double*));                    /* swaption/cap exponents (beta) */

/* Calculates the payoff at each exercise date for general midatlantic deal */
static LGMErr GenMidAtpayoff(
    double**    payoff,
    long        nx,
    double*     x,
    double*     reduction,
    LGMDealType dealtype,
    void*       dealPtr,
    LGMCalParm* CalReq,
    Date        tNow,
    String      ycName,
    LGM_TS*     tsPtr,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double*));                    /* swaption/cap exponents (beta) */

/* Calculates the payoff at each exercise date for Bermudan inverse floater */
static LGMErr BerInvFltpayoff(
    double**    payoff,
    long        nx,
    double*     x,
    double*     reduction,
    LGMDealType dealtype,
    void*       dealPtr,
    LGMCalParm* CalReq,
    Date        tNow,
    String      ycName,
    LGM_TS*     tsPtr,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double*));                    /* swaption/cap exponents (beta) */

/* Calculates the payoff at each exercise date for Bermudan Cap floater */
static LGMErr BerCapFltpayoff(
    double**    payoff,
    long        nx,
    double*     x,
    double*     reduc,
    LGMDealType dealtype,
    void*       dealPtr,
    LGMCalParm* CalReq,
    Date        tNow,
    String      ycName,
    LGM_TS*     tsPtr,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double*));                    /* swaption/cap exponents (beta) */

LGMErr GetVol_SABRValues(
    long               tNow,
    long               value_date,
    long               start,
    long               end,
    SrtCallTimeSwapPtr deal,
    double*            alpha,
    double*            beta,
    double*            sigma,
    double*            rho);

LGMErr GetVol_FudgeTimeSwap_new(
    long               tNow,
    long               value_date,
    long               start,
    long               end,
    double             strike,
    double             fra,
    SrtCallTimeSwapPtr deal,
    double*            result,
    double             alpha,
    double             beta,
    double             sigma,
    double             rho);

/****************************************************************************/
/**** Comments ****/
/* Routines to be completed later		*/
/* LGMErr ConvEvalSimMidAtwBarrier()	*/

/**** HOLIDAYS ***/
/* Wherever we need to use the add_unit function,
the currency code ccy is available */

LGMErr NewMidAtEval(
    ConvParams* parms,    /* convolution numerical constants */
                          /* info about today */
    Date        tNow,     /* calculation date */
    int         endofday, /* 1 = cannot exercise today, 0 = can exercise today */
    String      ycName,
    LGMCalParm* CalReq,
    LGM_TS*     tsPtr,                                          /* calibrated term structure */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double*),                     /* swaption/cap exponents (beta) */
                                                                /* information about the deal */
    LGMDealType dealtype,                                       /* type of deal */
    void*       dealPtr,                                        /* definition of the deal */
                                                                /* output */
    double*  LGMvalue,                                          /* value of the deal */
    long*    nExBdry, /* number of points in the exercise boundary */
    Date**   tExBdry, /* number of dates in the exercise boundary */
    double** xExBdry) /* array of exercise points (in x) */
{
    LGMErr             error;
    LGMMarkConv        conventions;
    double *           zeta, gamma;
    double             LGMval, *LGMvalVega;
    long               nEx, j, k;
    SrtSimMidAtPtr     MidAt;
    SrtGenMidAtPtr     GenExOnce;
    SrtCallInvFltPtr   BerInvFlt;
    SrtCallCapFltPtr   BerCapFlt;
    SrtCallTimeSwapPtr BerTimeSwap;
    Date*              tExArr;
    Date               tfirst;
    SrtCrvPtr          yldcrv;
    double *           Zeta1, *Zeta2, *Zeta12, *H1, *H2;
    double **          Zeta1Scenari, **Zeta2Scenari, **Zeta12Scenari, **H1Scenari, **H2Scenari;
    long               lStartTime, lEndTime;

    yldcrv = lookup_curve(ycName);

    error = NULL;
    if (dealtype >= LGMLastDealType)
        return ("unknown deal type");

    error = LGMCcyDefaults(yldcrv, &conventions); /* get currency code */
    if (error)
        return (error);

    /* Get number of exercise and zeta values at each valid exercise date */
    tfirst = tNow + 1 - 1;
    if (endofday == 1)
        tfirst = add_unit(
            tNow,
            1,
            SRT_BDAY,
            SUCCEEDING /*, conv.ccy*/); /* tfirst is the first possible exercise */

    switch (dealtype)
    {
        /* Simple Midatlantic option */
    case SimpleMidAt:
        MidAt = (SrtSimMidAt*)dealPtr; /* re-cast as correct deal type */

        /* find exercise dates on or after tfirst*/
        nEx = MidAt->nEx;
        for (j = 0; j < nEx && (MidAt->tEx[j] < tfirst); j++)
            ;
        nEx = nEx - j;
        if (nEx < 1)
        {
            *LGMvalue = 0;
            return (error);
        }
        MidAt->FirstExer = j; /* ensure evaluation begins at right exercise */

        /* allocate zeta array */

        if (CalReq->LGMOneTwoFactor == 2)
        {
            tExArr = (Date*)srt_calloc(nEx, sizeof(Date));
            if (tExArr == NULL)
                return ("allocation failed in MidAtEval");

            for (j = 0; j < nEx; j++)
                tExArr[j] = MidAt->tEx[j];

            Zeta1  = tsPtr->Zeta1;
            Zeta2  = tsPtr->Zeta2;
            Zeta12 = tsPtr->Zeta12;
            H1     = tsPtr->H1;
            H2     = tsPtr->H2;
            gamma  = tsPtr->gamma;

            Zeta1Scenari  = tsPtr->Zeta1Scenari;
            Zeta2Scenari  = tsPtr->Zeta2Scenari;
            Zeta12Scenari = tsPtr->Zeta12Scenari;
            H1Scenari     = tsPtr->H1Scenari;
            H2Scenari     = tsPtr->H2Scenari;

            lStartTime = time(NULL);

            error = LGMAutoCal2DTree(
                nEx,
                tExArr,
                Zeta1,
                Zeta12,
                Zeta2,
                H1,
                H2,
                gamma,
                tNow,
                ycName,
                GetVol,
                GetBeta,
                (void*)dealPtr,
                Midat2DPayoff,
                &LGMval,
                xExBdry);

            lEndTime = time(NULL);
            smessage("Midatlantic Pricing achieved in %.2f s", difftime(lEndTime, lStartTime));

            LGMvalVega = (double*)srt_calloc(2 * tsPtr->numZ, sizeof(double));

            /*
            error = LGMAutoCal2DVega(2*tsPtr->numZ, nEx, tExArr, Zeta1Scenari, Zeta12Scenari,
                                                            Zeta2Scenari, H1Scenari, H2Scenari,
            gamma, tNow, ycName, GetVol, GetBeta, (void*) dealPtr, Midat2DPayoff, LGMvalVega,
            xExBdry);
            */
            srt_free(LGMvalVega);

            *LGMvalue = LGMval;
        }
        else if (CalReq->LGMOneTwoFactor == 1)
        {
            zeta   = (double*)srt_calloc(nEx, sizeof(double));
            tExArr = (Date*)srt_calloc(nEx, sizeof(Date));
            if (zeta == NULL || tExArr == NULL)
                return ("allocation failed in LGMEval");

            /* find today's zetas at the exercise dates */
            for (k = 0; k < nEx; k++)
            {
                tExArr[k] = MidAt->tEx[MidAt->FirstExer + k];
                zeta[k]   = LGMZetaFromTS(tExArr[k], tNow, tsPtr);
            }

            error = Convolver(
                nEx,
                zeta,
                parms,
                tNow,
                ycName,
                GetVol,
                GetBeta,
                tsPtr,
                dealtype,
                (void*)dealPtr,
                MidAtpayoff,
                CalReq,
                &LGMval,
                xExBdry);

            *LGMvalue = LGMval;
        }
        break;

        /* General Exercise once (out of many) option */
    case GenMidAt:
        GenExOnce = (SrtGenMidAt*)dealPtr; /* re-cast as correct deal type */

        /* find exercise dates on or after tfirst*/
        nEx = GenExOnce->nEx;
        for (j = 0; j < nEx && (GenExOnce->tEx[j] < tfirst); j++)
            ;
        nEx = nEx - j;
        if (nEx < 1)
        {
            *LGMvalue = 0;
            return (error);
        }
        GenExOnce->FirstExer = j; /* Ensure that evaluation begins at right exercise */

        /* allocate zeta array */
        zeta   = (double*)srt_calloc(nEx, sizeof(double));
        tExArr = (Date*)srt_calloc(nEx, sizeof(Date));
        if (zeta == NULL || tExArr == NULL)
            return ("allocation failed in MidAtEval");

        /* find today's zetas at the exercise dates */
        for (k = 0; k < nEx; k++)
        {
            tExArr[k] = GenExOnce->tEx[j + k];
            zeta[k]   = LGMZetaFromTS(tExArr[k], tNow, tsPtr);
        }

        /* evaluate deal */
        error = Convolver(
            nEx,
            zeta,
            parms,
            tNow,
            ycName,
            GetVol,
            GetBeta,
            tsPtr,
            dealtype,
            (void*)dealPtr,
            GenMidAtpayoff,
            CalReq,
            &LGMval,
            xExBdry);
        *LGMvalue = LGMval;
        break;

        /* Bermudan inverse floater */
    case CallInvFlt:
        BerInvFlt = (SrtCallInvFlt*)dealPtr; /* re-cast as correct deal type */

        /* find exercise dates on or after tfirst*/
        nEx = BerInvFlt->nEx;
        for (j = 0; j < nEx && (BerInvFlt->tEx[j] < tfirst); j++)
            ;
        nEx = nEx - j;
        if (nEx < 1)
        {
            *LGMvalue = 0;
            return (error);
        }
        BerInvFlt->FirstEx = j; /* ensure eval begins at right exercise */

        /* allocate zeta array */
        zeta   = (double*)srt_calloc(nEx, sizeof(double));
        tExArr = (Date*)srt_calloc(nEx, sizeof(Date));
        if (zeta == NULL || tExArr == NULL)
            return ("allocation failed in LGMEval");

        /* find today's zetas at the exercise dates */
        for (k = 0; k < nEx; k++)
        {
            tExArr[k] = BerInvFlt->tEx[j + k];
            zeta[k]   = LGMZetaFromTS(tExArr[k], tNow, tsPtr);
        }

        /* evaluate deal */
        error = Convolver(
            nEx,
            zeta,
            parms,
            tNow,
            ycName,
            GetVol,
            GetBeta,
            tsPtr,
            dealtype,
            (void*)dealPtr,
            BerInvFltpayoff,
            CalReq,
            &LGMval,
            xExBdry);
        *LGMvalue = LGMval;
        break;

        /* Bermudan Cap floater */
    case CallCapFlt:
        BerCapFlt = (SrtCallCapFlt*)dealPtr; /* re-cast as correct deal type */

        /* find exercise dates on or after tfirst*/
        nEx = BerCapFlt->nEx;
        for (j = 0; j < nEx && (BerCapFlt->tEx[j] < tfirst); j++)
            ;
        nEx = nEx - j;
        if (nEx < 1)
        {
            *LGMvalue = 0;
            return (error);
        }
        BerCapFlt->FirstEx = j; /* ensure eval begins at right exercise */

        /* allocate zeta array */
        zeta   = (double*)srt_calloc(nEx, sizeof(double));
        tExArr = (Date*)srt_calloc(nEx, sizeof(Date));
        if (zeta == NULL || tExArr == NULL)
            return ("allocation failed in LGMEval");

        /* find today's zetas at the exercise dates */
        for (k = 0; k < nEx; k++)
        {
            tExArr[k] = BerCapFlt->tEx[j + k];
            zeta[k]   = LGMZetaFromTS(tExArr[k], tNow, tsPtr);
        }

        /* evaluate deal */
        error = Convolver(
            nEx,
            zeta,
            parms,
            tNow,
            ycName,
            GetVol,
            GetBeta,
            tsPtr,
            dealtype,
            (void*)dealPtr,
            BerCapFltpayoff,
            CalReq,
            &LGMval,
            xExBdry);
        *LGMvalue = LGMval;
        break;

        /* Callable Time Swap */
    case CallTimeSwap:
        BerTimeSwap = (SrtCallTimeSwap*)dealPtr; /* re-cast as correct deal type */

        /* find exercise dates on or after tfirst*/
        nEx = BerTimeSwap->nEx;
        for (j = 0; j < nEx && (BerTimeSwap->tEx[j] < tfirst); j++)
            ;
        nEx = nEx - j;
        if (nEx < 1)
        {
            *LGMvalue = 0;
            return (error);
        }
        parms->nx            = BerTimeSwap->nx;
        BerTimeSwap->FirstEx = j; /* ensure eval begins at right exercise */

        /* allocate zeta array */
        zeta   = (double*)srt_calloc(nEx, sizeof(double));
        tExArr = (Date*)srt_calloc(nEx, sizeof(Date));
        if (zeta == NULL || tExArr == NULL)
            return ("allocation failed in LGMEval");

        /* find today's zetas at the exercise dates */
        for (k = 0; k < nEx; k++)
        {
            tExArr[k] = BerTimeSwap->tEx[j + k];
            zeta[k]   = LGMZetaFromTS(tExArr[k], tNow, tsPtr);
        }

        /* evaluate deal */
        /*		error = ConvolverTimeSwap(nEx, zeta, parms, tNow, ycName, GetVol, GetBeta,
                                                        tsPtr, dealtype, (void*)dealPtr,
           BerTimeSwapPayoff, CalReq, &LGMval, xExBdry);	*/
        error = ComputeBermudanAndEuropeanTimeSwap(
            nEx,
            zeta,
            parms,
            tNow,
            ycName,
            GetVol,
            GetBeta,
            tsPtr,
            dealtype,
            (void*)dealPtr,
            BerTimeSwapPayoff,
            CalReq,
            &LGMval,
            xExBdry);
        *LGMvalue = LGMval;
        break;

    }               /* end switch(dealtype)
              
/* free and return */
    *nExBdry = nEx; /* number of points in the exercise boundary */
    if (*tExBdry != NULL)
        srt_free(*tExBdry);
    *tExBdry = tExArr; /* dates in the exercise boundary */
                       /* xExBdry set in Convolver */
    if (CalReq->LGMOneTwoFactor == 1)
    {
        srt_free(zeta);
    }
    return (error);
}

/**************************************************************************************/
/*********************** Convolver ****************************************************/
/* LGM replacement for tree */
/* Uses convolution to evaluate deals within a LGM framework **************************/
/* Outputs:
        *answer = value of the deal,
        xExBdry[0,1,...,nEx-1] are the x values of the exercise boundary at each exercise
                                                                if multiple exercise points, gives
   the nearest to the money
*/
/* Convolver requires three sets of arguments
        1) the number of exercises and the zeta's at the exercise points
        2a) information about the market (including today, a yield curve, and a calibrated
                        LGM term structure
        2b) information about caplet/swaption vols (for deals whose payoffs include options)
        3) information about the deal, and a function to calculate the payoff for all
                        the exercise dates
The arguments of the payoff function are:
        (*payofffunc)(double payoff[][], long nx, double x[], double reduction[]
                                        LGMDealType dealtype, void *dealPtr,
                                        Date EvalDate, String ycName, LGM_TS tsPtr)
 */

LGMErr Convolver(
    /* info about convolutions */
    long        nEx,   /* number of exercises */
    double*     zeta,  /* [0, 1, ..., nEx-1], values of zeta at the exercise dates */
    ConvParams* parms, /* convolution numerical constants */

    /* info about today's discount curve and swaption/cap vols */
    Date   EvalDate,                                            /* tNow for deal evaluation */
    String ycName,                                              /* yield curve name */
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double*),                     /* swaption/cap exponents (beta) */
    LGM_TS* tsPtr,                                              /* calibrated LGM term structure */

    /* information about the deal */
    LGMDealType dealtype,
    void*       dealPtr,
    LGMErr (*payofffunc)(),

    LGMCalParm* CalReq,

    /* output */
    double*  answer,  /* value of the deal */
    double** xExBdry) /* array of exercise points (in x) */
{
    /* declarations */
    double *varr, *x, *Qarr, *Qplus, *differ, **payofftable, *payoff;
    double *vmax, *reduction;
    double* x_exerarr;
    double* weights;
    double  Cweight[4];
    double  gridwidth, maxh, stencil, h, dx;
    double  var, totalvar, totalstddev, previouszeta;
    double  alpha, gamma, Qinterp, Pinterp;
    double  a1, a2, a3, exactkink, xkink, y, theta, kshift;
    double  corr, ratio;
    long    nx, n0, nz, m, i, j, k;
    long    killkinks, integershift, knew, ipos;
    long    kk, ifoundkink, icor;
    LGMErr  error;
    int     integrateagain;

    integrateagain = 0;
startagain:
    /* eliminate trivial case */
    if (nEx < 1)
    {
        *answer = 0.0;
        return (NULL);
    }
    /* unload the convolution parameters and set them to their default values if needed */
    gridwidth = parms->gridwidth;
    if (gridwidth < 2.0 || gridwidth > 8.0)
        gridwidth = 6.0;

    nx = parms->nx;
    if (nx < 32 || nx > 576)
        nx = 192;
    n0 = ((long)(nx / 2)); /* make sure nx is even */
    nx = 2 * n0;

    maxh = parms->h;
    if (maxh < 0.02 || maxh > 0.5)
        maxh = 0.0625;

    stencil = parms->stencil;
    if (stencil < 2.0 || stencil > 8.0)
        stencil = 6.0;

    killkinks = parms->killkinks;
    if (killkinks != 0)
        killkinks = 1;

    /* Calculate the v[j] array (standard deviations) */
    varr      = NULL;
    vmax      = NULL;
    reduction = NULL;
    varr      = (double*)srt_calloc(nEx, sizeof(double));
    vmax      = (double*)srt_calloc(nEx, sizeof(double));
    reduction = (double*)srt_calloc(nEx, sizeof(double));

    if (varr == NULL || vmax == NULL || reduction == NULL)
    {
        srt_free(varr);
        srt_free(vmax);
        srt_free(reduction);
        return ("first allocation failed in convolver");
    }

    totalvar = zeta[nEx - 1];
    if (totalvar < 1.e-12)
    {
        totalvar = 1.e-12; /* Clearly an intrinsic value calculation the hard way */
        for (j = 0; j < nEx; j++)
            zeta[j] = totalvar * ((double)(j + 1)) / ((double)nEx);
    }
    totalstddev = sqrt(totalvar);

    /* find varr[j], the std dev from tEx[j-1] to tEx[j] */
    /* if the change in variance is non-positive, signal with varr[j] = -1 */
    previouszeta = 0.;
    for (j = 0; j < nEx; j++)
    {
        var = zeta[j] - previouszeta;
        if (var < 1.e-7 * totalvar)
        {
            varr[j] = -1.0;
            zeta[j] = previouszeta;
        }
        else
            varr[j] = sqrt(var);
        previouszeta = zeta[j];
    }

    vmax[0] = max(0, varr[0]); /* vmax[j] is max of varr[k] for k<=j */
    for (j = 1; j < nEx; j++)
        vmax[j] = max(vmax[j - 1], varr[j]);

    /* STEP 1: Set up x and z grids, and allocate space */
    /* x[k] = (k-n0)*dx*reduction[j], k=0,1,..., nx, j=0,1,2,...,nEx-1 */
    dx = 2.0 * totalstddev * gridwidth / ((double)nx);

    /* z[i] = (i-m)*h, i = 0, 1, ..., 2m=nz */
    /* integration intervals from z[i] to z[i+1], for i = 0, 1, ..., 2m-1 = nz-1 */
    h = 0.9233333 * dx / vmax[nEx - 1]; /* 0.9233 to prevent staircase from roundoff */
    h = min(h, maxh);

    m  = ((long)(1.0 + stencil / h));
    nz = 2 * m;

    /* find out where grid spacing can be reduced */
    /* can be reduced if reduced grid
            a. keeps totalstddev standard deviations (from zero to tEx[j]) on grid
            b. keeps delta_z less than delta_x
    */
    reduction[nEx - 1] = 1.0;
    for (j = nEx - 2; j >= 0; j--)
    {
        reduction[j] = reduction[j + 1];
        if (varr[j] > 0)
        {
            while (sqrt(zeta[j]) < 0.5 * totalstddev * reduction[j] &&
                   h < 0.5 * dx * reduction[j] / vmax[j])
                reduction[j] = 0.5 * reduction[j];
        }
    }

    /* Allocate arrays */

    x           = NULL;
    Qarr        = NULL;
    Qplus       = NULL;
    payofftable = NULL;
    weights     = NULL;
    if ((xExBdry) && (*xExBdry != NULL))
        srt_free(*xExBdry);

    x           = (double*)srt_calloc(nx + 1, sizeof(double));
    Qarr        = (double*)srt_calloc(nx + 1, sizeof(double));
    Qplus       = (double*)srt_calloc(nx + 1, sizeof(double));
    differ      = (double*)srt_calloc(nx + 1, sizeof(double));
    payofftable = dmatrix(0, nEx - 1, 0, nx);
    weights     = (double*)srt_calloc(nz + 1, sizeof(double));
    x_exerarr   = (double*)srt_calloc(nEx, sizeof(double));

    if (x == NULL || Qarr == NULL || Qplus == NULL || differ == NULL || payofftable == NULL ||
        weights == NULL || x_exerarr == NULL)
    {
        srt_free(varr);
        srt_free(vmax);
        srt_free(reduction);
        srt_free(x);
        srt_free(Qarr);
        srt_free(Qplus);
        srt_free(differ);
        srt_free(weights);
        srt_free(x_exerarr);
        free_dmatrix(payofftable, 0, nEx - 1, 0, nx);
        return ("allocation failed in convolver");
    }

    if (xExBdry)
        *xExBdry = x_exerarr;

    /* STEP 2: Generate weights and payoff table */
    GenConvWeights(weights, h, m);
    for (k = 0; k <= nx; k++)
        x[k] = dx * ((double)(k - n0));

    /* call payoff function to set up calculation */
    error = (*payofffunc)(
        payofftable,
        nx,
        x,
        reduction, /* output and grid */
        dealtype,
        dealPtr,
        CalReq, /* deal */
        EvalDate,
        ycName,
        tsPtr,
        GetVol,
        GetBeta); /* today's market & tsPtr */

    /* Backward induction */
    for (k = 0; k <= nx; k++) /* Value is zero after last exercise has passed */
        Qarr[k] = 0;

    for (j = (nEx - 1); j >= 0; j--) /* For each exercise date (starting from the last) */
    {
        if (j < nEx - 1 && reduction[j] < 0.90 * reduction[j + 1]) /* reduce grid size if needed */
        {
            ratio = reduction[j] / reduction[j + 1];
            for (knew = 0; knew <= nx; knew++) /* transfer results */
                Qplus[knew] = Qarr[knew];
            for (knew = 0; knew <= nx; knew++) /* interpolate to get Qarr on new grid */
            {
                k = IntPart(ratio * ((double)(knew - n0)) + (double)n0, &gamma);
                if (dealtype == CallTimeSwap)
                {
                    Qarr[knew] = InterpfromArrayTS(k, gamma, nx, Qplus);
                }
                else
                {
                    Qarr[knew] = InterpfromArray(k, gamma, nx, Qplus);
                }
            }
        }
        payoff = payofftable[j];
        for (k = 0; k <= nx; k++)
        {
            Qplus[k]  = 0;
            differ[k] = Qarr[k] - payoff[k];
        }
        if (varr[j] > 0) /* integrate */
        {
            alpha = h * varr[j] /
                    (dx * reduction[j]); /* Qplus[k] = sum of weight[i]*max{Qarr[k*], payoff[k*]} */
            for (i = 0; i <= nz; i++)    /*  at k* = k + (i-m)*alpha */
            {
                kshift       = alpha * ((double)(i - m)); /* set up interpolation */
                integershift = IntPart(kshift, &gamma);
                for (k = 0; k <= nx; k++)
                {
                    knew = k + integershift; /* interpolate to get value at k+integershift+gamma */
                    Qinterp  = InterpfromArray(knew, gamma, nx, Qarr);
                    Pinterp  = InterpfromArray(knew, gamma, nx, payoff);
                    Qplus[k] = Qplus[k] + weights[i] * max(Qinterp, Pinterp);
                }
            }
        }
        else
        {
            for (k = 0; k <= nx; k++)
                Qplus[k] = max(Qarr[k], payoff[k]);
        }
        /* Correct values for any kinks */
        /* should use continue statements to straighten out this nest of vipers */
        ifoundkink = 0; /* Search to find each crossing (kink) */
        for (kk = 1; kk <= nx - 2; kk++)
        {
            if ((differ[kk] > 0 && differ[kk + 1] > 0) || (differ[kk] < 0 && differ[kk + 1] < 0) ||
                (differ[kk - 1] > 0 && differ[kk] == 0 && differ[kk + 1] > 0) ||
                (differ[kk - 1] < 0 && differ[kk] == 0 && differ[kk + 1] < 0))
                continue; /* No kink go to next point */

            a1 = differ[kk];                  /* otherwise there's a kink between kk, kk+1 */
            a2 = differ[kk + 1] - differ[kk]; /* locate it exactly */
            a3 = 0.25 * (differ[kk + 2] + differ[kk - 1] - differ[kk] - differ[kk + 1]);
            a2 = a2 - a3;

            if (fabs(a3 / a2) < 1.e-06)
                exactkink = -(a1 / a2) * (1.0 + a1 * a3 / (a2 * a2));
            else if (a2 > 0)
                exactkink = (-a2 + sqrt(a2 * a2 - 4.0 * a1 * a3)) / (2.0 * a3);
            else
                exactkink = (-a2 - sqrt(a2 * a2 - 4.0 * a1 * a3)) / (2.0 * a3);
            exactkink = exactkink + (double)kk; /* kink is at exactkink */

            xkink = (exactkink - (double)n0) * dx * reduction[j]; /* store kink for output */
            if (ifoundkink == 0 || fabs(xkink) < fabs(x_exerarr[j]))
                x_exerarr[j] = xkink;
            ifoundkink = 1;

            if (killkinks == 0 || varr[j] <= 0 || alpha <= 1.e-4)
                continue; /* continue to the next point unless we need to kill the kink */

            for (k = 0; k <= nx; k++) /* Correct Qplus[k] for kink */
            {
                y = (exactkink - (double)k) / alpha + (double)m;
                i = IntPart(y, &theta);

                if (i >= 2 && i <= nz - 3)
                {
                    GetCorrection(Cweight, theta);
                    corr   = 0.0;
                    kshift = ((double)k) + alpha * ((double)(i - 1 - m));
                    for (icor = 0; icor <= 3; icor++)
                    {
                        knew    = IntPart(kshift, &gamma);
                        Qinterp = fabs(InterpfromArray(knew, gamma, nx, differ));
                        ipos    = i - 1 + icor;
                        corr    = corr + Cweight[icor] * weights[ipos] * Qinterp;
                        kshift  = kshift + alpha;
                    }
                    Qplus[k] = Qplus[k] + corr;
                }
            } /* end k loop */
        }     /* end kk loop */

        if (ifoundkink == 0) /* if no kink was found */
        {
            x_exerarr[j] = n0 * dx; /* put exer boundary at edge of grid */
            if (fabs(differ[0]) < fabs(differ[nx]))
                x_exerarr[j] = -x_exerarr[j];
        }

        for (k = 0; k <= nx; k++) /* update Qarray */
            Qarr[k] = Qplus[k];
    } /* end loop over exercise dates */
    *answer = Qarr[n0];

    /* free up space allocated */
    srt_free(varr);
    srt_free(vmax);
    srt_free(reduction);
    srt_free(differ);
    srt_free(x);
    srt_free(Qarr);
    srt_free(Qplus);
    free_dmatrix(payofftable, 0, nEx - 1, 0, nx);
    srt_free(weights);
    if (!xExBdry)
        srt_free(x_exerarr);

    /* debug */
    if (integrateagain == 1)
        goto startagain;

    return (NULL);
}

/*****************************************************************************************/
/* Integer truncator for convolver */
static long IntPart(double x, double* remainder)
{
    long intpart;
    intpart    = (long)x;
    *remainder = x - (double)intpart;
    if (*remainder < 0)
    {
        *remainder = *remainder + 1.0;
        intpart    = intpart - 1;
    }
    return (intpart);
}

/*****************************************************************************************/
/* Interpolator for convolver
Function: Interpolates value at k+dk from an array of values
TheArray[k], k=0, ..., nk. */
static double InterpfromArray(long k, double dk, long nk, double* TheArray)
{
    double interpval, x;

    if (k > 0 && k < (nk - 1))
        interpval = TheArray[k] + dk * (TheArray[k + 1] - TheArray[k]) -
                    0.25 * dk * (1.0 - dk) *
                        (TheArray[k - 1] - TheArray[k] - TheArray[k + 1] + TheArray[k + 2]);

    else if (k <= 0)
    {
        x         = dk + (double)k;
        interpval = TheArray[0] + x * (TheArray[1] - TheArray[0]) -
                    0.5 * x * (1.0 - x) * (TheArray[0] - 2.0 * TheArray[1] + TheArray[2]);
    }

    else
    {
        x = dk + (double)(k + 1 - nk);
        interpval =
            TheArray[nk - 1] + x * (TheArray[nk] - TheArray[nk - 1]) -
            0.5 * x * (1.0 - x) * (TheArray[nk] - 2.0 * TheArray[nk - 1] + TheArray[nk - 2]);
    }

    return (interpval);
}

/*****************************************************************************************/
/* Interpolator for convolver
Function: Interpolates value at k+dk from an array of values
TheArray[k], k=0, ..., nk. For time swap, extrapolates flat; this solves a problem for very
small maturities  */
static double InterpfromArrayTS(long k, double dk, long nk, double* TheArray)
{
    double interpval, x;

    if (k > 0 && k < (nk - 1))
        interpval = TheArray[k] + dk * (TheArray[k + 1] - TheArray[k]) -
                    0.25 * dk * (1.0 - dk) *
                        (TheArray[k - 1] - TheArray[k] - TheArray[k + 1] + TheArray[k + 2]);

    else if (k <= 0)
    {
        x         = dk + (double)k;
        interpval = TheArray[0] + x * (TheArray[1] - TheArray[0]) -
                    0.5 * x * (1.0 - x) * (TheArray[0] - 2.0 * TheArray[1] + TheArray[2]);
        /* TRY FOR EXTRAPOLATION FLAT */
        interpval = TheArray[0];
    }

    else
    {
        x = dk + (double)(k + 1 - nk);
        interpval =
            TheArray[nk - 1] + x * (TheArray[nk] - TheArray[nk - 1]) -
            0.5 * x * (1.0 - x) * (TheArray[nk] - 2.0 * TheArray[nk - 1] + TheArray[nk - 2]);
        /* TRY FOR EXTRAPOLATION FLAT */
        interpval = TheArray[nk];
    }

    return (interpval);
}

/**************************************************************************************/
/* Computes integration weights for convolver */
static void GenConvWeights(double* w, double h, long m)
{
    long   i, nz;
    double z;
    double sum, sum2, sum4, a, b;

    nz = 2 * m;
    for (i = 0; i <= m; i++)
    {
        z    = h * ((double)(i - m));
        w[i] = h * LGMsafeGauss(z);
    }

    for (i = 0; i < m; i++)
        w[nz - i] = w[i];

    /* End effects */
    z         = -h * ((double)m);
    w[0]      = w[0] * 11.0 / 24.0 + norm(z);
    w[2]      = w[2] * 25.0 / 24.0;
    w[nz]     = w[0];
    w[nz - 2] = w[2];

    /* Renormalize */
    sum  = 0.0;
    sum2 = 0.0;
    sum4 = 0.0;
    for (i = 0; i <= nz; i++)
    {
        z    = h * ((double)(i - m));
        sum  = sum + w[i];
        sum2 = sum2 + w[i] * z * z;
        sum4 = sum4 + w[i] * z * z * z * z;
    }

    a = (sum4 - sum2) / (sum4 * sum - sum2 * sum2);
    b = (sum - sum2) / (sum4 * sum - sum2 * sum2);

    for (i = 0; i <= nz; i++)
    {
        z    = h * ((double)(i - m));
        w[i] = w[i] * (a + b * z * z);
    }
    return;
}

/********************************************************************/
/* routine to calculate the correction weights for the kink removal */
static void GetCorrection(double* w, double x)
{
    double xx, yy;
    if (x > 1.0)
        x = 1.0;
    if (x < 0.0)
        x = 0.0;
    xx   = x * x;
    yy   = (1.0 - x) * (1.0 - x);
    w[0] = yy * (2.0 - yy) / 24.;
    w[1] = (1.0 - yy * (13.0 + 2.0 * x - 3.0 * xx)) / 24.0;
    w[2] = (1.0 - xx * (12.0 + 4.0 * x - 3.0 * xx)) / 24.0;
    w[3] = xx * (2.0 - xx) / 24.0;

    return;
}

/********************************************************************************/
/********************************************************************************/
/* Payoff routines: Calculates the payoff at each exercise date */
static LGMErr GenEurpayoff(
    double**    payoff,
    long        nx,
    double*     x,
    double*     reduction,
    LGMDealType dealtype,
    void*       dealPtr,
    LGMCalParm* CalReq,
    Date        tNow,
    String      ycName,
    LGM_TS*     tsPtr,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double*))                     /* swaption/cap exponents (beta) */
{
    SrtGenEurPtr deal;
    Date         tEx, tPay;
    double       zeta, Gpay, df, aD;
    long         i, nPay, k;
    LGMErr       error;

    error = NULL;

    if (dealtype != GenEur)
        return ("wrong deal type in GenEurpayoff");
    deal = (SrtGenEur*)dealPtr; /* re-cast as correct deal type */

    tEx  = deal->tEx;
    zeta = LGMZetaFromTS(tEx, tNow, tsPtr); /* Get zeta at exercise date */

    for (k = 0; k <= nx; k++) /* Initialize array */
        payoff[0][k] = 0;

    nPay = deal->nPay;
    for (i = 0; i < nPay; i++) /* For each payment received on this exercise */
    {
        tPay = deal->tPay[i];
        Gpay = LGMGFromTS(tPay, tsPtr);      /* Get G at paydate */
        df   = swp_f_df(tNow, tPay, ycName); /* Get dis fac at paydate */
        if (df == SRT_DF_ERROR)
            return (error);

        aD = df * deal->Payment[i];              /* PV of payment */
        aD = aD * exp(-.5 * Gpay * Gpay * zeta); /* prefactor */

        for (k = 0; k <= nx; k++) /* for each gridpoint */
            payoff[0][k] = payoff[0][k] + aD * exp(Gpay * x[k] * reduction[0]); /* add payment */

    } /* end sum over payments */
}

/**********************************************************************/
static LGMErr MidAtpayoff(
    double**    payoff,
    long        nx,
    double*     x,
    double*     reduction,
    LGMDealType dealtype,
    void*       dealPtr,
    LGMCalParm* CalReq,
    Date        tNow,
    String      ycName,
    LGM_TS*     tsPtr,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double*))                     /* swaption/cap exponents (beta) */
{
    SrtSimMidAtPtr MidAt;
    double *       zetaex, *Gpay, *Gset;
    double *       aD, *aredD, *StrikeD;
    double **      prefac, **postfac;
    double         df, argum;
    Date           tPay, tset, tEx;
    long           j, j0, jex, nEx, i, i0, ipay, ifirst, nPay, k;
    LGMErr         error;

    error = NULL;

    if (dealtype != SimpleMidAt)
        return ("wrong deal type in MidAtpayoff");
    MidAt = (SrtSimMidAt*)dealPtr; /* re-cast as correct deal type */

    /* Allocate space */
    j0   = MidAt->FirstExer;    /* exercises j = j0, ..., nEx-1 */
    nEx  = MidAt->nEx - j0;     /* Number of exercises (excluding ones before today+eod) */
    i0   = MidAt->FirstPay[j0]; /* first possible payment */
    nPay = MidAt->nPay - i0;    /* number of payments that could be received */

    zetaex  = NULL;
    StrikeD = NULL;
    aredD   = NULL;
    Gset    = NULL;
    aD      = NULL;
    Gpay    = NULL;
    prefac  = NULL;
    postfac = NULL;

    zetaex  = (double*)srt_calloc(nEx, sizeof(double));  /* zeta(tEx_j), j=0,...,nEx-1 */
    StrikeD = (double*)srt_calloc(nEx, sizeof(double));  /* Strike[j]*D(tset_j), j=0,...,nEx-1 */
    aredD   = (double*)srt_calloc(nEx, sizeof(double));  /* Ared[j]*D(tfirst_j), j=0,...,nEx-1 */
    Gset    = (double*)srt_calloc(nEx, sizeof(double));  /* G(tset_j), j=0,...,nEx-1 */
    Gpay    = (double*)srt_calloc(nPay, sizeof(double)); /* G(tPay_i), i=0,...,nPay-1 */
    aD      = (double*)srt_calloc(nPay, sizeof(double)); /* a[i]D(tPay_i), i=0,...,nPay-1 */
    prefac  = dmatrix(0, nx, 0, nPay - 1); /* exp(Gpay[i]*x[k]), k=0,..,nx, i=0,...,nPay-1 */
    postfac = dmatrix(
        0, nEx - 1, 0, nPay - 1); /* exp(-.5*zetaex[j]*Gpay[i]^2), j=0,...,nEx-1, i=0,...,nPay-1 */

    if (zetaex == NULL || StrikeD == NULL || aredD == NULL || Gset == NULL || Gpay == NULL ||
        aD == NULL || prefac == NULL || postfac == NULL)
    {
        srt_free(zetaex);
        srt_free(StrikeD);
        srt_free(aredD);
        srt_free(Gset);
        srt_free(aD);
        srt_free(Gpay);
        free_dmatrix(prefac, 0, nx, 0, nPay - 1);
        free_dmatrix(postfac, 0, nEx - 1, 0, nPay - 1);
        return ("alloc. failed in LGMautocal");
    }

    /* compute discounted cash flows and G at the pay dates */
    for (i = 0; i < nPay; i++)
    {
        ipay = i0 + i; /* Note: payments are offset by -i0 from payments in MidAt structure*/
        tPay = MidAt->tPay[ipay];
        df   = swp_f_df(tNow, tPay, ycName);
        if (df == SRT_DF_ERROR)
            return (error);
        aD[i]   = df * (MidAt->Payment[ipay]);
        Gpay[i] = LGMGFromTS(tPay, tsPtr);
    }

    /* compute zeta at exercise dates, G and StrikeD at settlement dates, and aredD at first
     * paydates */
    for (j = 0; j < nEx; j++)
    {
        jex       = j0 + j; /* Note: exer's are offset by -j0 from exer's in MidAt structure*/
        tEx       = MidAt->tEx[jex];
        zetaex[j] = LGMZetaFromTS(tEx, tNow, tsPtr);

        tset = MidAt->tStart[jex];
        df   = swp_f_df(tNow, tset, ycName);
        if (df == SRT_DF_ERROR)
            return (error);
        StrikeD[j] = df * MidAt->Strike[jex];
        Gset[j]    = LGMGFromTS(tset, tsPtr);

        i    = MidAt->FirstPay[jex];
        tPay = MidAt->tPay[i];
        df   = swp_f_df(tNow, tPay, ycName);
        if (df == SRT_DF_ERROR)
            return (error);
        aredD[j] = df * MidAt->RedFirstPay[jex];
    }

    /* compute postfactors */
    for (i = 0; i < nPay; i++)
    {
        argum = -0.5 * Gpay[i] * Gpay[i];
        for (j = 0; j < nEx; j++)
            postfac[j][i] = exp(argum * zetaex[j]);
    }

    /* compute payoff table */
    for (j = 0; j < nEx; j++)
    {
        jex = j0 + j;
        if (j == 0 || reduction[j] > 1.1 * reduction[j - 1]) /* recompute prefactors if needed */
        {
            ifirst = MidAt->FirstPay[jex] - i0;
            for (i = ifirst; i < nPay; i++)
            {
                for (k = 0; k <= nx; k++)
                    prefac[k][i] = exp(Gpay[i] * x[k] * reduction[j]);
            }
        }
        for (k = 0; k <= nx; k++)
        {
            i = MidAt->FirstPay[jex] - i0;
            payoff[j][k] =
                -StrikeD[j] *
                exp(Gset[j] * (reduction[j] * x[k] - 0.5 * Gset[j] * zetaex[j])); /* Strike */
            payoff[j][k] = payoff[j][k] +
                           (aD[i] - aredD[j]) * prefac[k][i] * postfac[j][i]; /* first payment */
            i++;
            for (; i < nPay; i++)
                payoff[j][k] =
                    payoff[j][k] + aD[i] * prefac[k][i] * postfac[j][i]; /* other payments */
        }
    }

    if (MidAt->PayRec == SRT_PAYER) /* Reverse sign for payers */
    {
        for (j = 0; j < nEx; j++)
            for (k = 0; k <= nx; k++)
                payoff[j][k] = -payoff[j][k];
    }

    /* Free and return */
    srt_free(zetaex);
    srt_free(StrikeD);
    srt_free(aredD);
    srt_free(Gset);
    srt_free(aD);
    srt_free(Gpay);
    free_dmatrix(prefac, 0, nx, 0, nPay - 1);
    free_dmatrix(postfac, 0, nEx - 1, 0, nPay - 1);
    return (error);
}

/**********************************************************************/
static LGMErr GenMidAtpayoff(
    double**    payoff,
    long        nx,
    double*     x,
    double*     reduction,
    LGMDealType dealtype,
    void*       dealPtr,
    LGMCalParm* CalReq,
    Date        tNow,
    String      ycName,
    LGM_TS*     tsPtr,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double*))                     /* swaption/cap exponents (beta) */
{
    SrtGenMidAtPtr deal;
    Date           tEx, tPay;
    double         zeta, Gpay, df, aD;
    long           j, j0, jex, nEx, i, nPay, k;
    LGMErr         error;

    error = NULL;

    if (dealtype != GenMidAt)
        return ("wrong deal type in GenMidAtpayoff");
    deal = (SrtGenMidAt*)dealPtr; /* re-cast as correct deal type */

    j0  = deal->FirstExer;
    nEx = deal->nEx - j0;

    for (j = 0; j < nEx; j++) /* For each exercise */
    {
        jex  = j0 + j;
        tEx  = deal->tEx[jex];
        zeta = LGMZetaFromTS(tEx, tNow, tsPtr); /* Get zeta at exercise date */

        for (k = 0; k <= nx; k++) /* Initialize array */
            payoff[j][k] = 0;

        nPay = deal->nPay[jex];
        for (i = 0; i < nPay; i++) /* For each payment received on this exercise */
        {
            tPay = deal->tPay[jex][i];
            Gpay = LGMGFromTS(tPay, tsPtr);      /* Get G at paydate */
            df   = swp_f_df(tNow, tPay, ycName); /* Get dis fac at paydate */
            if (df == SRT_DF_ERROR)
                return (error);

            aD = df * deal->Payment[jex][i];         /* PV of payment */
            aD = aD * exp(-.5 * Gpay * Gpay * zeta); /* prefactor */

            for (k = 0; k <= nx; k++) /* for each gridpoint */
                payoff[j][k] =
                    payoff[j][k] + aD * exp(Gpay * x[k] * reduction[j]); /* add payment */

        } /* end sum over payments */
    }     /* end loop over exercises */
    return error;
}

/**********************************************************************/

static double GetNormalVol(
    double      mat,
    Date        start,
    Date        end,
    double      fwd,
    double      strike,
    SRT_Boolean isCap,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*),
    LGMErr (*GetBeta)(Date, Date, double*))
{
    LGMErr error1;
    Err    error2;
    double beta;
    double inp_vol;
    double norm_vol;
    double std, nstd;
    double bs_price;

    if (mat < 1.0 / 365)
    {
        mat = 1.0 / 365;
    }

    if (fwd <= 0.0 || strike <= 0.0)
    {
        return 1.0e-08;
    }

    error1 = GetVol(start, end, fwd, isCap, &inp_vol);
    error1 = GetBeta(start, end, &beta);

    std  = pow(fwd, beta) * inp_vol * sqrt(mat);
    nstd = (strike - fwd) / std;

    error1 = GetVol(start, end, strike, isCap, &inp_vol);

    if (fabs(beta) < 1.0e-04)
    {
        return inp_vol;
    }

    if (fabs(nstd) > 5.0)
    {
        strike = fwd + 5.0 * std;
    }

    /*	Here we assume that input is lognormal */

    bs_price = srt_f_optblksch(fwd, strike, inp_vol, mat, 1.0, SRT_CALL, SRT_PREMIUM);

    error2 = srt_f_optimpvol(bs_price, fwd, strike, mat, 1.0, SRT_CALL, SRT_NORMAL, &norm_vol);

    return norm_vol;
}

static LGMErr BerInvFltpayoff(
    /* Output */
    double** payoff,
    /* Grid in the state variable */
    long    nx,
    double* x,
    /* Unused */
    double* reduc,
    /* Deal information */
    LGMDealType dealtype,
    void*       dealPtr,
    /* Used to extract the required caplet vol method */
    LGMCalParm* CalReq,
    /* Today's market information */
    Date    tNow,
    String  ycName,
    LGM_TS* tsPtr,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*),
    LGMErr (*GetBeta)(Date, Date, double*))
{
    SrtCallInvFlt* ptr;
    LGMErr         error = NULL;
    double         gear;
    long           j0, i0, nEx, nCpn, i, j, k;
    double*        zEx    = NULL;
    Date*          tFix   = NULL;
    double*        zFix   = NULL;
    double*        dSet   = NULL;
    double*        Gset   = NULL;
    double*        dStart = NULL;
    double*        dPay   = NULL;
    double*        Gstart = NULL;
    double*        Gpay   = NULL;
    double*        beta   = NULL;
    double*        payCvg = NULL;
    double*        refCvg = NULL;
    Date           t0, t1, tSpot, tStart, tEnd;
    double         YrtoExpFromtEx, YrtoExpFromtNow;
    long           ifirst;
    SrtCrvPtr      yldcrv;
    long           spotlag;
    double         fwd_libor, caplet_strike, caplet_pv, caplet_pv_today, coupon_pv;
    double         total_vol, total_var, model_var, caplet_vol;
    double         discount;
    int            switch_fwd;

    if (dealtype != CallInvFlt)
    {
        return "wrong deal type in BerInvFltpayoff";
    }
    ptr = (SrtCallInvFlt*)dealPtr; /* re-cast as correct deal type */

    yldcrv  = lookup_curve(ycName);
    spotlag = get_spotlag_from_curve(yldcrv);

    /* ALLOCATE SPACE */

    j0   = ptr->FirstEx;
    nEx  = ptr->nEx - j0;
    i0   = ptr->iSet[j0];
    nCpn = ptr->nCpn - i0;

    zEx  = (double*)srt_calloc(nEx, sizeof(double));
    dSet = (double*)srt_calloc(nEx, sizeof(double));
    Gset = (double*)srt_calloc(nEx, sizeof(double));

    tFix   = (Date*)srt_calloc(nCpn, sizeof(Date));
    zFix   = (double*)srt_calloc(nCpn, sizeof(double));
    dStart = (double*)srt_calloc(nCpn, sizeof(double));
    Gstart = (double*)srt_calloc(nCpn, sizeof(double));
    dPay   = (double*)srt_calloc(nCpn, sizeof(double));
    Gpay   = (double*)srt_calloc(nCpn, sizeof(double));
    payCvg = (double*)srt_calloc(nCpn, sizeof(double));
    refCvg = (double*)srt_calloc(nCpn, sizeof(double));

    if (zEx == NULL || dSet == NULL || Gset == NULL || tFix == NULL || dStart == NULL ||
        dPay == NULL || Gstart == NULL || Gpay == NULL || payCvg == NULL || refCvg == NULL)
    {
        error = "alloc. failed in BerInvFltpayoff";
        goto cleanup;
    }

    /* FILL ARRAYS */

    for (j = 0; j < nEx; j++)
    {
        zEx[j]  = LGMZetaFromTS(ptr->tEx[j + j0], tNow, tsPtr);
        dSet[j] = swp_f_df(tNow, ptr->tSet[j + j0], ycName);
        Gset[j] = LGMGFromTS(ptr->tSet[j + j0], tsPtr);
    }

    for (i = 0; i < nCpn; i++)
    {
        t0        = ptr->tCpnStart[i + i0];
        t1        = ptr->tCpnPay[i + i0];
        tFix[i]   = add_unit(t0, -spotlag, SRT_BDAY, NO_BUSDAY_CONVENTION);
        zFix[i]   = LGMZetaFromTS(tFix[i], tNow, tsPtr);
        dStart[i] = swp_f_df(tNow, t0, ycName);
        dPay[i]   = swp_f_df(tNow, t1, ycName);
        Gstart[i] = LGMGFromTS(t0, tsPtr);
        Gpay[i]   = LGMGFromTS(t1, tsPtr);
        payCvg[i] = ptr->cvg[i + i0];
        refCvg[i] = ptr->lcvg[i + i0];
    }

    for (j = 0; j < nEx; j++)
    {
        for (k = 0; k <= nx; k++)
        {
            payoff[j][k] = 0.0;
        }
    }

    ptr->fwdVolCpn = nCpn;
    ptr->fwdVolEx  = nEx;
    ptr->fwdVolMat = dmatrix(0, nCpn, 0, nEx);

    if (ptr->fwdVolMat == NULL)
    {
        error = "alloc. failed in BerInvFltpayoff";
        goto cleanup;
    }

    for (i = 0; i < nCpn; i++)
    {
        ptr->fwdVolMat[i + 1][0] = tFix[i];
    }

    for (j = 0; j < nEx; j++)
    {
        ptr->fwdVolMat[0][j + 1] = ptr->tEx[j + j0];
    }

    /* PAYOFF CALCULATION */

    /* For caplets, model fwd (0) or 1+cvg*fwd (1) */
    if (CalReq->CapletVolMeth == MODEL || CalReq->CapletVolMeth == MARKET)
    {
        switch_fwd = 1;
    }
    else
    {
        switch_fwd = 0;
    }

    /* Loop on exercises */
    for (j = 0; j < nEx; j++)
    {
        /* First coupon upon exercise */
        ifirst = ptr->iSet[j + j0] - i0;

        /* Loop on coupons upon exercise */
        for (i = ifirst; i < nCpn; i++)
        {
            /* Time from exercise to coupon fix */
            YrtoExpFromtEx = DMAX(1.0e-04, (tFix[i] - ptr->tEx[j + j0]) * YEARS_IN_DAY);
            /* Time from now to coupon fix */
            YrtoExpFromtNow = (tFix[i] - tNow) * YEARS_IN_DAY;

            /* Gearing of coupon */
            gear = ptr->gear[i + i0];

            /* Caplet strike */
            caplet_strike = ptr->cap_str[i + i0];

            /* If caplet expiry is less than 25 days, we price it with a vol of 0 */
            /*	Also, if caplet strike is irrelevant, we price it with a vol of 0 */
            if (tFix[i] - ptr->tEx[j + j0] > 25 && caplet_strike < 1 && caplet_strike > 1.0e-08)
            {
                /* Forward vol */
                switch (CalReq->CapletVolMeth)
                {
                case MODEL:

                    caplet_vol = (Gstart[i] - Gpay[i]) * sqrt((zFix[i] - zEx[j]) / YrtoExpFromtEx);
                    break;

                case SLIDING2:

                    tSpot      = add_unit(tNow, spotlag, SRT_BDAY, SUCCEEDING);
                    tStart     = tSpot + tFix[i] - ptr->tEx[j + j0];
                    tEnd       = tStart + ptr->tCpnPay[i + i0] - ptr->tCpnStart[i + i0];
                    caplet_vol = GetNormalVol(
                        (tFix[i] - ptr->tEx[j + j0]) * YEARS_IN_DAY,
                        tStart,
                        tEnd,
                        (swp_f_df(tNow, tStart, ycName) / swp_f_df(tNow, tEnd, ycName) - 1.0) /
                            refCvg[i],
                        caplet_strike,
                        SRT_TRUE,
                        GetVol,
                        GetBeta);
                    break;

                case CONVERGING:

                    caplet_vol = GetNormalVol(
                        (tFix[i] - tNow) * YEARS_IN_DAY,
                        ptr->tCpnStart[i + i0],
                        ptr->tCpnPay[i + i0],
                        (dStart[i] / dPay[i] - 1.0) / refCvg[i],
                        caplet_strike,
                        SRT_TRUE,
                        GetVol,
                        GetBeta);
                    break;

                case MARKET:

                    /* Get the total normal vol */
                    total_vol = GetNormalVol(
                        (tFix[i] - tNow) * YEARS_IN_DAY,
                        ptr->tCpnStart[i + i0],
                        ptr->tCpnPay[i + i0],
                        (dStart[i] / dPay[i] - 1.0) / refCvg[i],
                        caplet_strike,
                        SRT_TRUE,
                        GetVol,
                        GetBeta);

                    /* Get today's PV */
                    caplet_pv_today = srt_f_optblknrm(
                        (dStart[i] / dPay[i] - 1.0) / refCvg[i],
                        caplet_strike,
                        total_vol,
                        YrtoExpFromtNow,
                        dPay[i] * refCvg[i],
                        SRT_CALL,
                        SRT_PREMIUM);

                    /* Get the total lognormal vol of 1 + cvg * Libor */
                    error = srt_f_optimpvol(
                        caplet_pv_today,
                        dStart[i] / dPay[i],
                        1.0 + refCvg[i] * caplet_strike,
                        YrtoExpFromtNow,
                        dPay[i],
                        SRT_CALL,
                        SRT_LOGNORMAL,
                        &total_vol);
                    if (error)
                        return error;

                    total_var = total_vol * total_vol * YrtoExpFromtNow;

                    /* Calculate the model vol */
                    model_var = (Gstart[i] - Gpay[i]) * (Gstart[i] - Gpay[i]) * zEx[j];

                    /* Get the forward vol */
                    caplet_vol = sqrt((total_var - model_var) / YrtoExpFromtEx);

                    break;
                }
            }
            else
            {
                caplet_vol = 1.0e-08;
            }

            if (switch_fwd == 1)
            {
                ptr->fwdVolMat[i + 1][j + 1] = caplet_vol * dStart[i] / dPay[i] / refCvg[i];
            }
            else
            {
                ptr->fwdVolMat[i + 1][j + 1] = caplet_vol;
            }

            /* Loop on nodes */
            for (k = 0; k <= nx; k++)
            {
                /* Forward libor fixed @ exercise j for coupon i */
                fwd_libor = dStart[i] / dPay[i];
                fwd_libor *=
                    exp((Gstart[i] - Gpay[i]) * x[k] * reduc[j + j0] -
                        0.5 * (Gstart[i] * Gstart[i] - Gpay[i] * Gpay[i]) * zEx[j]);
                fwd_libor -= 1.0;
                fwd_libor /= refCvg[i];

                /* Caplet */
                if (switch_fwd == 1)
                {
                    if (1.0 + refCvg[i] * fwd_libor > 0 && 1.0 + refCvg[i] * caplet_strike > 0)
                    {
                        caplet_pv = srt_f_optblksch(
                            1.0 + refCvg[i] * fwd_libor,
                            1.0 + refCvg[i] * caplet_strike,
                            caplet_vol,
                            YrtoExpFromtEx,
                            1.0,
                            SRT_CALL,
                            SRT_PREMIUM);
                    }
                    else
                    {
                        caplet_pv = refCvg[i] * (fwd_libor - caplet_strike);
                        if (caplet_pv < 0.0)
                        {
                            caplet_pv = 0.0;
                        }
                    }

                    caplet_pv *= gear;
                }
                else
                {
                    if (fwd_libor > 0 && caplet_strike > 0)
                    {
                        caplet_pv = srt_f_optblknrm(
                            fwd_libor,
                            caplet_strike,
                            caplet_vol,
                            YrtoExpFromtEx,
                            1.0,
                            SRT_CALL,
                            SRT_PREMIUM);
                    }
                    else
                    {
                        caplet_pv = fwd_libor - caplet_strike;
                        if (caplet_pv < 0.0)
                        {
                            caplet_pv = 0.0;
                        }
                    }

                    caplet_pv *= gear * refCvg[i];
                }

                /* Coupon */
                coupon_pv = ptr->a[i + i0] - gear * refCvg[i] * fwd_libor + caplet_pv;

                /* PV */
                discount = dPay[i] * exp(Gpay[i] * (x[k] * reduc[j + j0] - 0.5 * Gpay[i] * zEx[j]));
                payoff[j][k] += coupon_pv * discount;
            }
        }

        for (k = 0; k <= nx; k++)
        {
            /*	Subtract strike */
            payoff[j][k] -= ptr->strike[j + j0] * dSet[j] *
                            exp(Gset[j] * (x[k] * reduc[j + j0] - 0.5 * Gset[j] * zEx[j]));
        }
    }

    if (ptr->PayRec == SRT_PAYER)
    {
        for (j = 0; j < nEx; j++)
        {
            for (k = 0; k <= nx; k++)
            {
                payoff[j][k] *= -1;
            }
        }
    }

/* Free and return */
cleanup:

    srt_free(zEx);
    srt_free(dSet);
    srt_free(Gset);
    srt_free(tFix);
    srt_free(zFix);
    srt_free(dStart);
    srt_free(dPay);
    srt_free(Gstart);
    srt_free(Gpay);
    srt_free(payCvg);
    srt_free(refCvg);

    return error;
}

/* Computes a fra in LGM autocal */
LGMErr fra_in_LGMAutocal(
    String  ycName,
    LGM_TS* tsPtr,
    long    value_date,
    long    tNow,
    long    start,
    long    end,
    double  xred, /* state variable in LGM autocal */
    double* result)
{
    LGMErr error = NULL;
    double level, cvg;
    double dfstart, dfend;
    double zeta, G; /* parameters for the term structure of LGM autocal */

    /* computes the df start seen at value date */
    dfstart = swp_f_df(tNow, start, ycName);
    zeta    = LGMZetaFromTS(value_date, tNow, tsPtr); /* Get zeta at exercise date */
    G       = LGMGFromTS(start, tsPtr);
    dfstart *= exp(-0.5 * G * G * zeta);
    dfstart *= exp(G * xred);

    /* computes the df end seen at value date */
    dfend = swp_f_df(tNow, end, ycName);
    zeta  = LGMZetaFromTS(value_date, tNow, tsPtr); /* Get zeta at exercise date */
    G     = LGMGFromTS(end, tsPtr);
    dfend *= exp(-0.5 * G * G * zeta);
    dfend *= exp(G * xred);

    /* computes the coverage and level */
    cvg   = coverage(start, end, BASIS_ACT_360);
    level = cvg * dfend;

    /* Computes the fra */
    *result = (dfstart - dfend) / level;
    return error;
}

/**********************************************************************/
static LGMErr BerCapFltpayoff(
    double**    payoff,
    long        nx,
    double*     x,
    double*     reduc,
    LGMDealType dealtype,
    void*       dealPtr,
    LGMCalParm* CalReq,
    Date        tNow,
    String      ycName,
    LGM_TS*     tsPtr,
    LGMErr (*GetVol)(Date, Date, double, SRT_Boolean, double*), /* swaption/cap vols */
    LGMErr (*GetBeta)(Date, Date, double*))                     /* swaption/cap exponents (beta) */
{
    SrtCallCapFlt* ptr;
    LGMErr         error = NULL;
    double         gear;
    double         CAPLETMAT, CAPLETMAXPV, CAPLETMINPV;
    double         ATMVOL;

    /* arrays for intermediates */
    long     j0, i0, nEx, nCpn, i, j, k;
    double*  zEx    = NULL; /* zeta(tEx[j]), j = j0, ..., nEx-1 */
    double*  dPay   = NULL; /* d(tCpn[i]), i = i0,..., nCpn */
    double*  beta   = NULL; /* beta for cpn i, i = i0+1,..., nCpn */
    double*  cvgCpn = NULL; /* cvg for cpn i, i = i0+1,..., nCpn */
    double*  cvgR   = NULL; /* flting rate cvg for in cpn i, i = i0+1,..., nCpn */
    double*  Gpay   = NULL; /* G(tCpn[i]), i = i0, ..., nCpn */
    Date*    tFix   = NULL; /* fixing date tFix[i], i = i0+1,...,nCpn */
    double** postf  = NULL; /* exp(-0.5*Gpay[i]*Gpay[i]*zEx[j]), [j][i], i=iSet[j],...,nCpn */
    double** pref   = NULL; /* D[i]*exp(G[i]*x[k]), [i][k], i = iSet[j],...,nCpn, k=0,..., nx */

    Date      t0, t1;
    double    cpn;
    double    amaxi, amini, rate;
    double    zFix;
    long      ifirst;
    SrtCrvPtr yldcrv;
    long      spotlag;

    /* Initialisation */
    if (dealtype != CallCapFlt)
        return ("wrong deal type in BerCapFltpayoff");
    ptr  = (SrtCallCapFlt*)dealPtr; /* re-cast as correct deal type */
    gear = ptr->lvg;

    yldcrv  = lookup_curve(ycName);
    spotlag = get_spotlag_from_curve(yldcrv);

    /* Allocate space */
    j0   = ptr->FirstEx;   /* exercises j = j0, ..., ptr.nEx-1 */
    nEx  = ptr->nEx - j0;  /* number of exercises j = 0,..., nEx-1*/
    i0   = ptr->iSet[j0];  /* cpn dates i = i0, ..., ptr.nCpn */
    nCpn = ptr->nCpn - i0; /* number of cpns: i = 0, ..., nCpn */

    zEx    = (double*)srt_calloc(nEx, sizeof(double));      /* zeta[j0+j], j=0,...,nEx-1 */
    dPay   = (double*)srt_calloc(nCpn + 1, sizeof(double)); /* dPay[i0+i], i = 0,..., nCpn */
    beta   = (double*)srt_calloc(nCpn + 1, sizeof(double)); /* beta[i0+i], i = *,1,..., nCpn */
    cvgCpn = (double*)srt_calloc(nCpn + 1, sizeof(double)); /* cvg[i0+i], i = *,1,..., nCpn */
    cvgR   = (double*)srt_calloc(nCpn + 1, sizeof(double)); /* cvg[i0+i], i = *,1,..., nCpn */
    Gpay   = (double*)srt_calloc(nCpn + 1, sizeof(double)); /* Gpay[i0+i], i = 0,..., nCpn */
    tFix   = (Date*)srt_calloc(nCpn + 1, sizeof(Date));     /* tFix[i0+i], i = *,1,..., nCpn */

    postf = dmatrix(0, nEx - 1, 0, nCpn); /* [j][i] = [0,1,...,nEx-1][0,...,nCpn] */
    pref  = dmatrix(0, nCpn, 0, nx);      /* [i][k] = [0,1,...,nCpn][0,...,nx] */

    if (zEx == NULL || dPay == NULL || beta == NULL || cvgCpn == NULL || cvgR == NULL ||
        Gpay == NULL || tFix == NULL || postf == NULL || pref == NULL)
    {
        error = "alloc. failed in BerCapFltpayoff";
        goto cleanup;
    }

    t0      = ptr->tCpn[i0];
    dPay[0] = swp_f_df(tNow, t0, ycName);
    if (dPay[0] == SRT_DF_ERROR)
    {
        error = "dis factor failed";
        goto cleanup;
    }
    Gpay[0] = LGMGFromTS(t0, tsPtr);
    for (j = 0; j < nEx; j++)
    {
        zEx[j]      = LGMZetaFromTS(ptr->tEx[j + j0], tNow, tsPtr);
        postf[j][0] = exp(-0.5 * Gpay[0] * Gpay[0] * zEx[j]);
    }

    for (i = 1; i <= nCpn; i++)
    {
        t0      = ptr->tCpn[i + i0 - 1];
        t1      = ptr->tCpn[i + i0];
        dPay[i] = swp_f_df(tNow, t1, ycName);
        if (dPay[i] == SRT_DF_ERROR)
        {
            error = "dis factor failed";
            goto cleanup;
        }

        cvgCpn[i] = coverage(t0, t1, ptr->cpnBasis);
        cvgR[i]   = coverage(t0, t1, ptr->aBasis);
        Gpay[i]   = LGMGFromTS(t1, tsPtr);
        tFix[i]   = add_unit(t0, -spotlag, SRT_BDAY, NO_BUSDAY_CONVENTION);

        for (j = 0; j < nEx; j++)
            postf[j][i] = exp(-0.5 * Gpay[i] * Gpay[i] * zEx[j]);
    }

    /*PAYOFF INITIALISATION */
    for (j = 0; j < nEx; j++)
    {
        for (k = 0; k <= nx; k++)
            payoff[j][k] = 0.0; /*PAYOFF VALUE AT THE J THE EXERCISE DATE AND NODE K */
    }                           /*END OF PAYOFF INITIALISATION */

    for (j = 0; j < nEx; j++)
    {
        if (j == 0 || reduc[j] > SRT_MIN_REDUC_FACT * reduc[j - 1]) /* X IS RE-SCALED */
        {
            ifirst = ptr->iSet[j + j0] - i0; /* NEED TO RECOMPUTE X-DEP FACTORS */
            for (i = ifirst; i <= nCpn; i++) /* ONLY NEED TO FOR CPN DATES > ifirst */
            {
                for (k = 0; k <= nx; k++)
                    pref[i][k] = dPay[i] * exp(Gpay[i] * x[k] * reduc[j]);
            }

        } /*END OF if (j==0 || reduc[j]>SRT_MIN_REDUC_FACT*reduc[j-1]) */

        ifirst = ptr->iSet[j + j0] - i0;
        for (i = ifirst + 1; i <= nCpn; i++)
        {
            CAPLETMAT = DMAX(1e-4, (tFix[i] - ptr->tEx[j + j0]) * YEARS_IN_DAY);
            zFix      = LGMZetaFromTS(tFix[i], tNow, tsPtr);
            ATMVOL    = (Gpay[i - 1] - Gpay[i]) * sqrt((zFix - zEx[j]) / CAPLETMAT);

            for (k = 0; k <= nx; k++)
            {
                rate =
                    (pref[i - 1][k] * postf[j][i - 1] / (pref[i][k] * postf[j][i]) - 1.0) / cvgR[i];

                /*RE - INITIALISTION OF THE CPN AND SETTING OF THE STRIKE*/
                cpn         = 0;
                CAPLETMAXPV = 0;
                CAPLETMINPV = 0;

                amini = ptr->amin[i + i0];
                amaxi = ptr->amax[i + i0];

                CAPLETMINPV = srt_f_optblksch(
                    (1 + cvgCpn[i] * rate),
                    (1 + cvgCpn[i] * amini),
                    ATMVOL,
                    CAPLETMAT,
                    1.0, /*DISCOUNT FACTOR */
                    SRT_CALL,
                    SRT_PREMIUM);

                CAPLETMAXPV = srt_f_optblksch(
                    (1 + cvgCpn[i] * rate),
                    (1 + cvgCpn[i] * amaxi),
                    ATMVOL,
                    CAPLETMAT,
                    1.0, /*DISCOUNT FACTOR */
                    SRT_CALL,
                    SRT_PREMIUM);

                cpn = cvgCpn[i] * ptr->marg[i + i0] + gear * (CAPLETMINPV - CAPLETMAXPV);

                payoff[j][k] = payoff[j][k] + cpn * pref[i][k] * postf[j][i];
            }
        }
        for (k = 0; k <= nx; k++)
        {
            payoff[j][k] = payoff[j][k] + pref[nCpn][k] * postf[j][nCpn];
            payoff[j][k] = payoff[j][k] - ptr->strike[j + j0] * pref[ifirst][k] * postf[j][ifirst];
        }
    }

    if (ptr->PayRec == SRT_PAYER) /* Reverse sign for payers */
    {
        for (j = 0; j < nEx; j++)
        {
            for (k = 0; k <= nx; k++)
                payoff[j][k] = -payoff[j][k];
        }
    }

/* Free and return */
cleanup:
    srt_free(zEx);
    srt_free(dPay);
    srt_free(beta);
    srt_free(cvgCpn);
    srt_free(cvgR);
    srt_free(Gpay);
    srt_free(tFix);

    free_dmatrix(postf, 0, nEx - 1, 0, nCpn);
    free_dmatrix(pref, 0, nCpn, 0, nx);
    return (error);
}
