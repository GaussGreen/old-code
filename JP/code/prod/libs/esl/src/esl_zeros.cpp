/****************************************************************************/
/*      Zero price and zero bank.                                           */
/****************************************************************************/
/*      ZEROS.c                                                             */
/****************************************************************************/
/****************************************************************************
** IMPORTANT: This file contains conditionally compiled code depending upon
**            whether one is using IRX curves (i.e. -DESL_NEW_DATE
**            -DESL_NEW_CURVE) - IRX dates are a requisite for IRX curves.
**
**            The code is organized so that the conditionally compiled code is
**            placed near the top of the file, and common code (without any
**            conditional compilation) is placed near the end.
**
**            PLEASE respect this convention when modifying this file.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "esl_zeros.h"
#include "esl_date.h"
#include "esl_util.h"
#include "esl_error.h"
#include "esl_macros.h"

#include "string.h"

#ifdef ESL_NEW_CURVE
#include "irx/zerocurve.h"
#endif


#ifndef ESL_NEW_CURVE

/* =========compute zero at MatDate given zero rate structure ==============
   support both ESL_INTERP_LINEAR and ESL_INTERP_FLATFWD 
   special case3
   (1) MatDate < ZeroDates[0] or MatDate < ValueDate
        flat zero rate extrapolation of ZeroRates[0]
   (2) MatDate > ZeroDates[NbZero-1]:
        (a) ESL_INTERP_LINEAR
            flat zero rate extrapolation of ZeroRates[NbZero-1]
        (b) ESL_INTERP_FLATFWD
            flat fwd rate of ZeroRates[NbZero-2] and ZeroRates[NbZero-1]

  The code can handle the following special cases:
  (3) MatDate < ValueDate
  (4) ValueDate > ZeroRates[0] 
*/
double ZeroPriceBase(IRDate    MatDate   /** (I) Mat date of the zero       */
                 ,IRDate        ValueDate  /** (I) Value date                 */
                 ,int           NbZero     /** (I) Number of zeros in curve   */
                 ,IRDate const* ZeroDates  /** (I) maturity dates in zero crv */
                 ,double const* ZeroRates  /** (I) Zero rates                 */
                 ,ESL_INTERP_TYPE InterpType 
        )
{
    double Price = -999.99;
    double t,  ZR,
           t1, ZR1, Z1,
           t2, ZR2, Z2;
    int    idx;

    /* basic checks */
    if ((ZeroDates == NULL) || (ZeroRates == NULL) || NbZero <= 0) 
    {
        DR_Error("Null Zero Curve ");
        goto RETURN;
    }

    if (MatDate == ValueDate) 
    {
        return 1.0;
    }

    t   = Daysact(ValueDate, MatDate)/365.0;
    idx = GetDLOffset(NbZero, ZeroDates, MatDate, CbkHIGHER);

    /* a flat curve */
    if (NbZero == 1)
    {   
        Price = pow(1.0 + ZeroRates[0], -t);
        goto RETURN;
    }

    if (idx == 0 || (idx > 0 && ZeroDates[idx] == MatDate))
    {  /* (1) MatDate < first ZeroDate : flat zero rate extrapolation.
          (2) MatDate is at a curve point : no interpolation is needed.
        */
       Price = pow(1+ZeroRates[idx], -t);
       goto RETURN;
    }
    else if (idx<0) 
    {   /* MatDate > last ZeroDate: 
           ESL_INTERP_LINEAR: flat zero rate extrapolation
           ESL_INTERP_FLATFWD: use previous flat fwd rate
        */
        switch (InterpType)
        {
        case ESL_INTERP_LINEAR: /* Linear Zero Cpn */
            t1  = Daysact(ValueDate, ZeroDates[NbZero-1])/365.0;
            ZR1 = ZeroRates[NbZero-1];

            t2  = t;
            ZR2 = ZeroRates[NbZero-1];
            break;
        
        case ESL_INTERP_FLATFWD: /* Flat Fwd */    
            t1  = Daysact(ValueDate, ZeroDates[NbZero-2])/365.0;
            ZR1 = ZeroRates[NbZero-2];

            t2  = Daysact(ValueDate, ZeroDates[NbZero-1])/365.0;
            ZR2 = ZeroRates[NbZero-1];
            break;        
        
        default:
            DR_Error("Support only LINEAR and FLATFWD interp type");
            goto RETURN;
        }
    }
    else 
    {   /* MatDate is inside the curve */    
        t1  = Daysact(ValueDate, ZeroDates[idx-1])/365.0;
        ZR1 = ZeroRates[idx-1];

        t2  = Daysact(ValueDate, ZeroDates[idx])/365.0;
        ZR2 = ZeroRates[idx]; 
    }
    
    
    if (IS_EQUAL(t1, t2))
    {
        DR_Error("Bad Zero Curve \n");
        goto RETURN;
    }


    switch (InterpType)
    {
    case ESL_INTERP_LINEAR: /* Linear Zero Cpn */    
        if (linterp(t, &ZR,t1,  t2,ZR1, ZR2) == FAILURE) 
        {
            goto RETURN;
        }
        Price = pow(1.0 + ZR, -t);
        break;
        
    case ESL_INTERP_FLATFWD: /* Flat Fwd */

        Z1  = pow(1.0 + ZR1, -t1);
        Z2  = pow(1.0 + ZR2, -t2);
        Price = Z1 * pow(Z2/Z1, (t-t1)/(t2-t1));
        break;
        
    default:
        DR_Error("Support only LINEAR and FLATFWD interp type");
        goto RETURN;
    }


RETURN:
   
    if (Price < 0.0)
    {
        DR_Error("ZeroPrice: failed.");
    }

    assert(DR_IS_FINITE(Price)); /* not checking in release mode for speed */

    return (Price);

} 


/*****  ZeroPrice  ****************************************************/
/**
 *      DEPRECATED!!! - Use GetZeroPrice instead.
 *
 *      Calculates a deterministic zero bond price
 *      Returns -999.99 if failure
 *
 *      Supports
 *      (EslGetZeroInterpolation() == ESL_INTERP_LINEAR) Linear Zero Cpn:
 *      -- linear zero cpn interp
 *      -- flat zero cpn extrapolation on both sides
 *
 *      (EslGetZeroInterpolation() == ESL_INTERP_FLATFWD) Flat Fwd:
 *      -- flat fwd interp and extrapolation
 */
double ZeroPrice(IRDate          MatDate   /** (I) Mat date of the zero       */
                 ,IRDate        ValueDate  /** (I) Value date                 */
                 ,int           NbZero     /** (I) Number of zeros in curve   */
                 ,IRDate const* ZeroDates  /** (I) maturity dates in zero crv */
                 ,double const* ZeroRates  /** (I) Zero rates                 */)
{
    return ZeroPriceBase(MatDate, ValueDate, NbZero, ZeroDates, ZeroRates, EslGetZeroInterpolation());
}


#endif  /* ifndef ESL_NEW_CURVE */


/*****************************************************************************
** MIXED IRX/NON-IRX CURVE CODE
*****************************************************************************/

/* Answers a new empty zero curve that has been "properly" initialized.     */
/* IMPORTANT:  Caller must free returned T_CURVE via ZeroCurveFree.         */
T_CURVE* ZeroCurveMakeEmpty(int numElems)
{
#ifdef ESL_NEW_CURVE
    return irxZeroCurveMakeEmpty(numElems);
#else
    T_CURVE *crv = (T_CURVE*)malloc(sizeof(T_CURVE));
    if (crv == NULL)
    {
        DR_Error("Cannot allocate T_CURVE.");
    }
    else
    {
        crv->NbZero = 0;
        crv->ValueDate = crv->Today = INVALID_DATE;
        crv->SpotDays = -1;
        crv->InterpType = ESL_INTERP_LINEAR  /* SRM3_LINEAR_INTERP */;
        crv->MMB = -1;
        strcpy(crv->SwapDCC, "365");
        crv->SwapFreq = ESL_FREQ_ANNUAL;
    }

    return crv;
#endif
}

int ZeroCurvePopulateFromRates(
            IRDate          baseDate,
            int             numElems,
            const IRDate   *zeroDates,
            const double   *zeroRates,
            T_CURVE* crv)
{
#ifdef ESL_NEW_CURVE
return irxZeroCurveConstructFromRates(crv, baseDate, numElems, zeroDates, zeroRates,  IRX_ANNUAL_RATE, IRX_ACT_365F);
#else
    int datesToSkip = 0;
    int i = 0;
    int status = FAILURE;

    if (crv != NULL)
    {
        crv->ValueDate = baseDate;

        /* Today's date and spot days are not available */
        crv->Today = crv->ValueDate;
        crv->SpotDays  = 0;

        crv->InterpType = ESL_INTERP_LINEAR;
        crv->MMB = 360;
        strcpy(crv->SwapDCC, "365");
        crv->SwapFreq = ESL_FREQ_ANNUAL;

        while (datesToSkip < numElems && zeroDates[datesToSkip] < baseDate) {
            ++datesToSkip;
        }

        crv->NbZero = numElems - datesToSkip;

        if (crv->NbZero >= MAXNBDATE) {
            crv->NbZero = MAXNBDATE;
        }

        for (i=0; i < crv->NbZero; ++i) {
            crv->ZeroDate[i] = zeroDates[datesToSkip + i];
            crv->Zero[i] = zeroRates[datesToSkip + i];
        }
        status = SUCCESS;
    }
    else
    {
        DR_Error("ZeroCurvePopulateFromRates: Input curve is NULL.");
    }

    return status;
#endif
}

/* Zero curve (T_CURVE) constructor.                                        */
/* IMPORTANT:  Input zero rates are assumed to be annually compounded.      */
/*             Also, DCC is hard-coded as ACT/365F.                         */
/* Caller is responsible for freeing returned T_CURVE via ZeroCurveFree.    */
T_CURVE* ZeroCurveMakeFromRates(
            IRDate          baseDate,
            int             numElems,
            const IRDate   *zeroDates,
            const double   *zeroRates)
{
#ifdef ESL_NEW_CURVE
    return irxZeroCurveMakeFromRates(
                baseDate,
                numElems,
                zeroDates,
                zeroRates,
                IRX_ANNUAL_RATE,
                IRX_ACT_365F);
#else

    T_CURVE *crv = ZeroCurveMakeEmpty(numElems);

    if (ZeroCurvePopulateFromRates(baseDate,
                                numElems,
                                zeroDates,
                                zeroRates,
                                crv) != SUCCESS) {
        DR_Error("ZeroCurvePopulateFromRates failed.");
        crv = NULL;
    }

    return crv;
#endif
}

/* Answers a copy of the input T_CURVE.                                     */ 
/* Caller is responsible for freeing returned T_CURVE via ZeroCurveFree.    */
T_CURVE* ZeroCurveCopy(const T_CURVE *src)
{
    T_CURVE *copy = NULL;

#ifdef ESL_NEW_CURVE
    copy = irxZeroCurveCopy(src);
#else
    if (src != NULL && (copy = ZeroCurveMakeEmpty(src->NbZero)) != NULL) {
        int i;

        for (i=0; i < src->NbZero; ++i) {
            copy->ZeroDate[i] = src->ZeroDate[i];
            copy->Zero[i] = src->Zero[i];
        }

        copy->NbZero = src->NbZero;
        copy->Today = src->Today;
        copy->ValueDate = src->ValueDate;
        copy->SpotDays = src->SpotDays;
        copy->InterpType = src->InterpType;
        copy->MMB = src->MMB;
        strcpy(copy->SwapDCC, src->SwapDCC);
        copy->SwapFreq = src->SwapFreq;
    }
#endif  /* ESL_NEW_CURVE */

    return copy;
}


/*  Properly dispose of a T_CURVE created via ZeroCurveMakeEmpty,           */
/*  ZeroCurveMakeFromRates, or ZeroCurveCopy.                               */ 
void ZeroCurveFree(T_CURVE *crv)
{
#ifdef ESL_NEW_CURVE
    irxZeroCurveFree(crv);
#else
    if (crv != NULL)
    {
        free(crv);
    }
#endif
}

void ZeroCurveInit(T_CURVE *crv)
{
#ifdef ESL_NEW_CURVE
    irxZeroCurveInit(crv);
#else
    crv->NbZero = 0;
    crv->InterpType = -1;
#endif
}

IRDate GetZeroCurveBaseDate(const T_CURVE *crv)
{
#ifdef ESL_NEW_CURVE
    return irxZeroCurveBaseDate(crv);
#else
    return crv->ValueDate;
#endif
}

IRDate GetZeroCurveLastDate(const T_CURVE* crv)
{
#ifdef ESL_NEW_CURVE
    return irxZeroCurveLastDate(crv);
#else
    return crv != NULL? crv->ZeroDate[crv->NbZero-1]: -1;
#endif
}


int IsZeroCurveEmpty(const T_CURVE* crv)
{
    if (crv == NULL)
        return TRUE;

#ifdef ESL_NEW_CURVE
    return crv->numItems <= 1L;
#else
    return crv->NbZero == 0L;
#endif

}

int GetZeroCurveNumPoints(const T_CURVE* crv)
{
    if (crv == NULL)
        return 0;

#ifdef ESL_NEW_CURVE
    return crv->numItems - 1;
#else
    return crv->NbZero;
#endif
}

int printZeroRates(FILE* fp, const T_CURVE* zc, const char* formatstr)
{
    int i;
    
#ifdef ESL_NEW_CURVE
    IRDate basedate = GetZeroCurveBaseDate(zc);

    for (i=1; i < zc->numItems;++i)
    {
        IRDate enddate = zc->startDates[i];
        double df = zc->prices[i];
        double zr;
       
        if (irxDiscountToRate(df, basedate, enddate, IRX_ACT_365F, IRX_ANNUAL_RATE, &zr) != SUCCESS)
        {
            DR_Error("cannot convert discount to rate");
            return FAILURE;
        }
        fprintf(fp,formatstr, YMDDateFromIRDate(enddate), zr*100);
    }
#else
    for (i=0; i< zc->NbZero;++i)
    {
        IRDate enddate = zc->ZeroDate[i];
        fprintf(fp,formatstr, YMDDateFromIRDate(enddate), zc->Zero[i]*100);
    }
#endif
    return SUCCESS;
}


/******************************************************************************
 *      Calculates a deterministic zero bond price
 *      Returns -999.99 if failure
 *
 *      Supports
 *      -- smooth-forward interp
 *      -- flat-forward interp
 *
 */
double   GetZeroPrice(IRDate          MatDate  /** (I) Mat date of the zero  */
                     ,T_CURVE const* crv)     /** (I) Zero curve            */
{
#ifdef ESL_NEW_CURVE
    double Price;
    if (irxZeroPrice(crv, MatDate, &Price) != SUCCESS)
        Price = -999.99;

    if (Price < 0.0)
        DR_Error("ZeroPrice: failed.");
    return Price;
#else
    return ZeroPriceBase(MatDate, crv->ValueDate, crv->NbZero, crv->ZeroDate, crv->Zero, (ESL_INTERP_TYPE)crv->InterpType); 
#endif
}


/**
 *     Calculate zero rate & price for a specific maturity from the zero
 *     curve stored in t_curve style (deterministic).
 *
 */
int  GetZeroPriceRate(
         double   *OutZeroRate,  /** (O) Zero rate                   */
         double   *OutZeroPrice, /** (O) Zero price                  */
         IRDate    MatDate,      /** (I) Maturity of the zero        */
         T_CURVE const* crv)
{
    IRDate startDate = GetZeroCurveBaseDate(crv);
    double t;

    if (MatDate == startDate)
    {
        *OutZeroPrice = 1.0;
        *OutZeroRate = 0.0;
        return SUCCESS;
    }
#if 0
    /* NOTE: Some tests (e.g. binturboopflows_t deal80.dat) fail this check.  */
    /*       This will be an issue when using IRX curves as irxZeroPrice will */
    /*       fail if MatDate <= zero curve base date.                         */
    else if (MatDate < startDate)
    {
        DR_Error("GetZeroPriceRate: maturity date (%ld) is before zero curve "
                 "base date (%ld).\n",
                 YMDDateFromIRDate(MatDate),
                 YMDDateFromIRDate(startDate));
        return FAILURE;
    }
#endif

    t = Daysact(startDate, MatDate) / 365.0;

#ifdef ESL_NEW_CURVE
    if (irxZeroPrice(crv, MatDate, OutZeroPrice) != SUCCESS)
#else
    *OutZeroPrice = GetZeroPrice(MatDate, crv);
    if (*OutZeroPrice < TINY)
#endif
        return FAILURE;

    *OutZeroRate = pow(*OutZeroPrice, -1.0/t ) - 1.0;

    return SUCCESS;
}


/**
 *  With IRX zero curve, one cannot set interpolation type explicitly.  *
 *  Rather, one needs to build the zero curve with the desired
 *  interpolation type.
 */

/* IMPORTANT: ALL NEW CODE SHOULD NOT ACCESS THESE FLAGS DIRECTLY!!!    *
**            INSTEAD, USE THE FUNCTIONS BELOW.                         */
ESL_INTERP_TYPE ZeroInterpTypeFlag = ESL_INTERP_FLATFWD; 
ESL_INTERP_TYPE ZeroInterpTypeFlagStub = ESL_INTERP_FLATFWD;    


/* To enable comparison testing with non-IRX enabled products that use linear *
** zero interpolation.                                                        */
static int forceFlatFwd()
{
    static int _forceFlatFwd = -1;
    if (_forceFlatFwd < 0)
        _forceFlatFwd = getenv("FORCE_FLAT_FWD") != NULL ? 1 : 0;
    return _forceFlatFwd;
}


ESL_INTERP_TYPE  EslGetZeroInterpolation()
{
#ifdef ESL_NEW_CURVE
    /* NOTE:  We really should not even be calling this function when using   *
     *        IRX curves, and we certainly must not use the value returned    *
     *        for any calculations.  We just answer the "default" so that we  *
     *        don't have to change all client (e.g. product) code.            */ 
    return ESL_INTERP_FLATFWD;
#else
    if (forceFlatFwd())
        return ESL_INTERP_FLATFWD;
    return ZeroInterpTypeFlag;
#endif
}



void EslSetZeroInterpolation(ESL_INTERP_TYPE t)
{
#ifdef ESL_NEW_CURVE
    if (!(t != ESL_INTERP_FLATFWD || forceFlatFwd()))
    {
        DR_Error("Illegal attempt to explicitly set interpolation type for IRX zero curve.");
        abort();
    } 
#else
    ZeroInterpTypeFlag = ZeroInterpTypeFlagStub = t;
#endif
}


/*****  ExtendTCurve  ********************************************************/
/**
 *      Flat extend zero curve between fromDate and toDate dates
 *
 *      FIXME:  Not implemented for IRX curves.  See, however,
 *              irxZeroCurveExtendToToday.
 */
int  ExtendTCurve(T_CURVE *tc, IRDate fromDate, IRDate toDate)
{
    static const char* routine = "ExtendTCurve";

#ifdef ESL_NEW_CURVE
    if (toDate == fromDate || toDate == 0)
        return (int)irxZeroCurveSetAndExtendToToday(tc, fromDate);
    else
    {
        DR_Error("%s: Generalized version of this function is not supported "
                 "for IRX curves.", routine);
        return FAILURE;
    }
#else
    int      i, nbZero = tc->NbZero;
    int      idx = 0;
    double   T2Vrate   = tc->Zero[0];
    double   stmZero;
    double   AZero, ARate;
    IRDate   baseDate = tc->ValueDate;

    /* no need to extend if value date equals today */  /* JAC - What about toDate? */
    if (baseDate == fromDate )
    {
        /*
        if (toDate > tc->ZeroDate[tc->NbZero-1])
        {
            tc->NbZero += 1;
            tc->ZeroDate[tc->NbZero-1] = toDate;
            tc->Zero[tc->NbZero-1] = tc->Zero[tc->NbZero-2];
        }
        */

        return SUCCESS;
    }

    if ( Daysact(tc->ZeroDate[nbZero-1], fromDate) >= 0 )
    {
        DR_Error ("%s: the last zero Date %ld <= from date %ld\n",
                  routine,
                  YMDDateFromIRDate(tc->ZeroDate[nbZero-1]),
                  YMDDateFromIRDate(fromDate));
        return FAILURE;
    }
    
    stmZero   = pow((1.0 + T2Vrate), -Daysact(fromDate, baseDate)/365.0);

    if ( Daysact(baseDate, fromDate) > 0)
    {
        stmZero = GetZeroPrice(fromDate, tc);
        stmZero = 1. / stmZero;
        idx     = GetDLOffset(nbZero, tc->ZeroDate, fromDate, CbkHIGHER);

        if ( Daysact(fromDate, tc->ZeroDate[idx]) == 0)
            idx++;
    }
    
    /* adjust zero rates to be from from date */
   
    for (i=idx; i<nbZero; i++)
    {
        double T = Daysact(tc->ValueDate, tc->ZeroDate[i])/365.;
        AZero = pow( 1 + tc->Zero[i] , -T);

        if (AZero <= TINY)
        {
            DR_Error(   "%s: zero < TINY for date %ld\n",
                        routine, 
                        YMDDateFromIRDate(tc->ZeroDate[i]));
        }

        AZero *= stmZero;   /* convert it to a zero from fromDate */

        ARate = pow(AZero, -365.0/Daysact(fromDate, tc->ZeroDate[i])) - 1.0;

        tc->Zero[i] = ARate;
    }
    
    /* set today and value date to fromDate */
    tc->Today = tc->ValueDate = fromDate;
    if (toDate != 0 && tc->ZeroDate[tc->NbZero-1] < toDate)
    {
        tc->NbZero += 1;
        tc->ZeroDate[tc->NbZero-1] = toDate;
        tc->Zero[tc->NbZero-1] = tc->Zero[tc->NbZero-2];
    }

    return SUCCESS;
#endif
}


/*****************************************************************************
** GENERIC (NON-CONDITIONAL) CODE
*****************************************************************************/
void EslSetZeroInterpolationLinear()
{
    EslSetZeroInterpolation(ESL_INTERP_LINEAR);
}


void EslSetZeroInterpolationFlatFwd()
{
    EslSetZeroInterpolation(ESL_INTERP_FLATFWD);
}

/************************************************************************************
**
** IMPORTANT:  READ COMMENTS AT TOP OF FILE BEFORE ADDING/MODIFYING CODE TO THIS FILE. 
**
************************************************************************************/
