#ifndef ESL_ZEROS_DOT_H
#define ESL_ZEROS_DOT_H

#include "esl_date.h"
#include "esl_types.h"

#ifdef  __cplusplus
extern "C" {
#endif


/* Answers a new empty zero curve that has been "properly" initialized.     */
/* IMPORTANT:  Caller must free returned T_CURVE via ZeroCurveFree.         */
T_CURVE* ZeroCurveMakeEmpty(int numElems    /** (I) num zeros */
                           );


int ZeroCurvePopulateFromRates(
            IRDate          baseDate,
            int             numElems,
            const IRDate   *zeroDates,
            const double   *zeroRates,
            T_CURVE* crv);

/* Zero curve (T_CURVE) constructor.                                        */
/* IMPORTANT:  Input zero rates are assumed to be annually compounded.      */
/*             Also, DCC is hard-coded as ACT/365F.                         */
/* Caller is responsible for freeing returned T_CURVE via ZeroCurveFree.    */
T_CURVE* ZeroCurveMakeFromRates(
    IRDate          baseDate,       /** (I) value date & today              */
    int             numElems,       /** (I) num zero dates and rates        */
    const IRDate   *zeroDates,      /** (I) zero dates                      */  
    const double   *zeroRates       /** (I) zero rates (ANNUAL)             */          
    );

/* Answers a copy of the input T_CURVE.                                     */
/* Caller is responsible for freeing returned T_CURVE via ZeroCurveFree.    */
T_CURVE* ZeroCurveCopy(const T_CURVE *src   /** (I) zero curve to copy      */
                       );

/*  Properly dispose of a T_CURVE created via ZeroCurveMakeEmpty,           */
/*  ZeroCurveMakeFromRates, or ZeroCurveCopy.                               */ 
void ZeroCurveFree(T_CURVE *crv     /** (I) zero curve to free              */
                   );


/* Simply answers the base date of the supplied zero curve. */
IRDate GetZeroCurveBaseDate(const T_CURVE* zc);

/* Simply answers the last date of the supplied zero curve. */
IRDate GetZeroCurveLastDate(const T_CURVE* zc);



/*****  GetZeroPrice  ****************************************************/
/**
 *      Calculates a deterministic zero bond price
 *      Returns -999.99 if failure
 *
 *      IF NOT USING IRX CURVES, supports
 *          (EslGetZeroInterpolation == ESL_INTERP_LINEAR) Linear Zero Cpn:
 *          -- linear zero cpn interp
 *          -- flat zero cpn extrapolation on both sides
 *
 *          (EslGetZeroInterpolation() == ESL_INTERP_FLATFWD) Flat Fwd:
 *          -- flat fwd interp and extrapolation
 */


double   GetZeroPrice(
            IRDate                  MatDate  /** (I) Mat date of the zero       */
            ,T_CURVE const*        crv);     /** (I) Zero curve                 */

#ifndef ESL_NEW_CURVE
/* DEPRECATED!!! - Use GetZeroPrice (above) instead.                          */
double ZeroPrice(IRDate         MatDate    /** (I) Mat date of the zero       */
                 ,IRDate        ValueDate  /** (I) Value date                 */
                 ,int           NbZero     /** (I) Number of zeros in curve   */
                 ,IRDate const* ZeroDates  /** (I) maturity dates in zero crv */
                 ,double const* ZeroRates  /** (I) Zero rates                 */
        );

/* DEPRECATED!!!  Use GetZeroPriceRate instead.                          */
int  Get_Zero(double   *OutZeroRate   /** (O) Zero rate                   */
             ,double   *OutZeroPrice  /** (O) Zero price                  */
             ,int       NbZero        /** (I) Number of zeros in the ZC   */
             ,double   *Zero          /** (I) Z rates on an ACT/365 basis */
             ,long     *ZeroDate      /** (I) Zero maturity dates         */
             ,long      CurrentDate   /** (I) Current date                */
             ,long      Maturity      /** (I) Maturity of the zero        */
         );
#endif


/*****  GetZero  *************************************************************/
/**
 *     Calculate zero rate & price for a specific maturity from the zero
 *     curve stored in t_curve style (deterministic).
 *
 *      NOTE - deprecated - use GetZeroPriceRate instead.
 */
int  GetZero(double         *ZeroRate      /** (O) Zero rate                   */
            ,double         *ZeroPrice     /** (O) Zero price                  */
            ,const T_CURVE*  zc            /** (I) Zero curve                  */
            ,IRDate           Maturity      /** (I) Maturity of the zero        */
         );


int  GetZeroPriceRate(
         double   *ZeroRate,  /** (O) Zero rate                   */
         double   *ZeroPrice, /** (O) Zero price                  */
         IRDate     MatDate,   /** (I) Maturity of the zero        */
         T_CURVE const* crv);

/*****  ExtendTCurve  ********************************************************/
/**
 *      Flat extend zero curve between fromDate and toDate dates
 *
 *      FIXME:  Not implemented for IRX curves.  See, however,
 *              irxZeroCurveExtendToToday.
 */
int ExtendTCurve(T_CURVE *tc, IRDate fromDate, IRDate toDate);





/* IMPORTANT:  IRX curves do not support linear interpolation and, by       *
**             default, are built with flat forward interpolation.  One     *
**             cannot simply set the interpolation type for IRX curves.     *
**             This needs to be specifed at the time of curve construction. */ 

/* VERY IMPORTANT:  ALL NEW CODE SHOULD NOT ACCESS THESE FLAGS DIRECTLY!!!  *
**                  INSTEAD, PLEASE USE THE FUNCTIONS BELOW.                */
extern ESL_INTERP_TYPE ZeroInterpTypeFlag;
extern ESL_INTERP_TYPE ZeroInterpTypeFlagStub;

ESL_INTERP_TYPE  EslGetZeroInterpolation();
ESL_INTERP_TYPE  EslGetZeroInterpolationStub();

/* Sets both interpolation type and stub */
void EslSetZeroInterpolation(ESL_INTERP_TYPE t);

/* Sets just interpolation stub */
void EslSetZeroInterpolationStub(ESL_INTERP_TYPE t);

/* Convenience functions. */
void EslSetZeroInterpolationLinear();
void EslSetZeroInterpolationFlatFwd();


#ifdef  __cplusplus
}
#endif

#endif /* ESL_ZEROS_DOT_H */
