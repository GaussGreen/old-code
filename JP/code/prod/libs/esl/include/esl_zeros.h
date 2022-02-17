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


void ZeroCurveInit(T_CURVE *crv);

/* Simply answers the base date of the supplied zero curve. */
IRDate GetZeroCurveBaseDate(const T_CURVE* zc);

/* Simply answers the last date of the supplied zero curve. */
IRDate GetZeroCurveLastDate(const T_CURVE* zc);

int IsZeroCurveEmpty(const T_CURVE* zc);

int GetZeroCurveNumPoints(const T_CURVE* zc);


int printZeroRates(FILE* fp, const T_CURVE* zc, const char* formatstr);

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

double ZeroPriceBase(IRDate    MatDate     /** (I) Mat date of the zero       */
                 ,IRDate        ValueDate  /** (I) Value date                 */
                 ,int           NbZero     /** (I) Number of zeros in curve   */
                 ,IRDate const* ZeroDates  /** (I) maturity dates in zero crv */
                 ,double const* ZeroRates  /** (I) Zero rates                 */
                 ,ESL_INTERP_TYPE InterpType 
                 );

#endif


/*****  GetZeroPriceRate  *************************************************************/
/**
 *     Calculate zero rate & price for a specific maturity from the zero
 *     curve stored in t_curve style (deterministic).
 *
 */
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


ESL_INTERP_TYPE  EslGetZeroInterpolation();
/* Sets both interpolation type and stub */
void EslSetZeroInterpolation(ESL_INTERP_TYPE t);
/* Convenience functions. */
void EslSetZeroInterpolationLinear();
void EslSetZeroInterpolationFlatFwd();


#ifdef  __cplusplus
}
#endif

#endif /* ESL_ZEROS_DOT_H */
