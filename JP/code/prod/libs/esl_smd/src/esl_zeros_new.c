/****************************************************************************/
/*      Zero price and zero bank.                                           */
/****************************************************************************/
/*      ZEROS.c                                                             */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "esl_zeros.h"
#include "esl_date.h"
#include "esl_util.h"
#include "esl_error.h"

#include "irx/zerocurve.h"


/**
 *      Calculates a deterministic zero bond price
 *      Returns -999.99 if failure
 *
 *      Supports
 *      -- smooth-forward interp
 *      -- flat-forward interp
 *
 */

double   GetZeroPrice(
            long                  MatDate  /** (I) Mat date of the zero       */
            ,T_CURVE const* crv)
{
    double Price;

    if (irxZeroPrice(crv, MatDate, &Price) != SUCCESS)
        Price = -999.99;

    if (Price < 0.0)
        DR_Error("ZeroPrice: failed.");
    return (Price);
}



/**
 *     Calculate zero rate & price for a specific maturity from the zero
 *	   curve stored in t_curve style (deterministic).
 *
 */
int  GetZeroPriceRate(
         double   *OutZeroRate,  /** (O) Zero rate	                 */
         double   *OutZeroPrice, /** (O) Zero price                  */
         long      MatDate,      /** (I) Maturity of the zero        */
         T_CURVE const* crv)

{
    double t;

    t = Daysact (crv->baseDate, MatDate) / 365.0;

    if (irxZeroPrice(crv, MatDate, OutZeroPrice) != SUCCESS)
        return (FAILURE);

    *OutZeroRate = pow(*OutZeroPrice, -1.0/t ) - 1.0;

    return (SUCCESS);
}

