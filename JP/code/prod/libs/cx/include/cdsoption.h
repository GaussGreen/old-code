/*
*****************************************************************************
** CDS option pricer analytics.
** 
** This uses simple option pricing using the forward spread of the option,
** with support for the multi-Q option pricer to incorporate smile parameters.
**
** The model supports both single name CDS options and index CDS options.
*****************************************************************************
*/

#ifndef CX_CDSOPTION_H
#define CX_CDSOPTION_H

#include <crxflow/include/crxdata.h>
#include "cx.h"

/*f
***************************************************************************
** CDS option calculator
***************************************************************************
*/
CrxTCdsOptionCalc* CrxCdsOptionCalc(
CrxTCdsOption*  deal,             /* (I) */
long            distType,         /* (I) */
TDate           today,            /* (I) */
double          vol,              /* (I) */
CrxTQDist*      qdist,            /* (I) */
TCurve*         discCurve,        /* (I) */
CxTCreditCurve* sprdCurve,        /* (I) */
double          recoveryRate      /* (I) */
);

/*f
***************************************************************************
** CDS option implied volatility calculator
***************************************************************************
*/
CrxTCdsOptionCalc* CrxCdsOptionVolCalc(
CrxTCdsOption*  deal,             /* (I) */
long            distType,         /* (I) */
TDate           today,            /* (I) */
double          price,            /* (I) */
CrxTQDist*      qdist,            /* (I) */
TCurve*         discCurve,        /* (I) */
CxTCreditCurve* sprdCurve,        /* (I) */
double          recoveryRate      /* (I) */
);

#endif
