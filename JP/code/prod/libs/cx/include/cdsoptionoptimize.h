/*
*****************************************************************************
** CDS option optimization functions.
** 
** Computes the optimal values of the multi-Q parameters to fit the market
** data.
*****************************************************************************
*/

#ifndef CX_CDSOPTIONOPTIMIZE_H
#define CX_CDSOPTIONOPTIMIZE_H

#include "cdsoption.h"

/**
***************************************************************************
** Multi-Q parameter optimization initialization.
**
** Initializes the state object with the initial guess of the ATM vol and
** Q-distribution.
***************************************************************************
*/
CrxTCdsOptionQOptimization* CrxCdsOptionQOptimizationInit(
/** Initial guess for the Q-distribution */
    CrxTQDist*      qdist,
/** Initial guess for the ATM volatility */
    double          vol);

/**
***************************************************************************
** Multi-Q parameter optimization.
**
** Designed so that we can pass the results back as an input to allow the
** function to be called in a loop for gradual improvements.
***************************************************************************
*/
CrxTCdsOptionQOptimization* CrxCdsOptionQOptimization(
/** Initial state of the optimization */
    CrxTCdsOptionQOptimization* initState,
/** Base deal. Describes the main properties of the deal under consideration.
    In addition, we allow the strike and optionType to be varied. */
    CrxTCdsOption*              baseDeal,
/** Model parameters controlling the optimization */
    CrxTCdsOptionQOptModel*     model,
/** Today - volatility starts today. */
    TDate                       today,
/** Number of options to be used in the optimization */
    int                         numOptions,
/** Strike spreads of the options. Measured in decimals e.g. 100bp is
    represented as 0.01 */
    double*                     strikes,
/** Option types of the options, i.e. call or put */
    long*                       optionTypes,
/** Weight to be applied to this particular option. Typical value = 1 */
    double*                     weights,
/** Observed Black-Scholes log-normal volatility for this strike. */
    double*                     strikeVols,
/** Interest rate discount curve. */
    TCurve*                     discCurve,
/** Credit spread curve. */
    CxTCreditCurve*             sprdCurve,
/** Credit recovery rate in case of default. */
    double                      recoveryRate,
/** Name of logfile. Can be NULL. The previous contents of the logfile
    will be overwritten. */
    char*                       logfilename
);

#endif
