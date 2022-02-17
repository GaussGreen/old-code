/*
***************************************************************************
** HEADER FILE: mqdata.h
**
** Data manipulation functions for the MQDATA structure.
***************************************************************************
*/

#ifndef _CRX_MQDATA_H
#define _CRX_MQDATA_H


#include <crxflow/include/crxdata.h>

#ifdef __cplusplus
extern "C"
{
#endif

#include <crxmultiq/include/crmultiq.h>
            
/**
***************************************************************************
** Constructor for MQDATA
***************************************************************************
*/
MQDATA* CrxMqdataMake(
/** Expiration time */
double          optExpy,
/** Forward rate */
double          fwdRate,
/** Market at the money volatility */
double          sigATM,
/** Target at the money option price */
double          optATM,
/** Target at the money option type */
long            optType,
/** Defines array size for kL, dL, xL, bL, qL.
    Number of left intervals */
long            nbQL,
/** Array of size nbQL.
    Left hand values of k */
double*         kL,
/** Array of size nbQL.
    Left hand values of d */
double*         dL,
/** Array of size nbQL+1.
    Left hand values of x */
double*         xL,
/** Array of size nbQL.
    Left hand values of b */
double*         bL,
/** Array of size nbQL.
    Left hand values of q */
double*         qL,
/** Defines array size for kR, dR, xR, bR, qR.
    Number of right intervals */
long            nbQR,
/** Array of size nbQR.
    Right hand values of k */
double*         kR,
/** Array of size nbQR.
    Right hand values of d */
double*         dR,
/** Array of size nbQR+1.
    Right hand values of x */
double*         xR,
/** Array of size nbQR.
    Right hand values of b */
double*         bR,
/** Array of size nbQR.
    Right hand values of q */
double*         qR,
/** Volatility of forward distribution */
double          sigMQ,
/** Calibrated at the money option price */
double          optMQ,
/** 1/E[H(x)] */
double          C,
/** F(x,C,K) = K(CH(x)-1) + 1 */
double          K,
long            calibType,
long            sSteps,
double          sDelta,
/** Calibration tolerance for forward */
double          fwdTol,
/** Calibration tolerance for at the money volatility */
double          atmTol,
long            calcFwd,
double          muMQ
);

/**
***************************************************************************
** Memory allocator for MQDATA
***************************************************************************
*/
MQDATA* CrxMqdataMakeEmpty(
/** Defines array size for kL, dL, xL, bL, qL.
    Number of left intervals */
long            nbQL,
/** Defines array size for kR, dR, xR, bR, qR.
    Number of right intervals */
long            nbQR
);

/**
***************************************************************************
** Copy constructor for MQDATA
***************************************************************************
*/
MQDATA* CrxMqdataCopy(MQDATA* src);

/**
***************************************************************************
** Destructor for MQDATA
***************************************************************************
*/
void CrxMqdataFree(MQDATA *p);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _CRX_MQDATA_H */
