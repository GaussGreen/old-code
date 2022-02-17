/*
***************************************************************************
** HEADER FILE: recovery.h
**
** Defines the recovery curve - this is a time dependent curve of expected
** recoveries.
**
** $Header$
***************************************************************************
*/

#ifndef CX_RECOVERY_H
#define CX_RECOVERY_H

#include "cx.h"  /* basic data types */

/**
***************************************************************************
** Constructs a recovery curve from a single recovery rate.
***************************************************************************
*/
CxTRecoveryCurve* CxRecoveryCurveMakeFromRecoveryRate
(double recoveryRate);

/**
***************************************************************************
** Ensure that the dates in the recovery curve are in the correct order,
** plus various other validations.
***************************************************************************
*/
int CxRecoveryCurveValidate
(CxTRecoveryCurve* curve);

/**
***************************************************************************
** Interpolates a value from a curve.
***************************************************************************
*/
int CxRecoveryCurveInterp
(const CxTRecoveryCurve* curve,        /* (I) Loss curve */
 TDate                   date,         /* (I) Interpolation date */
 double                 *interpValue); /* (O) Interpolated value */

/**
***************************************************************************
** Computes an unweighted average recovery between two dates.
***************************************************************************
*/
int CxRecoveryCurveAverage
(const CxTRecoveryCurve *curve,         /* (I) Recovery curve */
 TDate                   startDate,     /* (I) Start date */
 TDate                   endDate,       /* (I) End date */
 double                 *averageValue); /* (O) Average value */

#endif /* CX_LOSSCURVE_H */

