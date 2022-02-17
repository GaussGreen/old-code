/*
***************************************************************************
** FILE NAME: creditcurve.h
**
** Credit curve functions.
***************************************************************************
*/

#ifndef CX_CREDIT_CURVE_H
#define CX_CREDIT_CURVE_H

#include "cx.h"
#include <alib/zcurve.h>

/**
***************************************************************************
** Validates (and potentially changes) a credit curve as part of its
** construction.
***************************************************************************
*/
int CxCreditCurveValidate
(CxTCreditCurve *cc);

/**
***************************************************************************
** Constructs a credit curve. Easier to use this in a spreadsheet
** environment than first constructing a TCurve.
***************************************************************************
*/
CxTCreditCurve* CxCreditCurveNew
(CxTCreditCurveType type,
 TDate             baseDate,
 int               numDates,
 TDate            *dates,
 double           *spreads,
 long              rateType,
 CxTDayCountConv   dayCountConv,
 TDateInterval    *ivl);

/**
***************************************************************************
** Routine for creating a smooth curve out of a flat forward curve.
***************************************************************************
*/
CxTCreditCurve* CxCreditCurveMakeSmoother
(CxTCreditCurve *cc,
 TDateList     *newDates);

/**
***************************************************************************
** Routine for converting the compounding basis of a credit curve.
**
** The curve is amended in place.
***************************************************************************
*/
int CxCreditCurveConvertRateType
(CxTCreditCurve *cc, 
 long            newRateType);

/**
***************************************************************************
** Returns the hazard rate for a particular date from a credit curve.
**
** Since survival probabilities in the credit curve are valued at the end
** of day, this calculation involves computing the survival probability from
** one day before and the current date, computing the one day hazard,
** and then converting this into an annualized exponential rate.
**
** For dates on or before the start date of the curve, the hazard rate is
** zero.
***************************************************************************
*/
int CxCreditCurveHazardRate
(CxTCreditCurve *creditCurve,
 TDate           date,
 double         *hazardRate);

/**
***************************************************************************
** Returns the survival probability computed at the end of the day.
** This is computed relative to the base date of the credit curve.
***************************************************************************
*/
int CxCreditCurveSurvivalProb
(CxTCreditCurve *creditCurve,
 TDate           date,
 double         *survivalProb);

/**
***************************************************************************
** Returns the conditional survival probability computed at the end of the
** day for endDate, conditional on survival to the end of the day for
** startDate.
***************************************************************************
*/
int CxCreditCurveConditionalSurvivalProb
(CxTCreditCurve *creditCurve,
 TDate           startDate,
 TDate           endDate,
 double         *survivalProb);


/**
***************************************************************************
** Slides the credit curve to a new date in the future.
***************************************************************************
*/
CxTCreditCurve* CxCreditCurveSlide
(CxTCreditCurve *creditCurve,
 TDate           slideDate);


#endif
 
 
 
