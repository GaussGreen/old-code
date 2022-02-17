/*===========================================================
  Name    : hw_vfdk_LDHD_lattice.cpp
  Owner   : DK
  Created : DK  
  Comment : Don't touch if your name is not DK !
  Time    : Last update 13:36 06/28/2004
=============================================================*/
#include <iostream>
#include <strstream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifndef WIN32
#include <sys/param.h>
#endif 
#define ACC 1.0e-04
#define JMAX 5000
#include "hw_vfdk_LDHD_lattice.h"
#include "Swaption.h"
#include "dk_utils.h"
#include "hw_vfdk_analytics.h"
#include "DKMaille.h"
#include "DKMaille2D.h"
#include "dk_smoothing.h"

#include "linalg.h"

DKMaille<double> g_dCurveDates;
DKMaille<double> g_dCurve;
DKMaille<double> g_dForeignCurveDates;
DKMaille<double> g_dForeignCurve;

/// JMP : Extraction de données du lattice sous Xl
ARM_Matrix DomVolGlobalExtraction;
ARM_Matrix ForVolGlobalExtraction;
ARM_Matrix FxVolGlobalExtraction;
ARM_Matrix VolAnalyticsGlobalExtraction;
ARM_Matrix LatticeVolsGlobalExtraction;
ARM_Vector NbNodesPerSliceGlobalExtraction;
ARM_Matrix AnalyticFxOptionGlobalExtraction;
int FxOptIdxGlobal=0;

double g_dOptimal;
double g_dTimeBoost;
unsigned int g_dSpotDate;
unsigned int g_uiNumFactors;
unsigned int g_uiProductModelCode;
unsigned int g_uiBranching;
unsigned int g_ppy;
unsigned int g_ppy1;
double g_dTreeCutOffDate;
double g_dTreeFirstExpiry;
double g_dOptionExpiry;
double g_dStrike;
double g_dType;
DKMaille<double> g_dDiscountCurveBase;
DKMaille<double> g_dDiscountCurveBase_FixedIncomeExotics;
DKMaille<double> g_dDiscountCurveForeign;
DKMaille<double> g_dZCCurveBase;
DKMaille<double> g_dZCCurveForeign;
DKMaille<double> g_dDiscountCurveBaseDates;
DKMaille<double> g_dDiscountCurveBaseDates_FixedIncomeExotics;
DKMaille<double> g_dDiscountCurveForeignDates;
DKMaille<double> g_dDiscountCurveBaseNoBasis;
DKMaille<double> g_dDiscountCurveBaseNoBasis_FixedIncomeExotics;
DKMaille<double> g_dDiscountCurveForeignNoBasis;
DKMaille<double> g_dNoticeDates_FixedIncomeExotics;
DKMaille<double> g_dModelParameters_FixedIncomeExotics;
double g_RHO;
double g_MR1;
double g_MR2;
double g_S1;
double g_S2;
double g_MR3;
double g_S3;
double g_RHO1;
double g_RHO2;
double g_dLimitXX1;
double g_dLimitXX2;
double g_dLimitXX3;
double g_GridScaling1;
double g_GridScaling2;
double g_GridScaling3;
double g_dSurvivalCalculation;
double g_dRedemption;
double g_dRedemptionOption;
double g_dRedemptionStrike;
double g_dRedemptionDate;
double g_dRedemptionNoticeDate;
double g_dNumSlicesBeforeFirstNotice;
double g_dQParameter_Base;
double g_dQParameter_Foreign;
DKMaille<double> g_dEquilibriumForwardBase;
DKMaille<double> g_dEquilibriumForwardForeign;

DKMaille<double> g_dEquilibriumForwardShortBase;
DKMaille<double> g_dEquilibriumForwardShortForeign;


DKMaille<double> g_dEquilibriumForwardBase_Derivative;
DKMaille<double> g_dEquilibriumForwardForeign_Derivative;
unsigned int g_uiI1Limit;
unsigned int g_uiI2Limit;
unsigned int g_uiI3Limit;
DKMaille<double> g_dLimitY1;
DKMaille<double> g_dLimitY2;
DKMaille<double> g_dLimitY3;
double g_dCorrelation;
double g_dCorrelation_1;
double g_dCorrelation_2;
double g_dCorrelationInitial;
double g_dRelativeFactor;
double g_dNormalisation_11;
double g_dNormalisation_12;
double g_dNormalisation_21;
double g_dNormalisation_22;
double g_dNormalisation_Inverse_11;
double g_dNormalisation_Inverse_12;
double g_dNormalisation_Inverse_21;
double g_dNormalisation_Inverse_22;
double g_d3Normalisation_11;
double g_d3Normalisation_12;
double g_d3Normalisation_21;
double g_d3Normalisation_22;
double g_d3Normalisation_13;
double g_d3Normalisation_23;
double g_d3Normalisation_31;
double g_d3Normalisation_32;
double g_d3Normalisation_33;
double g_d3Normalisation_Inverse_11;
double g_d3Normalisation_Inverse_12;
double g_d3Normalisation_Inverse_21;
double g_d3Normalisation_Inverse_22;
double g_d3Normalisation_Inverse_13;
double g_d3Normalisation_Inverse_23;
double g_d3Normalisation_Inverse_31;
double g_d3Normalisation_Inverse_32;
double g_d3Normalisation_Inverse_33;
double g_dADPLimit;
DKMaille<double> g_dSpecialDates;
DKMaille<double> g_dDates;
DKMaille<double> g_LATTICE_DATE;
DKMaille<double> g_S1_TD;
DKMaille<double> g_S1_DERIVATIVE_TD;
DKMaille<double> g_S2_TD;
DKMaille<double> g_S2_DERIVATIVE_TD;
DKMaille<double> g_S3_TD;
DKMaille<double> g_S3_DERIVATIVE_TD;
DKMaille<double> g_dStdDevBaseX;
DKMaille<double> g_dStdDevBaseZ;
DKMaille<double> g_dStdDevForeignX;
DKMaille<double> g_dStdDevForeignZ;
DKMaille<double> g_dSpotFXVolDatesTD;
DKMaille<double> g_dSpotFXVolTD;
DKMaille<double> g_LATTICE_NORM;
DKMaille<double> g_LATTICE_NORM_FOREIGN;
DKMaille<double> g_LATTICE_DECAY;
DKMaille<double> g_LATTICE_DECAY_FOREIGN;
DKMaille<double> g_LATTICE_PERIOD_TIME_SPAN;
DKMaille<double> g_LATTICE_ZCB_PRICE;
DKMaille<double> g_LATTICE_ZCB_PRICE_FOREIGN;
DKMaille<double> g_LATTICE_RATE_STEP1;
DKMaille<double> g_LATTICE_RATE_STEP2;
DKMaille<double> g_LATTICE_RATE_STEP3;
DKMaille<long> g_LATTICE_MAX_STATE1;
DKMaille<long> g_LATTICE_MAX_STATE2;
DKMaille<long> g_LATTICE_MAX_STATE3;
DKMaille<long> g_LATTICE_MIN_STATE1;
DKMaille<long> g_LATTICE_MIN_STATE2;
DKMaille<long> g_LATTICE_MIN_STATE3;
DKMaille<unsigned int> g_LATTICE_NUMBER_OF_STATES1;
DKMaille<unsigned int> g_LATTICE_NUMBER_OF_STATES2;
DKMaille<unsigned int> g_LATTICE_NUMBER_OF_STATES3;
DKMaille< DKMaille<double> > g_LATTICE_X1_RATE;
DKMaille< DKMaille<double> > g_LATTICE_X2_RATE;
DKMaille< DKMaille<double> > g_LATTICE_X3_RATE;
DKMaille< DKMaille<double> > g_dStdDevCutOff;
DKMaille< DKMaille<double> > g_dUpper;
DKMaille< DKMaille<double> > g_dLower;
DKMaille< DKMaille<double> > g_dNextSliceLimitUp;
DKMaille< DKMaille<double> > g_dCurrentSliceLimitUp;
DKMaille< DKMaille<double> > g_dNextSliceLimitDown;
DKMaille< DKMaille<double> > g_dCurrentSliceLimitDown;
DKMaille< DKMaille<double> > g_dCurrentCenter;
DKMaille< DKMaille<double> > g_dNextCenter;
DKMaille< DKMaille<int> > g_iLimit;
DKMaille< DKMaille<int> > g_iCurrentShift;
DKMaille< DKMaille<int> > g_iNextShift;
DKMaille< DKMaille<double> > g_LATTICE_G_RATE;
DKMaille< DKMaille<double> > g_LATTICE_SHORT_RATE;
DKMaille< DKMaille<double> > g_LATTICE_ARROW_DEBREU;
DKMaille< DKMaille<double> > g_LATTICE_ARROW_DEBREU_PROB;
DKMaille< DKMaille<double> > g_LATTICE_G_RATE_FOREIGN;
DKMaille< DKMaille<double> > g_LATTICE_SHORT_RATE_FOREIGN;
DKMaille< DKMaille<double> > g_LATTICE_G_RATE_FX;
DKMaille< DKMaille<double> > g_LATTICE_SHORT_RATE_FX;
DKMaille< DKMaille<double> > g_LATTICE_SHORT_RATE_FX_PRIME;
DKMaille< DKMaille<bool> > g_LATTICE_NODE_EXISTENCE;
DKMaille< DKMaille<int> > g_LATTICE_TARGET_1;
DKMaille< DKMaille<int> > g_LATTICE_TARGET_2;
DKMaille< DKMaille<int> > g_LATTICE_TARGET_3;
DKMaille< DKMaille<double> > g_LATTICE_UP_1;
DKMaille< DKMaille<double> > g_LATTICE_UP_2;
DKMaille< DKMaille<double> > g_LATTICE_UP_3;
DKMaille< DKMaille<double> > g_LATTICE_MID_1;
DKMaille< DKMaille<double> > g_LATTICE_MID_2;
DKMaille< DKMaille<double> > g_LATTICE_MID_3;
DKMaille< DKMaille<double> > g_LATTICE_DOWN_1;
DKMaille< DKMaille<double> > g_LATTICE_DOWN_2;
DKMaille< DKMaille<double> > g_LATTICE_DOWN_3;
// DKMaille2D< DKMaille<double> > g_LATTICE_GREEN;
DKMaille2D< DKMaille<unsigned int> > g_LATTICE_CONNECTION;
DKMaille2D< DKMaille<unsigned int> > g_LATTICE_GREEN;
unsigned int g_dConvert=100000000;
double g_dIsSwaptionCalibrationWithBasis;
DKMaille<double> g_dShiftShortGRate;
DKMaille<double> g_dShiftShortGRateBase;
DKMaille<double> g_dShiftShortGRateForeign;
DKMaille<double> g_dFromForeignDriftToCash;
DKMaille< DKMaille<double> > g_d1DGreenBase;
DKMaille< DKMaille<double> > g_d1DGreenForeign;
DKMaille< DKMaille<double> > g_d1DShortRateBase;
DKMaille< DKMaille<double> > g_d1DShortRateForeign;
DKMaille<double> g_dNoticeDates;
DKMaille2D<double> g_sigma;
DKMaille2D<double> g_meanRevRate;
// DKMaille2D<double> g_SigmaOriginal;
bool g_bStringModel;
unsigned int g_uiModelIdentifier;
unsigned int g_uiNumSteps;
unsigned int g_uiMinimum;
unsigned int g_uiNumSlices;
double g_dSpotFX;
double g_dSpotFXTilde;
double g_dDeltaFlag;
unsigned int g_uiDeltaTick;
double g_dLastDate;
double g_dFXForwardAtMaturity;
unsigned int g_uiSmoothing;
double g_dLocalFXVol;
double g_dLocalSwaptionVol;
DKMaille<unsigned int> g_uiVolIsCalculated;
DKMaille<double> g_dStrip;
DKMaille<double> g_dStripJulian;
DKMaille<double> g_dStripDomesticStdDev;
DKMaille<double> g_dStripForeignStdDev;
DKMaille<double> g_dStripSpotFXVol;
DKMaille<double> g_dDiffusion_11;
DKMaille<double> g_dDiffusion_12;
DKMaille<double> g_dDiffusion_13;
DKMaille<double> g_dDiffusion_21;
DKMaille<double> g_dDiffusion_22;
DKMaille<double> g_dDiffusion_23;
DKMaille<double> g_dDiffusion_31;
DKMaille<double> g_dDiffusion_32;
DKMaille<double> g_dDiffusion_33;
DKMaille<double> g_LATTICE_DATE_FULL;
DKMaille<double> g_S1_TD_FULL;
DKMaille<double> g_S2_TD_FULL;
DKMaille<double> g_S3_TD_FULL;
double g_dIsFundingLegStochastic;
unsigned int g_uiIsDK;
unsigned int g_uiIsSP;

void deleteAll()
{

    g_dDiscountCurveBase.clear();
    g_dDiscountCurveForeign.clear();
    g_dZCCurveBase.clear();
    g_dZCCurveForeign.clear();
    g_dDiscountCurveBaseDates.clear();
    g_dDiscountCurveForeignDates.clear();
    g_dLimitY1.clear();
    g_dLimitY2.clear();
    g_dLimitY3.clear();
    g_dSpecialDates.clear();
    g_dDates.clear();
    g_LATTICE_DATE.clear();
    g_S1_TD.clear();
    g_S1_DERIVATIVE_TD.clear();
    g_S2_TD.clear();
    g_S2_DERIVATIVE_TD.clear();
    g_S3_TD.clear();
    g_S3_DERIVATIVE_TD.clear();
    g_dStdDevBaseX.clear();
    g_dStdDevBaseZ.clear();
    g_dStdDevForeignX.clear();
    g_dStdDevForeignZ.clear();
    g_dSpotFXVolDatesTD.clear();
    g_dSpotFXVolTD.clear();
    g_LATTICE_NORM.clear();
    g_LATTICE_NORM_FOREIGN.clear();
    g_LATTICE_DECAY.clear();
    g_LATTICE_DECAY_FOREIGN.clear();
    g_LATTICE_PERIOD_TIME_SPAN.clear();
    g_LATTICE_ZCB_PRICE.clear();
    g_LATTICE_ZCB_PRICE_FOREIGN.clear();
    g_LATTICE_RATE_STEP1.clear();
    g_LATTICE_RATE_STEP2.clear();
    g_LATTICE_RATE_STEP3.clear();
    g_LATTICE_MAX_STATE1.clear();
    g_LATTICE_MAX_STATE2.clear();
    g_LATTICE_MAX_STATE3.clear();
    g_LATTICE_MIN_STATE1.clear();
    g_LATTICE_MIN_STATE2.clear();
    g_LATTICE_MIN_STATE3.clear();
    g_LATTICE_NUMBER_OF_STATES1.clear();
    g_LATTICE_NUMBER_OF_STATES2.clear();
    g_LATTICE_NUMBER_OF_STATES3.clear();
    g_LATTICE_X1_RATE.clear();
    g_LATTICE_X2_RATE.clear();
    g_LATTICE_X3_RATE.clear();
    g_dStdDevCutOff.clear();
    g_dUpper.clear();
    g_dLower.clear();
    g_dNextSliceLimitUp.clear();
    g_dCurrentSliceLimitUp.clear();
    g_dNextSliceLimitDown.clear();
    g_dCurrentSliceLimitDown.clear();
    g_dCurrentCenter.clear();
    g_dNextCenter.clear();
    g_iLimit.clear();
    g_iCurrentShift.clear();
    g_iNextShift.clear();
    g_LATTICE_G_RATE.clear();
    g_LATTICE_SHORT_RATE.clear();
    g_LATTICE_ARROW_DEBREU.clear();
    g_LATTICE_ARROW_DEBREU_PROB.clear();
    g_LATTICE_G_RATE_FOREIGN.clear();
    g_LATTICE_SHORT_RATE_FOREIGN.clear();
    g_LATTICE_G_RATE_FX.clear();
    g_LATTICE_SHORT_RATE_FX.clear();
    g_LATTICE_SHORT_RATE_FX_PRIME.clear();
    g_LATTICE_NODE_EXISTENCE.clear();
    g_LATTICE_TARGET_1.clear();
    g_LATTICE_TARGET_2.clear();
    g_LATTICE_TARGET_3.clear();
    g_LATTICE_UP_1.clear();
    g_LATTICE_UP_2.clear();
    g_LATTICE_UP_3.clear();
    g_LATTICE_MID_1.clear();
    g_LATTICE_MID_2.clear();
    g_LATTICE_MID_3.clear();
    g_LATTICE_DOWN_1.clear();
    g_LATTICE_DOWN_2.clear();
    g_LATTICE_DOWN_3.clear();
    // g_LATTICE_CONNECTION.clear();
    // g_LATTICE_GREEN.clear();
    g_dShiftShortGRate.clear();
    g_dShiftShortGRateBase.clear();
    g_dShiftShortGRateForeign.clear();
    g_dFromForeignDriftToCash.clear();
    g_d1DGreenBase.clear();
    g_d1DGreenForeign.clear();
    g_d1DShortRateBase.clear();
    g_d1DShortRateForeign.clear();
    g_dNoticeDates.clear();
    // g_sigma.clear();
    // g_meanRevRate.clear();
    // g_SigmaOriginal.clear();
    g_uiVolIsCalculated.clear();
    g_dDiffusion_11.clear();
    g_dDiffusion_12.clear();
    g_dDiffusion_13.clear();
    g_dDiffusion_21.clear();
    g_dDiffusion_22.clear();
    g_dDiffusion_23.clear();
    g_dDiffusion_31.clear();
    g_dDiffusion_32.clear();
    g_dDiffusion_33.clear();

    g_LATTICE_DATE_FULL.clear();
    g_S1_TD_FULL.clear();
    g_S2_TD_FULL.clear();
    g_S3_TD_FULL.clear();

}


double hwmr(unsigned uiFactorIndex, unsigned uiSliceIndex)
{
    return g_meanRevRate.at(uiFactorIndex-1, uiSliceIndex);
} // hwmr(...)

unsigned int NumSlicesBeforeFirstNotice2F(double dFirstNoticeDate)
{
    unsigned int uiSlices=0;
    uiSlices=(unsigned int)(g_dNumSlicesBeforeFirstNotice*dFirstNoticeDate);
    return uiSlices;
}

unsigned int NumSlicesBeforeFirstNotice3F(double dFirstNoticeDate)
{
    unsigned int uiSlices=0;
    uiSlices=(unsigned int)(g_dNumSlicesBeforeFirstNotice*dFirstNoticeDate);
    if(uiSlices<4) uiSlices=4;
    if(uiSlices>20) uiSlices=(unsigned int)2.*dFirstNoticeDate;
    return uiSlices;
}

unsigned int NumSlicesBeforeFirstNotice1F(double dFirstNoticeDate)
{
    unsigned int uiSlices=0;
    uiSlices=(unsigned int)(g_dNumSlicesBeforeFirstNotice*dFirstNoticeDate);
    return uiSlices;
}

void CalcSimpleLatticeDates(DKMaille<double> dLatticeGeometryData,DKMaille<double> &dates)
{
    dates.resize(dLatticeGeometryData.entries()-1);
    for(unsigned int ui=1;ui<dLatticeGeometryData.entries();ui++)
        dates.at(ui-1)=(dLatticeGeometryData.at(ui)-dLatticeGeometryData.at(1))/365.;
}


void NewTimeSliceGenerator(DKMaille<double> &dates, unsigned int uiTotal)
{
    double datesArray[2000];
    unsigned uiNumDates = g_dSpecialDates.entries();
    bool    bCramToFirstExpiry = false;
    double dLastDate = g_dSpecialDates.at(uiNumDates-1);
    double dStartDate=g_dSpecialDates.at(0);
    g_uiNumSteps=0;

    bCramToFirstExpiry = true;


    unsigned uiStepsToFirstExpiry = 0;
    unsigned uiStepsAfterFirstExpiry = 0;
    unsigned uiStepsInTree=0;

    if (bCramToFirstExpiry)
    {
        if((g_dTreeCutOffDate-g_dSpecialDates.at(0)>12.0)
                &&(g_dTreeCutOffDate-g_dSpecialDates.at(0)<=15.0)
                &&(g_ppy>6))
            uiStepsAfterFirstExpiry = (unsigned)(6.*(dLastDate-g_dSpecialDates.at(1)));
        else if((g_dTreeCutOffDate-g_dSpecialDates.at(0)>15.0)
                &&(g_dTreeCutOffDate-g_dSpecialDates.at(0)<=20.0)
                &&(g_ppy>4))
            uiStepsAfterFirstExpiry = (unsigned)(4.*(dLastDate-g_dSpecialDates.at(1)));
        else if((g_dTreeCutOffDate-g_dSpecialDates.at(0)>20.0)
                &&(g_ppy>4))
            uiStepsAfterFirstExpiry = (unsigned)(4.*(dLastDate-g_dSpecialDates.at(1)));
        else uiStepsAfterFirstExpiry = (unsigned) (g_ppy * (dLastDate-g_dSpecialDates.at(1)));
        uiStepsToFirstExpiry = (unsigned)(sqrt((g_dSpecialDates.at(1) - dStartDate)*365.25)*2.);
    }

    switch (g_uiNumFactors)
    {
    case 1:
        uiStepsToFirstExpiry = NumSlicesBeforeFirstNotice1F(g_dSpecialDates.at(1));
        break;
    case 2:
        uiStepsToFirstExpiry = NumSlicesBeforeFirstNotice2F(g_dSpecialDates.at(1));
        break;
    case 3:
        uiStepsToFirstExpiry = NumSlicesBeforeFirstNotice3F(g_dSpecialDates.at(1));
        break;
    default:
        throw("ERROR: 1 or 2 factor only implemented so far");
        break;
    } // switch

    for (unsigned int ii=0; ii < g_uiNumSteps; ii++)
    {
        dates[ii] = datesArray[ii];
    }
} // NewTimeSliceGenerator(...)

double LatticeDiscount(unsigned int ui,unsigned int uk)
{
    if(g_bStringModel&&g_uiProductModelCode!=2.&&g_uiProductModelCode!=3.)
        return g_LATTICE_NORM.at(ui)*exp(-g_LATTICE_DECAY.at(ui)*g_LATTICE_SHORT_RATE.at(ui).at(uk));
    else
        return exp(-g_LATTICE_SHORT_RATE.at(ui).at(uk)*g_LATTICE_PERIOD_TIME_SPAN.at(ui));
} // LatticeDiscount

double LatticePrimaryDiscount(unsigned int ui,double dShortRate)
{
    // Function to be called only for standard non-string model
    return exp(-dShortRate*g_LATTICE_PERIOD_TIME_SPAN.at(ui));
} // LatticePrimaryDiscount

double LatticeDiscountForeign(unsigned int ui,unsigned int uk)
{
    if(g_bStringModel&&g_uiProductModelCode!=2.&&g_uiProductModelCode!=3.)
        return g_LATTICE_NORM_FOREIGN.at(ui)*exp(-g_LATTICE_DECAY_FOREIGN.at(ui)*g_LATTICE_SHORT_RATE_FOREIGN.at(ui).at(uk));
    else
        return exp(-g_LATTICE_SHORT_RATE_FOREIGN.at(ui).at(uk)*g_LATTICE_PERIOD_TIME_SPAN.at(ui));
} // LatticeDiscountForeign

double LatticePrimaryDiscountForeign(unsigned int ui,double dShortRateForeign)
{
    // Function to be called only for standard non-string model
    return exp(-dShortRateForeign*g_LATTICE_PERIOD_TIME_SPAN.at(ui));
} // LatticePrimaryDiscountForeign

void CalcLatticeDates(DKMaille<double> &dates, unsigned int uiMinimum)
{
    double datesArray[2000];
    double dStepSize;
    double dNewStepSize;
    double dStepGap;
    double dPreStepSize;
    unsigned uiNumDates = g_dSpecialDates.entries();

    double dLastDate = g_dSpecialDates.at(uiNumDates-1);
    double dStartDate=g_dSpecialDates.at(0);
    g_uiNumSteps=0;
    unsigned uiStepsToFirstExpiry = 0;
    unsigned uiStepsAfterFirstExpiry = 0;
    unsigned uiStepsInTree=0;

    unsigned uiStepsAfterFirstExpiryAndBeforeOptimal = (unsigned) (g_ppy1 * (g_dOptimal-g_dSpecialDates.at(1))+0.5);

    unsigned uiStepsAfterOptimal = (unsigned) (g_ppy * (dLastDate-g_dOptimal)+0.5);

    switch (g_uiNumFactors)
    {
    case 1:
        uiStepsToFirstExpiry = NumSlicesBeforeFirstNotice1F(g_dSpecialDates.at(1));
        break;
    case 2:
        uiStepsToFirstExpiry = NumSlicesBeforeFirstNotice2F(g_dSpecialDates.at(1));
        break;
    case 3:
        uiStepsToFirstExpiry = NumSlicesBeforeFirstNotice3F(g_dSpecialDates.at(1));
        break;
    default:
        throw("ERROR: 1 or 2 factor only implemented so far");
        break;
    } // switch

    if ((uiStepsAfterFirstExpiry + uiStepsToFirstExpiry) < uiMinimum)
        uiStepsAfterFirstExpiry = uiMinimum - uiStepsToFirstExpiry-1;

    double 	dPreStepSize1=(dLastDate - g_dOptimal)/uiStepsAfterOptimal;
    dPreStepSize=(g_dOptimal - g_dSpecialDates.at(1))/uiStepsAfterFirstExpiryAndBeforeOptimal;

    datesArray[0] = 0.;
    register unsigned ii = 1;
    for(ii=1; ii < uiNumDates; ii++)
    {
        if((ii==1))
        {
            dStepSize = (g_dSpecialDates.at(1)-dStartDate)/(double)uiStepsToFirstExpiry;
        }
        else if(datesArray[g_uiNumSteps]<g_dOptimal)
            dStepSize = dPreStepSize;
        else
            dStepSize = dPreStepSize1;

        dStepGap = ((g_dSpecialDates.at(ii)-g_dSpecialDates.at(ii-1))/dStepSize);

        if(dStepSize<(g_dSpecialDates.at(ii)-g_dSpecialDates.at(ii-1))-1.e-08)
            dNewStepSize = (g_dSpecialDates.at(ii)-g_dSpecialDates.at(ii-1))/dStepGap;
        else
            dNewStepSize = (g_dSpecialDates.at(ii)-g_dSpecialDates.at(ii-1));

        while(datesArray[g_uiNumSteps]<(g_dSpecialDates.at(ii)-dStartDate)-1.e-08)
        {
            g_uiNumSteps++;
            if ((datesArray[g_uiNumSteps-1]+dNewStepSize>(g_dSpecialDates.at(ii)-dStartDate)-dNewStepSize/2.))
            {
                datesArray[g_uiNumSteps] = (g_dSpecialDates.at(ii)-dStartDate);
            }
            else
            {
                datesArray[g_uiNumSteps]=datesArray[g_uiNumSteps-1]+dNewStepSize;
            }
        } // for (g_uiNumSteps...)
    } // for (ii<=uiNumDates)

    // Copy over
    g_uiNumSteps++; // for the anchor
    dates.resize(g_uiNumSteps);
    for (ii=0; ii < g_uiNumSteps; ii++)
    {
        dates[ii] = datesArray[ii];
    }
} // CalcLatticeDates(...)

void LabelLatticeWithDates()
{

    unsigned uiLimit = g_uiNumSlices;

    g_LATTICE_DATE.resize(uiLimit);
    g_LATTICE_NORM.resize(uiLimit);
    g_LATTICE_DECAY.resize(uiLimit);
    if(g_uiProductModelCode>3)
    {
        g_LATTICE_NORM_FOREIGN.resize(uiLimit);
        g_LATTICE_DECAY_FOREIGN.resize(uiLimit);
    }
    g_LATTICE_PERIOD_TIME_SPAN.resize(uiLimit);

    // Label slices with dates
    unsigned ii = 0;
    for (ii = 0; ii < uiLimit; ii++)
    {
        // We use a double because we interpolated between dates to get here.
        g_LATTICE_DATE.at(ii)=g_dDates.at(ii);
    }

    double dt=0.;
    for (ii = 0; ii < uiLimit-1; ii++)
    {
        dt = g_LATTICE_DATE.at(ii+1) - g_LATTICE_DATE.at(ii);
        g_LATTICE_PERIOD_TIME_SPAN.at(ii)=dt;
    }
    g_LATTICE_PERIOD_TIME_SPAN.at(uiLimit-1)=dt;
} // LabelLatticeWithDates()


void InitLatticeModel()
{

    int r = g_uiNumFactors;
    int c = g_uiNumSlices;

    g_sigma.resize(r, c);
    g_meanRevRate.resize(r, c);

    // g_SigmaOriginal.resize(r,c);


    switch (g_uiNumFactors)
    {
    case 2:
        {
            if(g_uiProductModelCode==2)
            {
                g_S1=fabs(g_S1/g_dEquilibriumForwardBase.at(g_uiNumSlices-1));
                g_S2=fabs(g_S2/g_dEquilibriumForwardBase.at(g_uiNumSlices-1));
            }
            else if(g_uiProductModelCode==4)
            {
                g_S1=fabs(g_S1/g_dEquilibriumForwardBase.at(g_uiNumSlices-1));
                g_S2=fabs(g_S2/g_dEquilibriumForwardForeign.at(g_uiNumSlices-1));
            }
            else
            {
                throw("You should not be here");
            }

            double a_mr     = g_MR1;
            double sigma_1  = g_S1;
            double b_mr     = g_MR2;
            double sigma_2  = g_S2;

            g_dCorrelation  = g_RHO;
            g_dRelativeFactor = sigma_2/sigma_1;
            g_dCorrelationInitial=g_dCorrelation;

            double dEigenValue_1=0.;
            double dEigenValue_2=0.;
            double dRelative_1=0.;
            double dRelative_2=0.;
            double dNorm_1=0.;
            double dNorm_2=0.;

            dEigenValue_1=0.5*((1.+pow(g_dRelativeFactor,2.))+sqrt(pow((1.+pow(g_dRelativeFactor,2.)),2.)
                               -4.*pow(g_dRelativeFactor,2.)*(1.-pow(g_dCorrelation,2.))));
            dEigenValue_2=0.5*((1.+pow(g_dRelativeFactor,2.))-sqrt(pow((1.+pow(g_dRelativeFactor,2.)),2.)
                               -4.*pow(g_dRelativeFactor,2.)*(1.-pow(g_dCorrelation,2.))));

            dRelative_1=g_dRelativeFactor*g_dCorrelation/(dEigenValue_1-pow(g_dRelativeFactor,2.));
            dRelative_2=g_dRelativeFactor*g_dCorrelation/(dEigenValue_2-pow(g_dRelativeFactor,2.));

            dNorm_1=sqrt(1.+pow(dRelative_1,2.));
            dNorm_2=sqrt(1.+pow(dRelative_2,2.));

            g_dNormalisation_11=1./dNorm_1;
            g_dNormalisation_12=dRelative_1/dNorm_1;
            g_dNormalisation_21=1./dNorm_2;
            g_dNormalisation_22=dRelative_2/dNorm_2;

            g_dNormalisation_Inverse_11=(1./(g_dNormalisation_11/g_dNormalisation_12-g_dNormalisation_21/g_dNormalisation_22))/g_dNormalisation_12;
            g_dNormalisation_Inverse_12=-(1./(g_dNormalisation_11/g_dNormalisation_12-g_dNormalisation_21/g_dNormalisation_22))/g_dNormalisation_22;
            g_dNormalisation_Inverse_21=(1./(g_dNormalisation_12/g_dNormalisation_11-g_dNormalisation_22/g_dNormalisation_21))/g_dNormalisation_11;
            g_dNormalisation_Inverse_22=-(1./(g_dNormalisation_12/g_dNormalisation_11-g_dNormalisation_22/g_dNormalisation_21))/g_dNormalisation_21;

            // Keep records of original vols
            for (c = 0; c<g_uiNumSlices; ++c)
            {
                g_sigma.at(0,c)       = sigma_1;
                g_meanRevRate.at(0,c) = a_mr;
                g_sigma.at(1,c)       = sigma_2;
                g_meanRevRate.at(1,c) = b_mr;
            }

            // g_SigmaOriginal=g_sigma;

            g_dLimitY1.resize(g_uiNumSlices);
            g_dLimitY2.resize(g_uiNumSlices);

            // Set the standard deviations of the diagonal factors
            for (c = 0; c<g_uiNumSlices; ++c)
            {

                g_dLimitY1.at(c)=fabs(g_dNormalisation_11*g_dLimitXX1+g_dNormalisation_12*g_dRelativeFactor*g_dLimitXX2)/sqrt(pow(g_dNormalisation_11,2.)
                                 +pow(g_dNormalisation_12,2)*pow(g_dRelativeFactor,2.)+2*g_dRelativeFactor*g_dCorrelation*g_dNormalisation_11*g_dNormalisation_12);

                g_dLimitY2.at(c)=fabs(g_dNormalisation_21*g_dLimitXX1+g_dNormalisation_22*g_dRelativeFactor*g_dLimitXX2)/sqrt(pow(g_dNormalisation_21,2.)
                                 +pow(g_dNormalisation_22,2)*pow(g_dRelativeFactor,2.)+2*g_dRelativeFactor*g_dCorrelation*g_dNormalisation_21*g_dNormalisation_22);

                g_sigma.at(0,c)       = sigma_1 * sqrt (g_dNormalisation_11*g_dNormalisation_11
                                                        +g_dNormalisation_12*g_dNormalisation_12*g_dRelativeFactor*g_dRelativeFactor
                                                        +2.*g_dNormalisation_11*g_dNormalisation_12*g_dRelativeFactor*g_dCorrelation);
                g_sigma.at(1,c)       = sigma_1 * sqrt (g_dNormalisation_21*g_dNormalisation_21
                                                        +g_dNormalisation_22*g_dNormalisation_22*g_dRelativeFactor*g_dRelativeFactor
                                                        +2.*g_dNormalisation_21*g_dNormalisation_22*g_dRelativeFactor*g_dCorrelation);
            }
            // Set the correlation to zero
            g_dCorrelation=0.;
        } // if on 2-factors
        break;
    case 3:
        {
            if(g_uiProductModelCode==3)
            {
                g_S1=fabs(g_S1/g_dEquilibriumForwardBase.at(g_uiNumSlices-1));
                g_S2=fabs(g_S2/g_dEquilibriumForwardBase.at(g_uiNumSlices-1));
                g_S3=fabs(g_S3/g_dEquilibriumForwardBase.at(g_uiNumSlices-1));
            }
            else if(g_uiProductModelCode==5)
            {
                // g_S1=fabs(g_S1/g_dEquilibriumForwardBase.at(g_uiNumSlices-1));
                g_S1=g_S1_TD.at(g_uiNumSlices-1);
                // g_S2=fabs(g_S2/g_dEquilibriumForwardForeign.at(g_uiNumSlices-1));
                g_S2=g_S2_TD.at(g_uiNumSlices-1);
                // g_S3=g_S3; // no joy for FX skew
                g_S3=g_S3_TD.at(g_uiNumSlices-1);
            }
            else
            {
                throw("You should not be here");
            }

            double a_mr     = g_MR1;
            double sigma_1  = g_S1;
            double b_mr     = g_MR2;
            double sigma_2  = g_S2;
            double c_mr     = g_MR3;
            double sigma_3  = g_S3;

            g_dCorrelation  = g_RHO;
            g_dCorrelation_1 = g_RHO1;
            g_dCorrelation_2 = g_RHO2;

            DKMaille2D<double> dCovariance(3,3);
            dCovariance.at(0,0) = g_S1*g_S1; dCovariance.at(0,1) = g_RHO*g_S1*g_S2; dCovariance.at(0,2) = g_RHO1*g_S1*g_S3;
            dCovariance.at(1,0) = g_RHO*g_S1*g_S2; dCovariance.at(1,1) = g_S2*g_S2; dCovariance.at(1,2) = g_RHO2*g_S2*g_S3;
            dCovariance.at(2,0) = g_RHO1*g_S1*g_S3; dCovariance.at(2,1) = g_RHO2*g_S2*g_S3; dCovariance.at(2,2) = g_S3*g_S3;

            DKMaille<double> d(3);
            DKMaille2D<double> v(3,3);

            int nrot;

            LatticeJacobiTransformation(dCovariance,3,d,v,&nrot);
            EigenSort(d,v,3);

            double x11=v.at(0,0); double x12=v.at(1,0); double x13=v.at(2,0);
            double x21=v.at(0,1); double x22=v.at(1,1); double x23=v.at(2,1);
            double x31=v.at(0,2); double x32=v.at(1,2); double x33=v.at(2,2);

            g_d3Normalisation_11=x11; g_d3Normalisation_Inverse_11=x11;
            g_d3Normalisation_12=x12; g_d3Normalisation_Inverse_21=x12;
            g_d3Normalisation_13=x13; g_d3Normalisation_Inverse_31=x13;

            g_d3Normalisation_21=x21; g_d3Normalisation_Inverse_12=x21;
            g_d3Normalisation_22=x22; g_d3Normalisation_Inverse_22=x22;
            g_d3Normalisation_23=x23; g_d3Normalisation_Inverse_32=x23;

            g_d3Normalisation_31=x31; g_d3Normalisation_Inverse_13=x31;
            g_d3Normalisation_32=x32; g_d3Normalisation_Inverse_23=x32;
            g_d3Normalisation_33=x33; g_d3Normalisation_Inverse_33=x33;

            double dRatio2To1=d.at(1)/d.at(0);
            double dRatio3To1=d.at(2)/d.at(0);

            // Rotated factor vols

            double dVolY1=sqrt(pow(g_d3Normalisation_11*sigma_1,2.)
                               +pow(g_d3Normalisation_12*sigma_2,2.)
                               +pow(g_d3Normalisation_13*sigma_3,2.)
                               +2.*g_RHO*g_d3Normalisation_11*g_d3Normalisation_12*sigma_1*sigma_2
                               +2.*g_RHO1*g_d3Normalisation_11*g_d3Normalisation_13*sigma_1*sigma_3
                               +2.*g_RHO2*g_d3Normalisation_12*g_d3Normalisation_13*sigma_2*sigma_3);

            double dVolY2=sqrt(pow(g_d3Normalisation_21*sigma_1,2.)
                               +pow(g_d3Normalisation_22*sigma_2,2.)
                               +pow(g_d3Normalisation_23*sigma_3,2.)
                               +2.*g_RHO*g_d3Normalisation_21*g_d3Normalisation_22*sigma_1*sigma_2
                               +2.*g_RHO1*g_d3Normalisation_21*g_d3Normalisation_23*sigma_1*sigma_3
                               +2.*g_RHO2*g_d3Normalisation_22*g_d3Normalisation_23*sigma_2*sigma_3);

            double dVolY3=sqrt(pow(g_d3Normalisation_31*sigma_1,2.)
                               +pow(g_d3Normalisation_32*sigma_2,2.)
                               +pow(g_d3Normalisation_33*sigma_3,2.)
                               +2.*g_RHO*g_d3Normalisation_31*g_d3Normalisation_32*sigma_1*sigma_2
                               +2.*g_RHO1*g_d3Normalisation_31*g_d3Normalisation_33*sigma_1*sigma_3
                               +2.*g_RHO2*g_d3Normalisation_32*g_d3Normalisation_33*sigma_2*sigma_3);

            double dVolX1 =   sqrt(pow(g_d3Normalisation_Inverse_11*dVolY1,2.)
                                   +pow(g_d3Normalisation_Inverse_12*dVolY2,2.)
                                   +pow(g_d3Normalisation_Inverse_13*dVolY3,2.));

            double dVolX2 =   sqrt(pow(g_d3Normalisation_Inverse_21*dVolY1,2.)
                                   +pow(g_d3Normalisation_Inverse_22*dVolY2,2.)
                                   +pow(g_d3Normalisation_Inverse_23*dVolY3,2.));

            double dVolX3 =   sqrt(pow(g_d3Normalisation_Inverse_31*dVolY1,2.)
                                   +pow(g_d3Normalisation_Inverse_32*dVolY2,2.)
                                   +pow(g_d3Normalisation_Inverse_33*dVolY3,2.));

            double dCorrelation12=(g_d3Normalisation_Inverse_11*g_d3Normalisation_Inverse_21*dVolY1*dVolY1
                                   +g_d3Normalisation_Inverse_12*g_d3Normalisation_Inverse_22*dVolY2*dVolY2
                                   +g_d3Normalisation_Inverse_13*g_d3Normalisation_Inverse_23*dVolY3*dVolY3)/(dVolX1*dVolX2);

            double dCorrelation13=(g_d3Normalisation_Inverse_11*g_d3Normalisation_Inverse_31*dVolY1*dVolY1
                                   +g_d3Normalisation_Inverse_12*g_d3Normalisation_Inverse_32*dVolY2*dVolY2
                                   +g_d3Normalisation_Inverse_13*g_d3Normalisation_Inverse_33*dVolY3*dVolY3)/(dVolX1*dVolX3);

            double dCorrelation23=(g_d3Normalisation_Inverse_21*g_d3Normalisation_Inverse_31*dVolY1*dVolY1
                                   +g_d3Normalisation_Inverse_22*g_d3Normalisation_Inverse_32*dVolY2*dVolY2
                                   +g_d3Normalisation_Inverse_23*g_d3Normalisation_Inverse_33*dVolY3*dVolY3)/(dVolX2*dVolX3);

            g_dRelativeFactor = sigma_2/sigma_1;

            g_dCorrelationInitial=g_dCorrelation;

            double dEigenValue_1=0.;
            double dEigenValue_2=0.;
            double dRelative_1=0.;
            double dRelative_2=0.;
            double dNorm_1=0.;
            double dNorm_2=0.;

            dEigenValue_1=0.5*((1.+pow(g_dRelativeFactor,2.))+sqrt(pow((1.+pow(g_dRelativeFactor,2.)),2.)
                               -4.*pow(g_dRelativeFactor,2.)*(1.-pow(g_dCorrelation,2.))));
            dEigenValue_2=0.5*((1.+pow(g_dRelativeFactor,2.))-sqrt(pow((1.+pow(g_dRelativeFactor,2.)),2.)
                               -4.*pow(g_dRelativeFactor,2.)*(1.-pow(g_dCorrelation,2.))));

            dRelative_1=g_dRelativeFactor*g_dCorrelation/(dEigenValue_1-pow(g_dRelativeFactor,2.));
            dRelative_2=g_dRelativeFactor*g_dCorrelation/(dEigenValue_2-pow(g_dRelativeFactor,2.));

            dNorm_1=sqrt(1.+pow(dRelative_1,2.));
            dNorm_2=sqrt(1.+pow(dRelative_2,2.));

            g_dNormalisation_11=1./dNorm_1;
            g_dNormalisation_12=dRelative_1/dNorm_1;
            g_dNormalisation_21=1./dNorm_2;
            g_dNormalisation_22=dRelative_2/dNorm_2;

            g_dNormalisation_Inverse_11=(1./(g_dNormalisation_11/g_dNormalisation_12
                                             -g_dNormalisation_21/g_dNormalisation_22))/g_dNormalisation_12;
            g_dNormalisation_Inverse_12=-(1./(g_dNormalisation_11/g_dNormalisation_12
                                              -g_dNormalisation_21/g_dNormalisation_22))/g_dNormalisation_22;
            g_dNormalisation_Inverse_21=(1./(g_dNormalisation_12/g_dNormalisation_11
                                             -g_dNormalisation_22/g_dNormalisation_21))/g_dNormalisation_11;
            g_dNormalisation_Inverse_22=-(1./(g_dNormalisation_12/g_dNormalisation_11
                                              -g_dNormalisation_22/g_dNormalisation_21))/g_dNormalisation_21;

            // Keep records of original vols
            for (c = 0; c<g_uiNumSlices; ++c)
            {
                g_sigma.at(0,c)       = sigma_1;
                g_meanRevRate.at(0,c) = a_mr;
                g_sigma.at(1,c)       = sigma_2;
                g_meanRevRate.at(1,c) = b_mr;
                g_sigma.at(2,c)       = sigma_3;
                if(g_uiProductModelCode!=5) g_meanRevRate.at(2,c) = c_mr;
                else g_meanRevRate.at(2,c)=0.;
            }

            // g_SigmaOriginal=g_sigma;

            g_dLimitY1.resize(g_uiNumSlices);
            g_dLimitY2.resize(g_uiNumSlices);
            g_dLimitY3.resize(g_uiNumSlices);

            for (c = 0; c<g_uiNumSlices; ++c)
            {
                g_sigma.at(0,c)       = dVolY1;
                g_sigma.at(1,c)       = dVolY2;
                g_sigma.at(2,c)       = dVolY3;


                if(g_uiProductModelCode==5) g_dLimitY1.at(c)=(fabs(g_d3Normalisation_11*g_dLimitXX1*sigma_1/**exp(-g_MR1*g_LATTICE_DATE.at(c))*/)
                            +fabs(g_d3Normalisation_12*g_dLimitXX2*sigma_2/**exp(-g_MR2*g_LATTICE_DATE.at(c))*/)
                            +fabs(g_d3Normalisation_13*g_dLimitXX3*sigma_3))/dVolY1;
                else g_dLimitY1.at(c)=(fabs(g_d3Normalisation_11*g_dLimitXX1*sigma_1/**exp(-g_MR1*g_LATTICE_DATE.at(c))*/)
                                           +fabs(g_d3Normalisation_12*g_dLimitXX2*sigma_2/**exp(-g_MR2*g_LATTICE_DATE.at(c))*/)
                                           +fabs(g_d3Normalisation_13*g_dLimitXX3*sigma_3/**exp(-g_MR3*g_LATTICE_DATE.at(c))*/))/dVolY1;

                if(g_uiProductModelCode==5) g_dLimitY2.at(c)=(fabs(g_d3Normalisation_21*g_dLimitXX1*sigma_1/**exp(-g_MR1*g_LATTICE_DATE.at(c))*/)
                            +fabs(g_d3Normalisation_22*g_dLimitXX2*sigma_2/**exp(-g_MR2*g_LATTICE_DATE.at(c))*/)
                            +fabs(g_d3Normalisation_23*g_dLimitXX3*sigma_3))/dVolY2;
                else g_dLimitY2.at(c)=(fabs(g_d3Normalisation_21*g_dLimitXX1*sigma_1/**exp(-g_MR1*g_LATTICE_DATE.at(c))*/)
                                           +fabs(g_d3Normalisation_22*g_dLimitXX2*sigma_2/**exp(-g_MR2*g_LATTICE_DATE.at(c))*/)
                                           +fabs(g_d3Normalisation_23*g_dLimitXX3*sigma_3/**exp(-g_MR3*g_LATTICE_DATE.at(c))*/))/dVolY2;

                if(g_uiProductModelCode==5) g_dLimitY3.at(c)=(fabs(g_d3Normalisation_31*g_dLimitXX1*sigma_1/**exp(-g_MR1*g_LATTICE_DATE.at(c))*/)
                            +fabs(g_d3Normalisation_32*g_dLimitXX2*sigma_2/**exp(-g_MR2*g_LATTICE_DATE.at(c))*/)
                            +fabs(g_d3Normalisation_33*g_dLimitXX3*sigma_3))/dVolY3;
                else g_dLimitY3.at(c)=(fabs(g_d3Normalisation_31*g_dLimitXX1*sigma_1/**exp(-g_MR1*g_LATTICE_DATE.at(c))*/)
                                           +fabs(g_d3Normalisation_32*g_dLimitXX2*sigma_2/**exp(-g_MR2*g_LATTICE_DATE.at(c))*/)
                                           +fabs(g_d3Normalisation_33*g_dLimitXX3*sigma_3/**exp(-g_MR3*g_LATTICE_DATE.at(c))*/))/dVolY3;
            }
        } // if on 3-factors
        break;
    case 1:
        {
            // g_S1      =  fabs(g_S1/g_dEquilibriumForwardBase.at(g_uiNumSlices-1));
            g_S1 = g_S1_TD.at(g_uiNumSlices-1);
            g_dLimitY1.resize(g_uiNumSlices);
            double dForward=0.;
            for (c = 0; c<g_uiNumSlices; ++c)
            {
                g_dLimitY1.at(c)=g_dLimitXX1/**exp(-g_MR1*g_LATTICE_DATE.at(c))*/;
                g_sigma.at(0,c)       = g_S1;
                g_meanRevRate.at(0,c) = g_MR1;
            }
            // g_SigmaOriginal=g_sigma;
        }
        break;
    default:
        {
            throw("ERROR: InitLatticeModel(): You have to take care of youself !");
        }
    };
} // InitLatticeModel()

double BaseDiscountInterpolate(double dDate)
{
    double dTaux=rateinterpolation_dk_maille(2,dDate,g_dDiscountCurveBaseDates,
                 g_dZCCurveBase,g_dDiscountCurveBaseDates.entries()-1);
    return exp(-dTaux*dDate);
}

double ForeignDiscountInterpolate(double dDate)
{
    double dTaux=rateinterpolation_dk_maille(2,dDate,g_dDiscountCurveForeignDates,
                 g_dZCCurveForeign,g_dDiscountCurveForeign.entries()-1);
    return exp(-dTaux*dDate);
}

void SetLatticeMarketDf(unsigned int uiSlice)
{
    double dDate=0.;
    double dMarketDf=0.;
    if(uiSlice<=g_uiNumSlices-1)
    {
        dDate = g_LATTICE_DATE.at(uiSlice);
        dMarketDf = BaseDiscountInterpolate(dDate);
        g_LATTICE_ZCB_PRICE.at(uiSlice)=dMarketDf;
    }
    else
    {
        dDate=g_LATTICE_DATE.at(uiSlice-1)+g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice-1);
        dMarketDf = BaseDiscountInterpolate(dDate);
        g_LATTICE_ZCB_PRICE.at(uiSlice)=dMarketDf;
    }
} // SetLatticeMarketDf(...)


void SetLatticeMarketDfForeign(unsigned int uiSlice)
{
    double dDate=0.;
    double dMarketDf=0.;
    if(uiSlice<=g_uiNumSlices-1)
    {
        dDate = g_LATTICE_DATE.at(uiSlice);
        dMarketDf = ForeignDiscountInterpolate(dDate);
        g_LATTICE_ZCB_PRICE_FOREIGN.at(uiSlice)=dMarketDf;
    }
    else
    {
        dDate=g_LATTICE_DATE.at(uiSlice-1)+g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice-1);
        dMarketDf = ForeignDiscountInterpolate(dDate);
        g_LATTICE_ZCB_PRICE_FOREIGN.at(uiSlice)=dMarketDf;
    }
} // SetLatticeMarketDfForeign(...)


double lattice_sigma(unsigned uiFactorIndex, unsigned uiSliceIndex)
{
    return g_sigma.at(uiFactorIndex-1, uiSliceIndex);
} // lattice_sigma(...)

DKMaille<double>& CalcLatticeRateStep(unsigned int uiSlice, double dt)
{
    static DKMaille<double> rateStep(g_uiNumFactors);
    rateStep.resize(g_uiNumFactors);

    double dFactor=1.;

    g_GridScaling1=1.;
    g_GridScaling2=1.;
    g_GridScaling3=1.;

    switch (g_uiNumFactors)
    {
    case 3:
        rateStep[2] = g_GridScaling3*lattice_sigma(3, uiSlice)*sqrt(3.0*dt);
    case 2:
        rateStep[1] = g_GridScaling2*lattice_sigma(2, uiSlice)* sqrt(3.0*dt);
    case 1:
        rateStep[0] = g_GridScaling1*lattice_sigma(1, uiSlice)*sqrt(3.0*dt);
        break;
    default:
        throw("ERROR: 1,2 or 3 factors only!");
        break;
    }
    return rateStep;
} // CalcLatticeRateStep(...)

DKMaille<double> SetLatticeRateStep(unsigned int uiSlice, double dt)
{
    DKMaille<double> &rateSteps = CalcLatticeRateStep(uiSlice, dt);
    switch (g_uiNumFactors)
    {
    case 3:
        g_LATTICE_RATE_STEP3.at(uiSlice)=rateSteps[2];
    case 2:
        g_LATTICE_RATE_STEP2.at(uiSlice)=rateSteps[1];
    case 1:
        g_LATTICE_RATE_STEP1.at(uiSlice)=rateSteps[0];
        break;
    default:
        throw("ERROR: 1,2 or 3 factors only!");
        break;
    }
    return rateSteps;
} // SetLatticeRateStep


double LatticeStateToXRate(unsigned int uiSlice, int uiState, unsigned int uiFactorIndex)
{
    double dr;
    switch (uiFactorIndex)
    {
    case 3:
        dr = g_LATTICE_RATE_STEP3.at(uiSlice);
        break;
    case 2:
        dr = g_LATTICE_RATE_STEP2.at(uiSlice);
        break;
    case 1:
        dr = g_LATTICE_RATE_STEP1.at(uiSlice);
        break;
    default:
        throw("ERROR: 1,2 or 3 factors only!");
        break;
    };
    return dr * uiState;
}// LatticeStateToXRate



double LatticeFromYRateToXRate(DKMaille<double> &dYRate,unsigned int uiFactor)
{
    double dRes=0.;
    if(g_uiNumFactors==2)
    {
        if(uiFactor==1)
        {
            dRes=dYRate.at(0)*g_dNormalisation_Inverse_11
                 +dYRate.at(1)*g_dNormalisation_Inverse_12;

        }
        if(uiFactor==2)
        {
            dRes=dYRate.at(0)*g_dNormalisation_Inverse_21
                 +dYRate.at(1)*g_dNormalisation_Inverse_22;
        }
    }
    if(g_uiNumFactors==3)
    {
        if(uiFactor==1)
        {
            dRes=dYRate.at(0)*g_d3Normalisation_Inverse_11
                 +dYRate.at(1)*g_d3Normalisation_Inverse_12
                 +dYRate.at(2)*g_d3Normalisation_Inverse_13;

        }
        if(uiFactor==2)
        {
            dRes=dYRate.at(0)*g_d3Normalisation_Inverse_21
                 +dYRate.at(1)*g_d3Normalisation_Inverse_22
                 +dYRate.at(2)*g_d3Normalisation_Inverse_23;
        }
        if(uiFactor==3)
        {
            dRes=dYRate.at(0)*g_d3Normalisation_Inverse_31
                 +dYRate.at(1)*g_d3Normalisation_Inverse_32
                 +dYRate.at(2)*g_d3Normalisation_Inverse_33;
        }
    }

    return dRes;
} // LatticeFromYRateToXRate

double LatticeXToGRate(DKMaille<double> &dXRate,unsigned int uiSlice)
{
    double dG=0.;
    switch (g_uiNumFactors)
    {
    case 3:
        {
            double dXX1 = LatticeFromYRateToXRate(dXRate,1)*g_S1_TD.at(uiSlice)/g_S1;
            double dXX2 = LatticeFromYRateToXRate(dXRate,2)*g_S2_TD.at(uiSlice)/g_S2;
            double dXX3 = LatticeFromYRateToXRate(dXRate,3)*g_S3_TD.at(uiSlice)/g_S3;
            dG=dXX1+dXX2+dXX3;
        }
        break;
    case 2:
        {
            double dXX1 = LatticeFromYRateToXRate(dXRate,1)*g_S1_TD.at(uiSlice)/g_S1;
            double dXX2 = LatticeFromYRateToXRate(dXRate,2)*g_S2_TD.at(uiSlice)/g_S2;
            dG=dXX1+dXX2;
        }
        break;
    case 1:
        {
            dG=dXRate.at(0)*g_S1_TD.at(uiSlice)/g_S1;
        }
        break;
    }
    return dG;
}// LatticeXToGRate

// Passing the return value in as a reference parameter avoids costly
// calls to the DKMaille ctor/dtor and array copies.
void LatticeXToGRateHybrid(DKMaille<double> &dXRate, DKMaille<double> &dG,unsigned int uiSlice)
{
    // This assummes 2 or 3 factors in quanto or fx style
    switch (g_uiNumFactors)
    {
    case 3:
        {
            double dXX1 = LatticeFromYRateToXRate(dXRate,1)*g_S1_TD.at(uiSlice)/g_S1;
            double dXX2 = LatticeFromYRateToXRate(dXRate,2)*g_S2_TD.at(uiSlice)/g_S2;
            double dXX3 = LatticeFromYRateToXRate(dXRate,3)*g_S3_TD.at(uiSlice)/g_S3;
            dG.at(0)=dXX1;
            dG.at(1)=dXX2;
            dG.at(2)=dXX3;
        }
        break;
    case 2:
        {
            double dXX1 = LatticeFromYRateToXRate(dXRate,1)*g_S1_TD.at(uiSlice)/g_S1;
            double dXX2 = LatticeFromYRateToXRate(dXRate,2)*g_S2_TD.at(uiSlice)/g_S2;
            dG.at(0)=dXX1;
            dG.at(1)=dXX2;
        }
        break;
    case 1:
        throw("LatticeXToGRateHybrid:: You should not have come here");
        break;
    }
}// LatticeXToGRateHybrid

double DKSkewFunction(double dSkewParameter, double dGRate,double dForward)
{
    if(fabs(dForward-0.002)<1.e-09) dSkewParameter=0.;
    double x=0.;
    if(dSkewParameter<1.e-10)
    {
        x=dForward*(1.+dGRate);
    }
    else
    {
        x = dForward * (1.+(1./dSkewParameter)*(exp(dSkewParameter*dGRate)-1.));
    }
    return x;
}// DKSkewFunction


double DKSkewFunctionDerivative(double dSkewParameter, double dGRate,double dForward)
{
    if(fabs(dForward-0.002)<1.e-09) dSkewParameter=0.;
    double x=0.;
    if(dSkewParameter<1.e-10)
    {
        x = dForward;
    }
    else
    {
        x = dForward * exp(dSkewParameter*dGRate);
    }
    return x;
}// DKSkewFunctionDerivative


double LatticeGToShortRate(const double dGRate,unsigned int uiSlice)
{
    // In principle X and dGRate are different depending on the distribution. In the HW model X=dGRate
    // The below does not need modification for multi-factor model
    double x = DKSkewFunction(g_dQParameter_Base,dGRate,g_dEquilibriumForwardBase.at(uiSlice));
    return x;
}// LatticeGToShortRate


double LatticeGToShortRateDerivative(const double dGRate,unsigned int uiSlice)
{
    // In principle X and dGRate are different depending on the distribution. In the HW model X=dGRate
    // The below does not need modification for multi-factor model
    double x = DKSkewFunctionDerivative(g_dQParameter_Base,dGRate,g_dEquilibriumForwardBase.at(uiSlice));
    return x;
}// LatticeGToShortRateDerivative


// Return x value as a reference parameter to avoid expensive array copies
// and ctor/dtor calls.
void LatticeGToShortRateFXHybrid(DKMaille<double> &dGRate, DKMaille<double> &x, unsigned int uiSlice)
{
    switch(g_uiNumFactors)
    {
    case 3:
        {
            x.at(0)=DKSkewFunction(g_dQParameter_Base,dGRate.at(0),g_dEquilibriumForwardBase.at(uiSlice));
            x.at(1)=DKSkewFunction(g_dQParameter_Foreign,dGRate.at(1),g_dEquilibriumForwardForeign.at(uiSlice));
            // x.at(2)=DKSkewFunction(1.0,dGRate.at(2));
            // First implementation FX is lognormal
            // Third line not needed
        }
        break;
    case 2:
        {
            x.at(0)=DKSkewFunction(g_dQParameter_Base,dGRate.at(0),g_dEquilibriumForwardBase.at(uiSlice));
            x.at(1)=DKSkewFunction(g_dQParameter_Foreign,dGRate.at(1),g_dEquilibriumForwardForeign.at(uiSlice));
        }
        break;
    case 1:
        {
            throw("LatticeGToShortRateFXHybrid:: You should not be here");
        }
        break;
    }
}// LatticeGToShortRateFXHybrid




unsigned int LatticeIndexConvert(unsigned int uiSlice,DKMaille<long> &states)
{
    unsigned int uiIndex;
    switch(g_uiNumFactors)
    {
    case 1:
        {
            unsigned int uiIndex1=states.at(0)-g_LATTICE_MIN_STATE1.at(uiSlice);
            uiIndex=uiIndex1;
        }
        break;
    case 2:
        {
            unsigned int uiIndex1=states.at(0)-g_LATTICE_MIN_STATE1.at(uiSlice);
            unsigned int uiIndex2=states.at(1)-g_LATTICE_MIN_STATE2.at(uiSlice);
            uiIndex=uiIndex1*(g_LATTICE_NUMBER_OF_STATES2.at(uiSlice))+uiIndex2;
        }
        break;
    case 3:
        {
            unsigned int uiIndex1=states.at(0)-g_LATTICE_MIN_STATE1.at(uiSlice);
            unsigned int uiIndex2=states.at(1)-g_LATTICE_MIN_STATE2.at(uiSlice);
            unsigned int uiIndex3=states.at(2)-g_LATTICE_MIN_STATE3.at(uiSlice);
            uiIndex=uiIndex1*(g_LATTICE_NUMBER_OF_STATES2.at(uiSlice))*(g_LATTICE_NUMBER_OF_STATES3.at(uiSlice))
                    +uiIndex2*(g_LATTICE_NUMBER_OF_STATES3.at(uiSlice))+uiIndex3;
        }
        break;
    }
    return uiIndex;
}

void LabelLatticeNodeWithXRate(unsigned int uiSlice,unsigned int uiIndex,DKMaille<double> &dXRate)
{
    switch(g_uiNumFactors)
    {
    case 1:
        {
            g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex)=dXRate.at(0);
        }
        break;
    case 2:
        {
            g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex)=dXRate.at(0);
            g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex)=dXRate.at(1);
        }
        break;
    case 3:
        {
            g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex)=dXRate.at(0);
            g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex)=dXRate.at(1);
            g_LATTICE_X3_RATE.at(uiSlice).at(uiIndex)=dXRate.at(2);
        }
        break;
    }
}

void LabelLatticeNodeWithGRate(unsigned int uiSlice,unsigned int uiIndex,double dGRate)
{
    g_LATTICE_G_RATE.at(uiSlice).at(uiIndex)=dGRate;
}

void LabelLatticeNodeWithGRateHybrid(unsigned int uiSlice,unsigned int uiIndex,DKMaille<double> &dGRate,unsigned int uiFactorIndex)
{
    if(uiFactorIndex==1) g_LATTICE_G_RATE.at(uiSlice).at(uiIndex)=dGRate.at(0);
    if(uiFactorIndex==2) g_LATTICE_G_RATE_FOREIGN.at(uiSlice).at(uiIndex)=dGRate.at(1);
    if(uiFactorIndex==3) g_LATTICE_G_RATE_FX.at(uiSlice).at(uiIndex)=dGRate.at(2);
}

void LabelLatticeNodeWithShortRate(unsigned int uiSlice,unsigned int uiIndex,double dShortRate)
{
    if(dShortRate<5.) g_LATTICE_SHORT_RATE.at(uiSlice).at(uiIndex)=dShortRate;
    else g_LATTICE_SHORT_RATE.at(uiSlice).at(uiIndex)=5.;
}

void LabelLatticeNodeWithShortRateFXHybrid(unsigned int uiSlice,unsigned int uiIndex,DKMaille<double> &dShortRate,unsigned int uiFactorIndex)
{
    if(uiFactorIndex==1)
    {
        if(dShortRate.at(0)<5.) g_LATTICE_SHORT_RATE.at(uiSlice).at(uiIndex)=dShortRate.at(0);
        else g_LATTICE_SHORT_RATE.at(uiSlice).at(uiIndex)=5.;
    }
    if(uiFactorIndex==2)
    {
        if(dShortRate.at(0)<5.) g_LATTICE_SHORT_RATE_FOREIGN.at(uiSlice).at(uiIndex)=dShortRate.at(1);
        else g_LATTICE_SHORT_RATE_FOREIGN.at(uiSlice).at(uiIndex)=5.;
    }
    // if(uiFactorIndex==3) g_LATTICE_SHORT_RATE_FX.at(uiSlice).at(uiIndex)=g_dSpotFX*dShortRate.at(2);
    // third line not needed
}

void PrepareTheLatticeNode(unsigned int uiSlice, DKMaille<long> &states)
{
    // Necessary for all slices
    static DKMaille<double> dXRate(g_uiNumFactors);
    dXRate.resize(g_uiNumFactors);
    for(unsigned int uiFactor=0;uiFactor<g_uiNumFactors;uiFactor++)
        dXRate.at(uiFactor)= LatticeStateToXRate(uiSlice, states.at(uiFactor),uiFactor+1);
    unsigned int uiIndex=LatticeIndexConvert(uiSlice,states);
    LabelLatticeNodeWithXRate(uiSlice, uiIndex, dXRate);
    g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex)=true;
    // Necessary for last slice
    if(uiSlice==g_uiNumSlices-1)
    {
        double dGRate = 0.;
        DKMaille<double> dGRateHybrid(g_uiNumFactors);
        if(g_uiProductModelCode<4) dGRate = LatticeXToGRate(dXRate,uiSlice);
        else LatticeXToGRateHybrid(dXRate, dGRateHybrid, uiSlice);
        static DKMaille<double> shortRateHybrid(g_uiNumFactors);
        shortRateHybrid.resize(g_uiNumFactors);
        double shortRate = LatticeGToShortRate(dGRate,uiSlice);
        if(g_uiProductModelCode>3)
        {
            LatticeGToShortRateFXHybrid(dGRateHybrid, shortRateHybrid, uiSlice);
        }

        // Convert to DK key count
        // Set gRate & shortRate
        if(g_uiProductModelCode<4)
        {
            LabelLatticeNodeWithGRate(uiSlice, uiIndex, dGRate);
            LabelLatticeNodeWithShortRate(uiSlice, uiIndex, shortRate);
        }
        else
        {
            for(unsigned uiK=1;uiK<=g_uiNumFactors;uiK++)
            {
                LabelLatticeNodeWithGRateHybrid(uiSlice, uiIndex, dGRateHybrid,uiK);
            }
            LabelLatticeNodeWithShortRateFXHybrid(uiSlice, uiIndex, shortRateHybrid,1);
            LabelLatticeNodeWithShortRateFXHybrid(uiSlice, uiIndex, shortRateHybrid,2);
        }
    }
}// PrepareTheLatticeNode(...)

void LatticeInitMaxMinX(unsigned int uiSlice)
{
    switch(g_uiNumFactors)
    {
    case 3:
        {
            g_LATTICE_MAX_STATE3.at(uiSlice)=0;
            g_LATTICE_MIN_STATE3.at(uiSlice)=0;
        }
    case 2:
        {
            g_LATTICE_MAX_STATE2.at(uiSlice)=0;
            g_LATTICE_MIN_STATE2.at(uiSlice)=0;
        }
    case 1:
        {
            g_LATTICE_MAX_STATE1.at(uiSlice)=0;
            g_LATTICE_MIN_STATE1.at(uiSlice)=0;
        }
        break;
    }
}


void LatticeSetMaxMinX(DKMaille<long> &iTargetStates,unsigned int uiSlice)
{
    static long dMax, dMin, dVal1, dVal2;

    switch(iTargetStates.entries())
    {
    case 3:
        {
            dMax = g_LATTICE_MAX_STATE3.at(uiSlice);
            dMin = g_LATTICE_MIN_STATE3.at(uiSlice);
            dVal1 = iTargetStates.at(2)+1;
            dVal2 = iTargetStates.at(2)-1;
            if (dVal1 > dMax)
            {
                g_LATTICE_MAX_STATE3.at(uiSlice)=dVal1;
            }
            if (dVal2 < dMin)
            {
                g_LATTICE_MIN_STATE3.at(uiSlice)=dVal2;
            }
            g_LATTICE_NUMBER_OF_STATES3.at(uiSlice)=g_LATTICE_MAX_STATE3.at(uiSlice)-g_LATTICE_MIN_STATE3.at(uiSlice)+1;
        }
    case 2:
        {
            dMax = g_LATTICE_MAX_STATE2.at(uiSlice);
            dMin = g_LATTICE_MIN_STATE2.at(uiSlice);
            dVal1 = iTargetStates.at(1)+1;
            dVal2 = iTargetStates.at(1)-1;

            if (dVal1 > dMax)
            {
                g_LATTICE_MAX_STATE2.at(uiSlice)=dVal1;
            }
            if (dVal2 < dMin)
            {
                g_LATTICE_MIN_STATE2.at(uiSlice)=dVal2;
            }
            g_LATTICE_NUMBER_OF_STATES2.at(uiSlice)=g_LATTICE_MAX_STATE2.at(uiSlice)-g_LATTICE_MIN_STATE2.at(uiSlice)+1;
        }
    case 1:
        {
            dMax = g_LATTICE_MAX_STATE1.at(uiSlice);
            dMin = g_LATTICE_MIN_STATE1.at(uiSlice);
            dVal1 = iTargetStates.at(0)+1;
            dVal2 = iTargetStates.at(0)-1;

            if (dVal1 > dMax)
            {
                g_LATTICE_MAX_STATE1.at(uiSlice)=dVal1;
            }
            if (dVal2 < dMin)
            {
                g_LATTICE_MIN_STATE1.at(uiSlice)=dVal2;
            }
            g_LATTICE_NUMBER_OF_STATES1.at(uiSlice)=g_LATTICE_MAX_STATE1.at(uiSlice)-g_LATTICE_MIN_STATE1.at(uiSlice)+1;
        }
        break;
    default:
        throw("ERROR: Wrong number of XRates supplied to LatticeSetMaxMinX");
        break;
    } // switch
} // LatticeSetMaxMinX(...)

double LatticeShortRateForwardMeanNew(DKMaille<double> &dXRate,
                                      double dt,
                                      unsigned uiFactorIndex,
                                      unsigned uiSliceIndex,
                                      unsigned int uiNodeIndex)
{
    double dRes=0.;
    // single currency fixed income
    if(g_uiProductModelCode<4)
    {
        // condition on whether we are dealing with unrotated or rotated to the diagonal factors
        if(g_uiNumFactors==2)
        {
            // We are in rotated world
            // What comes in as dXRate is actually dYRate (orthogonal factors)
            // from Y1, Y2 to X1, X2
            // Rotation formulas
            double dXX1 = LatticeFromYRateToXRate(dXRate,1);
            double dXX2 = LatticeFromYRateToXRate(dXRate,2);
            if(uiFactorIndex==1) dRes=g_dDiffusion_11.at(uiSliceIndex)*dXX1
                                          +g_dDiffusion_12.at(uiSliceIndex)*dXX2;
            if(uiFactorIndex==2) dRes=g_dDiffusion_21.at(uiSliceIndex)*dXX1
                                          +g_dDiffusion_22.at(uiSliceIndex)*dXX2;
        }
        else if(g_uiNumFactors==3)
        {
            // We are in rotated world
            // What comes in as dXRate is actually dYRate (orthogonal factors)
            // from Y1, Y2 to X1, X2
            // Rotation formulas
            double dXX1 = LatticeFromYRateToXRate(dXRate,1);
            double dXX2 = LatticeFromYRateToXRate(dXRate,2);
            double dXX3 = LatticeFromYRateToXRate(dXRate,3);

            if(uiFactorIndex==1) dRes=g_dDiffusion_11.at(uiSliceIndex)*dXX1
                                          +g_dDiffusion_12.at(uiSliceIndex)*dXX2
                                          +g_dDiffusion_13.at(uiSliceIndex)*dXX3;

            if(uiFactorIndex==2) dRes=g_dDiffusion_21.at(uiSliceIndex)*dXX1
                                          +g_dDiffusion_22.at(uiSliceIndex)*dXX2
                                          +g_dDiffusion_23.at(uiSliceIndex)*dXX3;

            if(uiFactorIndex==3) dRes=g_dDiffusion_31.at(uiSliceIndex)*dXX1
                                          +g_dDiffusion_32.at(uiSliceIndex)*dXX2
                                          +g_dDiffusion_33.at(uiSliceIndex)*dXX3;
        }
        else
        {
            dRes=g_dDiffusion_11.at(uiSliceIndex)*dXRate.at(uiFactorIndex-1);
        }
    }
    else // we are in hybrid mode
    {
        // condition on whether we are dealing with unrotated or rotated to the diagonal factors
        if(g_uiNumFactors==2)
        {
            // We are in rotated world
            // What comes in as dXRate is actually dYRate (orthogonal factors)
            // from Y1, Y2 to X1, X2
            // Rotation formulas
            double dXX1 = LatticeFromYRateToXRate(dXRate,1);
            double dXX2 = LatticeFromYRateToXRate(dXRate,2);
            if(uiFactorIndex==1) dRes=g_dDiffusion_11.at(uiSliceIndex)*dXX1
                                          +g_dDiffusion_12.at(uiSliceIndex)*dXX2;
            if(uiFactorIndex==2) dRes=g_dDiffusion_21.at(uiSliceIndex)*dXX1
                                          +g_dDiffusion_22.at(uiSliceIndex)*dXX2;

        }
        else if(g_uiNumFactors==3)
        {
            // We are in rotated world
            // What comes in as dXRate is actually dYRate (orthogonal factors)
            // from Y1, Y2 to X1, X2
            // Rotation formulas

            // Drift of the IR processes
            double dXX1 = LatticeFromYRateToXRate(dXRate,1);
            double dXX2 = LatticeFromYRateToXRate(dXRate,2);
            double dXX3 = LatticeFromYRateToXRate(dXRate,3);

            // Drift of the FX process
            double dShortRateBase=log(LatticeDiscount(uiSliceIndex,uiNodeIndex))/(-dt);
            double dShortRateForeign=log(LatticeDiscountForeign(uiSliceIndex,uiNodeIndex))/(-dt);

            double dLogFXDrift=0.;

            // Volatility contribution term taken care of by string calibration in FillFX
            if(g_dADPLimit==1.0)
            {
                dLogFXDrift=(-pow(g_S3_TD.at(uiSliceIndex),2.)*dt*0.5+(dShortRateBase-dShortRateForeign)*dt)*g_S3/g_S3_TD.at(uiSliceIndex);
            }
            else
            {
                dLogFXDrift=(+(dShortRateBase-dShortRateForeign)*dt)*g_S3/g_S3_TD.at(uiSliceIndex);
            }
            if(uiFactorIndex==1) dRes=g_dDiffusion_11.at(uiSliceIndex)*dXX1
                                          +g_dDiffusion_12.at(uiSliceIndex)*dXX2
                                          +g_dDiffusion_13.at(uiSliceIndex)*dXX3
                                          +g_d3Normalisation_13*dLogFXDrift;

            if(uiFactorIndex==2) dRes=g_dDiffusion_21.at(uiSliceIndex)*dXX1
                                          +g_dDiffusion_22.at(uiSliceIndex)*dXX2
                                          +g_dDiffusion_23.at(uiSliceIndex)*dXX3
                                          +g_d3Normalisation_23*dLogFXDrift;

            if(uiFactorIndex==3) dRes=g_dDiffusion_31.at(uiSliceIndex)*dXX1
                                          +g_dDiffusion_32.at(uiSliceIndex)*dXX2
                                          +g_dDiffusion_33.at(uiSliceIndex)*dXX3
                                          +g_d3Normalisation_33*dLogFXDrift;
        }
        else
        {
            throw("LatticeShortRateForwardMeanNew:: You should never come here");
        }
    }
    return dRes;
} // LatticeShortRateForwardMeanNew(...)



void LabelLatticeLangevinDiffusion()
{
    double dRes=0.;
    // single currency fixed income
    double dt=0.;
    for(unsigned int uiSliceIndex=0;uiSliceIndex<g_uiNumSlices-1;uiSliceIndex++)
    {
        dt=g_LATTICE_DATE.at(uiSliceIndex+1)-g_LATTICE_DATE.at(uiSliceIndex);
        if(g_uiProductModelCode<4)
        {
            // condition on whether we are dealing with unrotated or rotated to the diagonal factors
            if(g_uiNumFactors==2)
            {
                // Rotation formulas
                double dMeanReversion1 = hwmr(1,uiSliceIndex);
                double dMeanReversion2 = hwmr(2,uiSliceIndex);
                g_dDiffusion_11.at(uiSliceIndex)=g_dNormalisation_11*(-dMeanReversion1-g_S1_DERIVATIVE_TD.at(uiSliceIndex)/g_S1_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_12.at(uiSliceIndex)=g_dNormalisation_12*(-dMeanReversion2-g_S2_DERIVATIVE_TD.at(uiSliceIndex)/g_S2_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_21.at(uiSliceIndex)=g_dNormalisation_21*(-dMeanReversion1-g_S1_DERIVATIVE_TD.at(uiSliceIndex)/g_S1_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_22.at(uiSliceIndex)=g_dNormalisation_22*(-dMeanReversion2-g_S2_DERIVATIVE_TD.at(uiSliceIndex)/g_S2_TD.at(uiSliceIndex))*dt;
            }
            else if(g_uiNumFactors==3)
            {
                // Rotation formulas
                double dMeanReversion1 = hwmr(1,uiSliceIndex);
                double dMeanReversion2 = hwmr(2,uiSliceIndex);
                double dMeanReversion3 = hwmr(3,uiSliceIndex);
                g_dDiffusion_11.at(uiSliceIndex)=g_d3Normalisation_11*(-dMeanReversion1-g_S1_DERIVATIVE_TD.at(uiSliceIndex)/g_S1_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_12.at(uiSliceIndex)=g_d3Normalisation_12*(-dMeanReversion2-g_S2_DERIVATIVE_TD.at(uiSliceIndex)/g_S2_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_13.at(uiSliceIndex)=g_d3Normalisation_13*(-dMeanReversion3-g_S3_DERIVATIVE_TD.at(uiSliceIndex)/g_S3_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_21.at(uiSliceIndex)=g_d3Normalisation_21*(-dMeanReversion1-g_S1_DERIVATIVE_TD.at(uiSliceIndex)/g_S1_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_22.at(uiSliceIndex)=g_d3Normalisation_22*(-dMeanReversion2-g_S2_DERIVATIVE_TD.at(uiSliceIndex)/g_S2_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_23.at(uiSliceIndex)=g_d3Normalisation_23*(-dMeanReversion3-g_S3_DERIVATIVE_TD.at(uiSliceIndex)/g_S3_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_31.at(uiSliceIndex)=g_d3Normalisation_31*(-dMeanReversion1-g_S1_DERIVATIVE_TD.at(uiSliceIndex)/g_S1_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_32.at(uiSliceIndex)=g_d3Normalisation_32*(-dMeanReversion2-g_S2_DERIVATIVE_TD.at(uiSliceIndex)/g_S2_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_33.at(uiSliceIndex)=g_d3Normalisation_33*(-dMeanReversion3-g_S3_DERIVATIVE_TD.at(uiSliceIndex)/g_S3_TD.at(uiSliceIndex))*dt;
            }
            else
            {
                g_dDiffusion_11.at(uiSliceIndex)=(-hwmr(1,uiSliceIndex)
                                                  -g_S1_DERIVATIVE_TD.at(uiSliceIndex)/g_S1_TD.at(uiSliceIndex)
                                                  //																						-g_dEquilibriumForwardBase_Derivative.at(uiSliceIndex)/g_dEquilibriumForwardBase.at(uiSliceIndex)
                                                 ) * dt;
            }
        }
        else // we are in hybrid mode
        {
            // condition on whether we are dealing with unrotated or rotated to the diagonal factors
            if(g_uiNumFactors==2)
            {
                // Rotation formulas
                double dMeanReversion1 = hwmr(1,uiSliceIndex);
                double dMeanReversion2 = hwmr(2,uiSliceIndex);
                g_dDiffusion_11.at(uiSliceIndex)=g_dNormalisation_11*(-dMeanReversion1-g_S1_DERIVATIVE_TD.at(uiSliceIndex)/g_S1_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_12.at(uiSliceIndex)=g_dNormalisation_12*(-dMeanReversion2-g_S2_DERIVATIVE_TD.at(uiSliceIndex)/g_S2_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_21.at(uiSliceIndex)=g_dNormalisation_21*(-dMeanReversion1-g_S1_DERIVATIVE_TD.at(uiSliceIndex)/g_S1_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_22.at(uiSliceIndex)=g_dNormalisation_22*(-dMeanReversion2-g_S2_DERIVATIVE_TD.at(uiSliceIndex)/g_S2_TD.at(uiSliceIndex))*dt;
            }
            else if(g_uiNumFactors==3)
            {
                // Rotation formulas
                // Drift of the IR processes
                double dMeanReversion1 = hwmr(1,uiSliceIndex);
                double dMeanReversion2 = hwmr(2,uiSliceIndex);
                g_dDiffusion_11.at(uiSliceIndex)=g_d3Normalisation_11*(-dMeanReversion1-g_S1_DERIVATIVE_TD.at(uiSliceIndex)/g_S1_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_12.at(uiSliceIndex)=g_d3Normalisation_12*(-dMeanReversion2-g_S2_DERIVATIVE_TD.at(uiSliceIndex)/g_S2_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_13.at(uiSliceIndex)=g_d3Normalisation_13*(-g_S3_DERIVATIVE_TD.at(uiSliceIndex)/g_S3_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_21.at(uiSliceIndex)=g_d3Normalisation_21*(-dMeanReversion1-g_S1_DERIVATIVE_TD.at(uiSliceIndex)/g_S1_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_22.at(uiSliceIndex)=g_d3Normalisation_22*(-dMeanReversion2-g_S2_DERIVATIVE_TD.at(uiSliceIndex)/g_S2_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_23.at(uiSliceIndex)=g_d3Normalisation_23*(-g_S3_DERIVATIVE_TD.at(uiSliceIndex)/g_S3_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_31.at(uiSliceIndex)=g_d3Normalisation_31*(-dMeanReversion1-g_S1_DERIVATIVE_TD.at(uiSliceIndex)/g_S1_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_32.at(uiSliceIndex)=g_d3Normalisation_32*(-dMeanReversion2-g_S2_DERIVATIVE_TD.at(uiSliceIndex)/g_S2_TD.at(uiSliceIndex))*dt;
                g_dDiffusion_33.at(uiSliceIndex)=g_d3Normalisation_33*(-g_S3_DERIVATIVE_TD.at(uiSliceIndex)/g_S3_TD.at(uiSliceIndex))*dt;
            }
            else
            {
                throw("LabelLatticeLangevinDiffusion:: You should never come here");
            }
        }
    }
} // LabelLatticeLangevinDiffusion(...)


long LatticeFromMeanToState(double mean, double dr)
{
    long iAbsStep = (long)(fabs(mean)/dr + 0.5);
    if (mean>0)
    {
        return iAbsStep;
    }
    else
    {
        return -iAbsStep;
    }
} // LatticeFromMeanToState

double hwsgm(unsigned uiFactorIndex, unsigned uiSliceIndex)
{
    return g_sigma.at(uiFactorIndex-1, uiSliceIndex);
} // hwsgm(...)

double LatticeShortRateForwardStd( double dt, unsigned uiFactorIndex, unsigned uiSliceIndex )
{
    // double a = hwmr(uiFactorIndex,uiSliceIndex);
    double s = hwsgm(uiFactorIndex,uiSliceIndex);
    double dVariance = 0.;
    dVariance = s*s*dt;
    return dVariance;

} // LatticeShortRateForwardStd(...)

DKMaille<double>& LatticeProbabilitySolution1F(double dMean,
        double dVariance,
        double dCurrentRate,
        double dTargetShortRate,
        double dRateStepsNextSlice)
{
    // dMean is the variation of the mean of the X brownian motion
    // dVariance is the variation of the square of the X brownian motion
    static DKMaille<double> prob(3);
    static DKMaille<double> analyticsprob(3);
    double dK = dTargetShortRate / dRateStepsNextSlice;
    double dX = dCurrentRate;
    double dDeltaX = dRateStepsNextSlice;
    double dM = dMean;
    double dV = dVariance;
    if(dVariance<1.e-12)
    {
        prob[2]=1./6.;
        prob[1]=2./3.;
        prob[0]=1./6.;
    }
    else
    {
        prob[2] = 0.5*( dK*(dK-1.) - (dM+dX)*(2.*dK-1.)/dDeltaX + (dV+(dM+dX)*(dM+dX))/(dDeltaX*dDeltaX) );
        prob[1] = - 2. * prob[2] + (dM+dX)/dDeltaX - (dK-1.);
        prob[0] = 1. - prob[1] - prob[2];
    }
    return prob;
} // LatticeProbabilitySolution1F

DKMaille<double>& LatticeProbabilitySolution1FArbitrary(double dMean,
        double dVariance,
        double dCurrentRate,
        double dTargetShortRateDown,
        double dTargetShortRateMiddle,
        double dTargetShortRateUp,
        double dRateStepsNextSlice)
{
    // dMean is the variation of the mean of the X brownian motion
    // dVariance is the variation of the square of the X brownian motion
    static DKMaille<double> prob(3);
    double dX = dCurrentRate;
    double dM = dMean;
    double dV = dVariance;
    double dXWide = dTargetShortRateUp-dTargetShortRateDown;
    double dXUp = dTargetShortRateUp-dTargetShortRateMiddle;
    double dXDown = dTargetShortRateMiddle-dTargetShortRateDown;
    if(dVariance<1.e-12)
    {
        prob[2]=1./6.;
        prob[1]=2./3.;
        prob[0]=1./6.;
    }
    else
    {
        prob[2] = ((dTargetShortRateDown+dTargetShortRateMiddle)*(dM+dX)-dTargetShortRateDown*dTargetShortRateMiddle
                   -dV-(dM+dX)*(dM+dX))/((-dXUp)*(dXWide));
        prob[1] = (-prob[2]*(dXWide)+(dM+dX)-dTargetShortRateDown)/((dXDown));
        prob[0] = 1. - prob[1] - prob[2];
    }
    return prob;
} // LatticeProbabilitySolution1FArbitrary


DKMaille<double>& LatticeFromTargetStateTo1FProbabilityArbitrary(unsigned int uiSliceIndexNext,
        unsigned int uiSliceIndex,
        int iTargetStateDown,
        int iTargetStateMiddle,
        int iTargetStateUp,
        unsigned int uiFactor,
        double dt,
        double dDiffMean,
        double nextDr,
        DKMaille<double> &dShortXRate)
{

    // but unfortunately this is how Mean and Std are defined in the original design
    unsigned int uiFactorIndex = uiFactor + 1 ;
    // calc TargetXRate
    double dTargetXRateDown = LatticeStateToXRate(uiSliceIndexNext, iTargetStateDown, uiFactorIndex);
    double dTargetXRateMiddle = LatticeStateToXRate(uiSliceIndexNext, iTargetStateMiddle, uiFactorIndex);
    double dTargetXRateUp = LatticeStateToXRate(uiSliceIndexNext, iTargetStateUp, uiFactorIndex);

    // calc mean of variation and variance of variation for 1F probability calculation
    double dDiffVar  = LatticeShortRateForwardStd( dt, uiFactorIndex, uiSliceIndex );
    // calculate corresponding 1F probability

    DKMaille<double>& dProbability1F=LatticeProbabilitySolution1FArbitrary(dDiffMean,
                                     dDiffVar,dShortXRate.at(uiFactor),dTargetXRateDown,dTargetXRateMiddle,dTargetXRateUp,nextDr);
    return dProbability1F;
} // LatticeFromTargetStateTo1FProbabilityArbitrary(...)


DKMaille<double>& LatticeProbabilitySolution1FEdge(double dMean,
        double dVariance,
        double dCurrentRate,
        double dTargetShortRate,
        double dRateStepsNextSlice,
        unsigned int uiSwitch)
{

    // Attention this assummes equal spacing between nodes. TO CHANGE !!!
    static DKMaille<double> prob(3);
    double dK = dTargetShortRate / dRateStepsNextSlice;
    double dX = dCurrentRate;
    double dDeltaX = dRateStepsNextSlice;
    double dM = dMean;
    double dV = dVariance;
    // Special treatment for zero dVariance
    if(dVariance<1.e-12)
    {
        if(uiSwitch==0)
        {
            prob[2]=1.;
            prob[1]=0.;
            prob[0]=0.;
        }
        if(uiSwitch==1)
        {
            prob[2]=0.;
            prob[1]=0.;
            prob[0]=1.;
        }
    }
    else
    {
        if(uiSwitch==0)
        {
            prob[0]=0.;
            prob[1]=((dK+1.)*(dK+1)-(dV+(dM+dX)*(dM+dX))/(dDeltaX*dDeltaX))/(1.+2.*dK);
            if(prob[1]<0.) prob[1]=0.;
            prob[2]=1.-prob[1];
        }
        if(uiSwitch==1)
        {
            prob[1]=((dV+(dM+dX)*(dM+dX))/(dDeltaX*dDeltaX)-(dK-1.)*(dK-1))/(1.+2.*dK);
            if(prob[1]<0.) prob[1]=0.;
            prob[0]=1.-prob[1];
            prob[2]=0.;
        }
    }
    return prob;
}
// LatticeProbabilitySolution1FEdge()

DKMaille<double>& LatticeFromTargetStateTo1FProbabilityEdge(unsigned int uiSliceIndexNext,
        unsigned int uiSliceIndex,
        int iTargetState,
        unsigned int uiFactor,
        double dt,
        double dDiffMean,
        double nextDr,
        DKMaille<double> &dShortXRate,
        unsigned int uiSwitch)
{
    // assign factor index from factor (I don't like this
    // but unfortunately this is how Mean and Std are defined in the original design
    unsigned int uiFactorIndex = uiFactor + 1 ;
    // calc TargetXRate
    double dTargetXRate = LatticeStateToXRate(uiSliceIndexNext, iTargetState, uiFactorIndex);
    // calc mean of variation and variance of variation for 1F probability calculation
    double dDiffVar  = LatticeShortRateForwardStd( dt, uiFactorIndex, uiSliceIndex );
    DKMaille<double>& dProbability1F=
        LatticeProbabilitySolution1FEdge(dDiffMean,dDiffVar,dShortXRate.at(uiFactor),dTargetXRate,nextDr,uiSwitch);
    return dProbability1F;

} // LatticeFromTargetStateTo1FProbabilityEdge(...)

void LatticeSetLimits(unsigned int uiSlice)
{
    unsigned int uiNextSliceIndex = uiSlice + 1;

    double dDate = g_LATTICE_DATE.at(uiSlice);
    double dNextDate = g_LATTICE_DATE.at(uiSlice+1);

    g_dStdDevCutOff.at(uiSlice).resize(g_uiNumFactors);
    g_dUpper.at(uiSlice).resize(g_uiNumFactors);
    g_dLower.at(uiSlice).resize(g_uiNumFactors);
    g_dNextSliceLimitUp.at(uiSlice).resize(g_uiNumFactors);
    g_dCurrentSliceLimitUp.at(uiSlice).resize(g_uiNumFactors);
    g_dNextSliceLimitDown.at(uiSlice).resize(g_uiNumFactors);
    g_dCurrentSliceLimitDown.at(uiSlice).resize(g_uiNumFactors);
    g_dCurrentCenter.at(uiSlice).resize(g_uiNumFactors);
    g_dNextCenter.at(uiSlice).resize(g_uiNumFactors);
    g_iLimit.at(uiSlice).resize(g_uiNumFactors);
    g_iCurrentShift.at(uiSlice).resize(g_uiNumFactors);
    g_iNextShift.at(uiSlice).resize(g_uiNumFactors);

    double LimitY1=0.;
    double LimitY2=0.;
    double LimitY3=0;

    LimitY1=g_dLimitY1.at(uiSlice);
    if(g_uiNumFactors>1) LimitY2=g_dLimitY2.at(uiSlice);
    if(g_uiNumFactors>2) LimitY3=g_dLimitY3.at(uiSlice);

    switch(g_uiNumFactors)
    {
    case 3:
        g_dStdDevCutOff.at(uiSlice).at(0)= LimitY1;
        g_dStdDevCutOff.at(uiSlice).at(1)= LimitY2;
        g_dStdDevCutOff.at(uiSlice).at(2)= LimitY3;
        break;
    case 2:
        g_dStdDevCutOff.at(uiSlice).at(0)= LimitY1;
        g_dStdDevCutOff.at(uiSlice).at(1)= LimitY2;
        break;
    case 1:
        g_dStdDevCutOff.at(uiSlice).at(0)= LimitY1;
        break;
    }


    // get short X rates for all factors
    switch(g_uiNumFactors)
    {
    case 3:
        g_dUpper.at(uiSlice).at(2)=g_LATTICE_MAX_STATE3.at(uiSlice)*g_LATTICE_RATE_STEP3.at(uiSlice);
        g_dLower.at(uiSlice).at(2)=g_LATTICE_MIN_STATE3.at(uiSlice)*g_LATTICE_RATE_STEP3.at(uiSlice);
    case 2:
        g_dUpper.at(uiSlice).at(1)=g_LATTICE_MAX_STATE2.at(uiSlice)*g_LATTICE_RATE_STEP2.at(uiSlice);
        g_dLower.at(uiSlice).at(1)=g_LATTICE_MIN_STATE2.at(uiSlice)*g_LATTICE_RATE_STEP2.at(uiSlice);
    case 1:
        g_dUpper.at(uiSlice).at(0)=g_LATTICE_MAX_STATE1.at(uiSlice)*g_LATTICE_RATE_STEP1.at(uiSlice);
        g_dLower.at(uiSlice).at(0)=g_LATTICE_MIN_STATE1.at(uiSlice)*g_LATTICE_RATE_STEP1.at(uiSlice);
        break;
    default:
        throw("ERROR: Only 1,2 or 3 factors permitted");
        break;
    } //switch
    // This is the where we calculate the expectation of the ShortX rate at t+dt for deciding the connection

    g_iLimit.at(uiSlice).at(0)=(int)g_uiI1Limit;
    if(g_uiNumFactors>1) g_iLimit.at(uiSlice).at(1)=(int)g_uiI2Limit;
    if(g_uiNumFactors>2) g_iLimit.at(uiSlice).at(2)=(int)g_uiI3Limit;

    if(g_uiNumFactors==1)
    {
        // IR 1 factor
        g_dCurrentCenter.at(uiSlice).at(0)=0.;
        g_dNextCenter.at(uiSlice).at(0)=0.;
    }
    if(g_uiNumFactors==2)
    {
        // IR or IR quanto is centered at zero on the IR dimension
        g_dCurrentCenter.at(uiSlice).at(0)=0.;
        g_dNextCenter.at(uiSlice).at(0)=0.;
        g_dCurrentCenter.at(uiSlice).at(1)=0.;
        g_dNextCenter.at(uiSlice).at(1)=0.;
    }

    if(g_uiNumFactors==3)
    {
        if(g_uiProductModelCode<4)
        {
            g_dCurrentCenter.at(uiSlice).at(0)=0.;
            g_dNextCenter.at(uiSlice).at(0)=0.;
            g_dCurrentCenter.at(uiSlice).at(1)=0.;
            g_dNextCenter.at(uiSlice).at(1)=0.;
            g_dCurrentCenter.at(uiSlice).at(2)=0.;
            g_dNextCenter.at(uiSlice).at(2)=0.;
        }
        else
        {
            double dFXCenter=log(g_LATTICE_ZCB_PRICE_FOREIGN.at(uiSlice)/
                                 g_LATTICE_ZCB_PRICE.at(uiSlice));
            double FXCenterForward=log(g_LATTICE_ZCB_PRICE_FOREIGN.at(uiSlice+1)/
                                       g_LATTICE_ZCB_PRICE.at(uiSlice+1));

            g_dCurrentCenter.at(uiSlice).at(0)=g_d3Normalisation_13*dFXCenter;
            g_dNextCenter.at(uiSlice).at(0)=g_d3Normalisation_13*FXCenterForward;
            g_dCurrentCenter.at(uiSlice).at(1)=g_d3Normalisation_23*dFXCenter;
            g_dNextCenter.at(uiSlice).at(1)=g_d3Normalisation_23*FXCenterForward;
            g_dCurrentCenter.at(uiSlice).at(2)=g_d3Normalisation_33*dFXCenter;
            g_dNextCenter.at(uiSlice).at(2)=g_d3Normalisation_33*FXCenterForward;
        }
    }

    unsigned uiFactor = 0;
    double dNextDr=0.;
    double dCurrentDr=0.;

    // get the current centers in state terms
    for (uiFactor = 0; uiFactor < g_uiNumFactors; ++uiFactor)
    {
        if(uiFactor==0) dNextDr=g_LATTICE_RATE_STEP1.at(uiNextSliceIndex);
        if(uiFactor==1) dNextDr=g_LATTICE_RATE_STEP2.at(uiNextSliceIndex);
        if(uiFactor==2) dNextDr=g_LATTICE_RATE_STEP3.at(uiNextSliceIndex);
        if(uiFactor==0) dCurrentDr=g_LATTICE_RATE_STEP1.at(uiSlice);
        if(uiFactor==1) dCurrentDr=g_LATTICE_RATE_STEP2.at(uiSlice);
        if(uiFactor==2) dCurrentDr=g_LATTICE_RATE_STEP3.at(uiSlice);
        g_iCurrentShift.at(uiSlice).at(uiFactor)=LatticeFromMeanToState(g_dCurrentCenter.at(uiSlice).at(uiFactor), dCurrentDr);
        g_iNextShift.at(uiSlice).at(uiFactor)=LatticeFromMeanToState(g_dNextCenter.at(uiSlice).at(uiFactor), dNextDr);
    }

    for (uiFactor = 0; uiFactor < g_uiNumFactors; ++uiFactor)
    {
        // set the excess limit
        g_dNextSliceLimitUp.at(uiSlice).at(uiFactor) = g_dNextCenter.at(uiSlice).at(uiFactor) +
                g_dStdDevCutOff.at(uiSlice).at(uiFactor) * g_sigma(uiFactor,uiSlice+1) * sqrt(dNextDate) ;
        g_dCurrentSliceLimitUp.at(uiSlice).at(uiFactor) = g_dCurrentCenter.at(uiSlice).at(uiFactor) +
                g_dStdDevCutOff.at(uiSlice).at(uiFactor)* g_sigma(uiFactor,uiSlice) * sqrt(dDate) ;
        g_dNextSliceLimitDown.at(uiSlice).at(uiFactor) = g_dNextCenter.at(uiSlice).at(uiFactor) -
                g_dStdDevCutOff.at(uiSlice).at(uiFactor) * g_sigma(uiFactor,uiSlice+1) * sqrt(dNextDate) ;
        g_dCurrentSliceLimitDown.at(uiSlice).at(uiFactor) = g_dCurrentCenter.at(uiSlice).at(uiFactor) -
                g_dStdDevCutOff.at(uiSlice).at(uiFactor) * g_sigma(uiFactor,uiSlice) * sqrt(dDate) ;
    }

} // void LatticeSetLimits()


DKMaille<long>& LatticeSearchNearestNodeOnNextSlice(unsigned int uiSlice,
        DKMaille<long> currentNodeStates,
        double dt,
        DKMaille<double> &nextDr,
        DKMaille<double> &dUnconditionalProbability)
{
    static DKMaille<long> iState(g_uiNumFactors);
    iState.resize(g_uiNumFactors);
    unsigned uiSliceIndex = uiSlice;
    double dDate = g_LATTICE_DATE.at(uiSlice);
    double dNextDate = g_LATTICE_DATE.at(uiSlice+1);
    double dLimitDate=g_dTreeFirstExpiry;
    unsigned int uiNodeIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
    double dADP=0.;
    if(g_dADPLimit!=1.0)
    {
        if(g_uiProductModelCode>3) dADP=g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uiNodeIndex);
        if(g_uiProductModelCode<=3) dADP=g_LATTICE_ARROW_DEBREU_PROB.at(uiSlice).at(uiNodeIndex);
    }
    else dADP=0.5;
    double dCutOffDate=g_dTreeCutOffDate;
    static DKMaille<double>     dShortXRate(g_uiNumFactors);
    static DKMaille<double>     mean(g_uiNumFactors);
    static DKMaille<double>     dDrift(g_uiNumFactors);

    dShortXRate.resize(g_uiNumFactors);
    mean.resize(g_uiNumFactors);
    dDrift.resize(g_uiNumFactors);

    // get short X rates for all factors
    switch(g_uiNumFactors)
    {
    case 3:
        dShortXRate[2]=g_LATTICE_X3_RATE.at(uiSlice).at(uiNodeIndex);
    case 2:
        dShortXRate[1]=g_LATTICE_X2_RATE.at(uiSlice).at(uiNodeIndex);
    case 1:
        dShortXRate[0]=g_LATTICE_X1_RATE.at(uiSlice).at(uiNodeIndex);
        break;
    default:
        throw("ERROR: Only 1,2 or 3 factors permitted");
        break;
    } //switch
    // This is the where we calculate the expectation of the ShortX rate at t+dt for deciding the connection

    unsigned int uiFactor=0;
    for (uiFactor = 1; uiFactor <= g_uiNumFactors; ++uiFactor)
    {
        dDrift[uiFactor-1]=LatticeShortRateForwardMeanNew(dShortXRate,dt,uiFactor,uiSliceIndex,uiNodeIndex);
        mean[uiFactor-1]=dDrift[uiFactor-1]+dShortXRate[uiFactor-1];
    }

    // This is where we test against pre-set limits of the lattice we want to build
    // The rule is that the if the Y1 or Y2 rate exceeds (dStdDevCutOff) * standard deviation * sqrt (t)
    // then we attempt to trim the lattice
    // We then check the probabilities

    static DKMaille<double> dProbability1F(g_uiBranching);
    dProbability1F.resize(g_uiBranching);
    //  DKMaille<double> dTargetXRate(g_uiNumFactors);

    bool bSignature=false;
    double dTargetMean=0.;
    unsigned int uiExcess=0;

    // Improve this to be more dynamic (ADP to follow volatility cone in 2-D) DK to do at a later release
    // DK note use Rayleigh method to decide the diffusion cone
    bool bForwardCondition=false;
    bool bEdgeCondition=false;

    bool bTopEdge=false;
    bool bLowEdge=false;

    unsigned int uiSwitch=0;

    for( uiFactor=0; uiFactor < g_uiNumFactors; uiFactor++ )
    {
        bool bEfficiency=false;
        if((g_LATTICE_DATE.at(uiSlice)>g_dTimeBoost)&&(g_uiNumFactors>2)&&dADP<g_dADPLimit) bEfficiency=true;

        // initialise signature and excess flags
        uiExcess=0;

        bForwardCondition=false;
        bEdgeCondition=false;

        // set the condition of excess
        if((mean[uiFactor] >= g_dNextSliceLimitUp.at(uiSlice).at(uiFactor))||(mean[uiFactor] <= g_dNextSliceLimitDown.at(uiSlice).at(uiFactor))) bForwardCondition=true;
        if((dShortXRate[uiFactor]==g_dUpper.at(uiSlice).at(uiFactor))||(dShortXRate[uiFactor]==g_dLower.at(uiSlice).at(uiFactor))) bEdgeCondition=true;

        if((dShortXRate[uiFactor]>= g_dCurrentSliceLimitUp.at(uiSlice).at(uiFactor))
                &&bForwardCondition
                &&bEdgeCondition) uiExcess=1;
        else if((dShortXRate[uiFactor] <= g_dCurrentSliceLimitDown.at(uiSlice).at(uiFactor))
                &&bForwardCondition
                &&bEdgeCondition) uiExcess=1;
        else uiExcess=0;

        if(uiSliceIndex==0) uiExcess=0;

        iState[uiFactor]=LatticeFromMeanToState(mean[uiFactor], nextDr[uiFactor]);

        if(uiExcess==0&&g_uiProductModelCode>3&&bEfficiency)
        {
            if(iState.at(uiFactor)>=g_iLimit.at(uiSlice).at(uiFactor)+g_iNextShift.at(uiSlice).at(uiFactor)
                    ||iState.at(uiFactor)<=-g_iLimit.at(uiSlice).at(uiFactor)+g_iNextShift.at(uiSlice).at(uiFactor))
                uiExcess=1;
        }


        // if target mean exceeds limit then bring the node down and calc 1F probability
        if(uiExcess == 0)
        {
            bSignature=true;
        }
        else
        {
            bSignature=false;
            // set the target state to the immediately lower absolute state
            if((mean[uiFactor] >= g_dNextSliceLimitUp.at(uiSlice).at(uiFactor))
                    ||(dShortXRate[uiFactor]==g_dUpper.at(uiSlice).at(uiFactor))
                    ||(dShortXRate[uiFactor]>= g_dCurrentSliceLimitUp.at(uiSlice).at(uiFactor))
                    ||(iState.at(uiFactor)>=g_iLimit.at(uiSlice).at(uiFactor)+g_iNextShift.at(uiSlice).at(uiFactor)))
            {
                iState[uiFactor]-=1;
                uiSwitch=0;
            }
            if((mean[uiFactor] <= g_dNextSliceLimitDown.at(uiSlice).at(uiFactor))
                    ||(dShortXRate[uiFactor]==g_dLower.at(uiSlice).at(uiFactor))
                    ||(dShortXRate[uiFactor] <= g_dCurrentSliceLimitDown.at(uiSlice).at(uiFactor))
                    ||(iState.at(uiFactor)<=-g_iLimit.at(uiSlice).at(uiFactor)+g_iNextShift.at(uiSlice).at(uiFactor)))
            {
                iState[uiFactor]+=1;
                uiSwitch=1;
            }
            // Set absolute limits on the lattice
            if(bEfficiency)
            {
                if(iState.at(uiFactor)>g_iLimit.at(uiSlice).at(uiFactor)+g_iNextShift.at(uiSlice).at(uiFactor))
                {
                    iState[uiFactor]=g_iLimit.at(uiSlice).at(uiFactor)+g_iNextShift.at(uiSlice).at(uiFactor);
                    uiSwitch=0;
                }
                if(iState.at(uiFactor)<-g_iLimit.at(uiSlice).at(uiFactor)+g_iNextShift.at(uiSlice).at(uiFactor))
                {
                    iState[uiFactor]=-g_iLimit.at(uiSlice).at(uiFactor)+g_iNextShift.at(uiSlice).at(uiFactor);
                    uiSwitch=1;
                }
            }

            // calculate corresponding 1F probability


            dProbability1F=LatticeFromTargetStateTo1FProbabilityArbitrary(uiSliceIndex+1,
                           uiSliceIndex,
                           iState[uiFactor]-1,
                           iState[uiFactor],
                           iState[uiFactor]+1,
                           uiFactor,
                           dt,
                           dDrift[uiFactor],
                           nextDr[uiFactor],
                           dShortXRate);



            // Calculate 1Fprobability vector signature and fill the unconditional probability passenger
            for(unsigned uiBranch = 0; uiBranch < g_uiBranching; uiBranch++)
            {
                // Test on the probabilities.
                // if((bEfficiency==false)&&(g_uiProductModelCode!=5)&&(uiFactor!=2))
                {
                    if(dProbability1F.at(uiBranch)<-0.00)
                        bSignature=true;
                }
                if(bSignature==false) dUnconditionalProbability.at(uiBranch+g_uiBranching*uiFactor)=dProbability1F.at(uiBranch);
            }
            if((bEfficiency))
            {
                // If probs negative it is because the edges are not stiff enough
                // Change probability calculation to fit variance only at the edge
                // No negative probs allowed in the lattice
                if(bSignature)
                {
                    // calculate corresponding 1F probability
                    dProbability1F=LatticeFromTargetStateTo1FProbabilityEdge(uiSliceIndex+1,
                                   uiSliceIndex,
                                   iState[uiFactor],
                                   uiFactor,
                                   dt,
                                   dDrift[uiFactor],
                                   nextDr[uiFactor],
                                   dShortXRate,
                                   uiSwitch);
                    // fill the passenger
                    for(unsigned uiBranch = 0; uiBranch < g_uiBranching; uiBranch++)
                    {
                        dUnconditionalProbability.at(uiBranch+g_uiBranching*uiFactor)=dProbability1F.at(uiBranch);
                    }
                    // deactivate flag
                    bSignature=false;
                }
            }
        }
        // if position in lattice does not exceed limit or if a resulting 1F probability
        // is negative then revert to normal branching
        if(bSignature)
        {
            // calculate the resulting target state as normal
            iState[uiFactor] = LatticeFromMeanToState(mean[uiFactor], nextDr[uiFactor]);
            // calculate corresponding 1F probability


            dProbability1F=LatticeFromTargetStateTo1FProbabilityArbitrary(uiSliceIndex+1,
                           uiSliceIndex,
                           iState[uiFactor]-1,
                           iState[uiFactor],
                           iState[uiFactor]+1,
                           uiFactor,
                           dt,
                           dDrift[uiFactor],
                           nextDr[uiFactor],
                           dShortXRate);

            // fill the unconditional probability passenger
            for(unsigned int uiBranch = 0; uiBranch < g_uiBranching; uiBranch++)
            {
                dUnconditionalProbability.at(uiBranch+g_uiBranching*uiFactor)=dProbability1F.at(uiBranch);
            }
        }
        if(uiFactor==0)
        {
            g_LATTICE_TARGET_1.at(uiSlice).at(uiNodeIndex)=iState.at(0);
            g_LATTICE_UP_1.at(uiSlice).at(uiNodeIndex)=dUnconditionalProbability.at(2);
            g_LATTICE_MID_1.at(uiSlice).at(uiNodeIndex)=dUnconditionalProbability.at(1);
            g_LATTICE_DOWN_1.at(uiSlice).at(uiNodeIndex)=dUnconditionalProbability.at(0);
        }
        if(uiFactor==1)
        {
            g_LATTICE_TARGET_2.at(uiSlice).at(uiNodeIndex)=iState.at(1);
            g_LATTICE_UP_2.at(uiSlice).at(uiNodeIndex)=dUnconditionalProbability.at(5);
            g_LATTICE_MID_2.at(uiSlice).at(uiNodeIndex)=dUnconditionalProbability.at(4);
            g_LATTICE_DOWN_2.at(uiSlice).at(uiNodeIndex)=dUnconditionalProbability.at(3);
        }
        if(uiFactor==2)
        {
            g_LATTICE_TARGET_3.at(uiSlice).at(uiNodeIndex)=iState.at(2);
            g_LATTICE_UP_3.at(uiSlice).at(uiNodeIndex)=dUnconditionalProbability.at(8);
            g_LATTICE_MID_3.at(uiSlice).at(uiNodeIndex)=dUnconditionalProbability.at(7);
            g_LATTICE_DOWN_3.at(uiSlice).at(uiNodeIndex)=dUnconditionalProbability.at(6);
        }
    }// end loop over factors
    return iState;
} // LatticeSearchNearestNodeOnNextSlice(...)

void LatticeResizeSlice(unsigned int uiSlice)
{
    // This sizes the next slice to the maximum possible number of nodes.
    // Note that some of those nodes may not be created
    unsigned int uiSize;
    if(g_uiNumFactors==3) uiSize=g_LATTICE_NUMBER_OF_STATES3.at(uiSlice)
                                     *g_LATTICE_NUMBER_OF_STATES2.at(uiSlice)*g_LATTICE_NUMBER_OF_STATES1.at(uiSlice);
    if(g_uiNumFactors==2) uiSize=g_LATTICE_NUMBER_OF_STATES2.at(uiSlice)*g_LATTICE_NUMBER_OF_STATES1.at(uiSlice);
    if(g_uiNumFactors==1) uiSize=g_LATTICE_NUMBER_OF_STATES1.at(uiSlice);

    unsigned int ui=0;

    switch(g_uiNumFactors)
    {
    case 3:
        {
            g_LATTICE_X3_RATE.at(uiSlice).resize(uiSize);
            g_LATTICE_TARGET_3.at(uiSlice).resize(uiSize);
            g_LATTICE_UP_3.at(uiSlice).resize(uiSize);
            g_LATTICE_MID_3.at(uiSlice).resize(uiSize);
            g_LATTICE_DOWN_3.at(uiSlice).resize(uiSize);
        }
    case 2:
        {
            g_LATTICE_X2_RATE.at(uiSlice).resize(uiSize);
            g_LATTICE_TARGET_2.at(uiSlice).resize(uiSize);
            g_LATTICE_UP_2.at(uiSlice).resize(uiSize);
            g_LATTICE_MID_2.at(uiSlice).resize(uiSize);
            g_LATTICE_DOWN_2.at(uiSlice).resize(uiSize);
        }
    case 1:
        {
            g_LATTICE_X1_RATE.at(uiSlice).resize(uiSize);
            g_LATTICE_TARGET_1.at(uiSlice).resize(uiSize);
            g_LATTICE_UP_1.at(uiSlice).resize(uiSize);
            g_LATTICE_MID_1.at(uiSlice).resize(uiSize);
            g_LATTICE_DOWN_1.at(uiSlice).resize(uiSize);
            g_LATTICE_G_RATE.at(uiSlice).resize(uiSize);
            g_LATTICE_SHORT_RATE.at(uiSlice).resize(uiSize, 0.0);
            g_LATTICE_ARROW_DEBREU.at(uiSlice).resize(uiSize);
            g_LATTICE_ARROW_DEBREU_PROB.at(uiSlice).resize(uiSize);
            g_LATTICE_NODE_EXISTENCE.at(uiSlice).resize(uiSize);

            if(g_uiProductModelCode>3)
            {
                g_LATTICE_G_RATE_FOREIGN.at(uiSlice).resize(uiSize);
                g_LATTICE_SHORT_RATE_FOREIGN.at(uiSlice).resize(uiSize);
                if(g_uiProductModelCode==5)
                {
                    g_LATTICE_G_RATE_FX.at(uiSlice).resize(uiSize);
                    g_LATTICE_SHORT_RATE_FX.at(uiSlice).resize(uiSize);
                }
            }

            if(g_uiNumFactors==1)
            {
                for(unsigned int uj=0;uj<3;uj++)
                {
                    g_LATTICE_GREEN.at(uiSlice,uj).resize(uiSize);
                    g_LATTICE_CONNECTION.at(uiSlice,uj).resize(uiSize);
                }
            }

            if(g_uiNumFactors==2)
            {
                for(unsigned int uj=0;uj<9;uj++)
                {
                    g_LATTICE_GREEN.at(uiSlice,uj).resize(uiSize);
                    g_LATTICE_CONNECTION.at(uiSlice,uj).resize(uiSize);
                }
            }

            if(g_uiNumFactors==3)
            {
                for(unsigned int uj=0;uj<27;uj++)
                {
                    g_LATTICE_GREEN.at(uiSlice,uj).resize(uiSize);
                    g_LATTICE_CONNECTION.at(uiSlice,uj).resize(uiSize);
                }
            }

            for(ui=0;ui<uiSize;ui++)
            {
                g_LATTICE_ARROW_DEBREU.at(uiSlice).at(ui)=0.;
                g_LATTICE_ARROW_DEBREU_PROB.at(uiSlice).at(ui)=0.;
            }
            // Initialize the existence holder
            for(ui=0;ui<uiSize;ui++)
                g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(ui)=false;
        }
    }

} // LatticeResizeSlice

int iIsNoticeDate(unsigned int uiLatticeDate,DKMaille<double> &dNoticeDates)
{
    int iRet=-1;
    for(unsigned int ui=0;ui<dNoticeDates.entries();ui++)
    {
        if(g_LATTICE_DATE.at(uiLatticeDate)==dNoticeDates.at(ui))
        {
            iRet=(int)ui;
        }
    }
    return iRet;
}

void LatticeFreeSlice(unsigned int uiSlice)
{
    // This frees unutilised memory in backwardation
    int iFlag=iIsNoticeDate(uiSlice,g_dNoticeDates);

    iFlag = uiSlice==0 ? iFlag = 0 : iFlag; /// JMP to avoid erasing data at time 0

    switch(g_uiNumFactors)
    {
    case 3:
        {
            if(g_dADPLimit!=1.0&&iFlag==-1) g_LATTICE_X3_RATE.at(uiSlice).clear();
            g_LATTICE_TARGET_3.at(uiSlice).clear();
            g_LATTICE_UP_3.at(uiSlice).clear();
            g_LATTICE_MID_3.at(uiSlice).clear();
            g_LATTICE_DOWN_3.at(uiSlice).clear();
        }
    case 2:
        {
            if(g_dADPLimit!=1.0&&iFlag==-1) g_LATTICE_X2_RATE.at(uiSlice).clear();
            g_LATTICE_TARGET_2.at(uiSlice).clear();
            g_LATTICE_UP_2.at(uiSlice).clear();
            g_LATTICE_MID_2.at(uiSlice).clear();
            g_LATTICE_DOWN_2.at(uiSlice).clear();
        }
    case 1:
        {
            if(g_dADPLimit!=1.0&&iFlag==-1) g_LATTICE_X1_RATE.at(uiSlice).clear();
            g_LATTICE_TARGET_1.at(uiSlice).clear();
            g_LATTICE_UP_1.at(uiSlice).clear();
            g_LATTICE_MID_1.at(uiSlice).clear();
            g_LATTICE_DOWN_1.at(uiSlice).clear();

            if(g_dADPLimit!=1.0&&iFlag==-1&&g_uiProductModelCode>3)
                g_LATTICE_G_RATE.at(uiSlice).clear();

            if(g_dSurvivalCalculation==0.0)
            {
                if(g_dADPLimit!=1.0&&iFlag==-1&&g_uiProductModelCode>3)
                    g_LATTICE_SHORT_RATE.at(uiSlice).clear();
            }

            if(g_dADPLimit!=1.0&&iFlag==-1&&g_uiProductModelCode>3)
                g_LATTICE_ARROW_DEBREU.at(uiSlice).clear();
            if(g_dADPLimit!=1.0&&iFlag==-1&&g_uiProductModelCode>3)
                g_LATTICE_ARROW_DEBREU_PROB.at(uiSlice).clear();
            if(g_dADPLimit!=1.0&&iFlag==-1&&g_uiProductModelCode>3)
                g_LATTICE_G_RATE_FOREIGN.at(uiSlice).clear();
            if(g_dADPLimit!=1.0&&iFlag==-1&&g_uiProductModelCode>3)
                g_LATTICE_SHORT_RATE_FOREIGN.at(uiSlice).clear();
            if(g_dADPLimit!=1.0&&iFlag==-1&&g_uiProductModelCode>3)
                g_LATTICE_G_RATE_FX.at(uiSlice).clear();
            if(g_dADPLimit!=1.0&&iFlag==-1&&g_uiProductModelCode>3)
                g_LATTICE_SHORT_RATE_FX.at(uiSlice).clear();
        }
    }
} // LatticeFreeSlice

void LatticeFreeGreen(unsigned int uiSlice)
{
    if(g_uiNumFactors==1)
    {
        for(unsigned int uj=0;uj<3;uj++)
        {
            g_LATTICE_GREEN.at(uiSlice,uj).clear();
            g_LATTICE_CONNECTION.at(uiSlice,uj).clear();
        }
    }

    if(g_uiNumFactors==2)
    {
        for(unsigned int uj=0;uj<9;uj++)
        {
            g_LATTICE_GREEN.at(uiSlice,uj).clear();
            g_LATTICE_CONNECTION.at(uiSlice,uj).clear();
        }
    }

    if(g_uiNumFactors==3)
    {
        for(unsigned int uj=0;uj<27;uj++)
        {
            g_LATTICE_GREEN.at(uiSlice,uj).clear();
            g_LATTICE_CONNECTION.at(uiSlice,uj).clear();
        }
    }

    int iFlag=iIsNoticeDate(uiSlice,g_dNoticeDates);

    if(g_dADPLimit!=1.0&&iFlag>-1&&g_uiProductModelCode>3)
        g_LATTICE_G_RATE.at(uiSlice).clear();
    if(g_dADPLimit!=1.0&&iFlag>-1&&g_uiProductModelCode>3)
        g_LATTICE_SHORT_RATE.at(uiSlice).clear();
    if(g_dADPLimit!=1.0&&iFlag>-1&&g_uiProductModelCode>3)
        g_LATTICE_ARROW_DEBREU.at(uiSlice).clear();
    if(g_dADPLimit!=1.0&&iFlag>-1&&g_uiProductModelCode>3)
        g_LATTICE_ARROW_DEBREU_PROB.at(uiSlice).clear();
    if(g_dADPLimit!=1.0&&iFlag>-1&&g_uiProductModelCode>3)
        g_LATTICE_G_RATE_FOREIGN.at(uiSlice).clear();
    if(g_dADPLimit!=1.0&&iFlag>-1&&g_uiProductModelCode>3)
        g_LATTICE_SHORT_RATE_FOREIGN.at(uiSlice).clear();
    if(g_dADPLimit!=1.0&&iFlag>-1&&g_uiProductModelCode>3)
        g_LATTICE_G_RATE_FX.at(uiSlice).clear();
    if(g_dADPLimit!=1.0&&iFlag>-1&&g_uiProductModelCode>3)
        g_LATTICE_SHORT_RATE_FX.at(uiSlice).clear();

    if(g_dADPLimit!=1.0&&iFlag>-1&&g_uiProductModelCode==3)
    {
        g_LATTICE_X1_RATE.at(uiSlice).clear();
        g_LATTICE_X2_RATE.at(uiSlice).clear();
        g_LATTICE_X3_RATE.at(uiSlice).clear();
    }

    if(g_dADPLimit!=1.0&&iFlag>-1&&g_uiProductModelCode==2)
    {
        g_LATTICE_X1_RATE.at(uiSlice).clear();
        g_LATTICE_X2_RATE.at(uiSlice).clear();
    }

    if(g_dADPLimit!=1.0&&iFlag>-1&&g_uiProductModelCode==1)
    {
        g_LATTICE_X1_RATE.at(uiSlice).clear();
    }

} // LatticeFreeGreen

double LatticeCurrentDf(unsigned int uiSlice, double dShiftShortGRate)
{
    double df = 0.0;
    static DKMaille<double>      dX(g_uiNumFactors);
    static DKMaille<double>   dXNow(g_uiNumFactors);
    static DKMaille<long>   currentNodeStates(g_uiNumFactors);

    dX.resize(g_uiNumFactors);
    dXNow.resize(g_uiNumFactors);
    currentNodeStates.resize(g_uiNumFactors);

    double dt = g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice);

    switch(g_uiNumFactors)
    {
    case 1:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                currentNodeStates.at(0)=i1;
                unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                {
                    double dAD = g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uiIndex);
                    dX[0] = g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                    dXNow[0] = dX[0] ;
                    double dG = LatticeXToGRate(dXNow,uiSlice);
                    double dShortRate =  LatticeGToShortRate(dShiftShortGRate + dG,uiSlice);
                    double discountFactorForThePeriod = LatticePrimaryDiscount(uiSlice,dShortRate);
                    df += discountFactorForThePeriod * dAD;
                }
            } // while loop over nodes
        }
        break;
    case 2:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                    // Get the states from the next slice that we will need to use in connecting
                    if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                    {
                        double dAD = g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uiIndex);
                        dX[0] = g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                        dX[1] = g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex);
                        dXNow[0] = dX[0] ;
                        dXNow[1] = dX[1] ;
                        double dG = LatticeXToGRate(dXNow,uiSlice);
                        double dShortRate =  LatticeGToShortRate(dShiftShortGRate + dG,uiSlice);
                        double discountFactorForThePeriod = LatticePrimaryDiscount(uiSlice,dShortRate);
                        df += discountFactorForThePeriod * dAD;
                    }
                }
            } // while loop over nodes
        }
        break;
    case 3:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    for(int i3=g_LATTICE_MIN_STATE3.at(uiSlice);i3<=g_LATTICE_MAX_STATE3.at(uiSlice);i3++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                        // Get the states from the next slice that we will need to use in connecting
                        if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                        {
                            double dAD = g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uiIndex);
                            dX[0] = g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                            dX[1] = g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex);
                            dX[2] = g_LATTICE_X3_RATE.at(uiSlice).at(uiIndex);
                            dXNow[0] = dX[0] ;
                            dXNow[1] = dX[1] ;
                            dXNow[2] = dX[2] ;
                            double dG = LatticeXToGRate(dXNow,uiSlice);
                            double dShortRate =  LatticeGToShortRate(dShiftShortGRate + dG,uiSlice);
                            double discountFactorForThePeriod = LatticePrimaryDiscount(uiSlice,dShortRate);
                            df += discountFactorForThePeriod * dAD;
                        }
                    }
                }
            } // while loop over nodes
        }
        break;
    }

    return df;
}// LatticeCurrentDf(...)

double LatticeCurrentDfDerivative(unsigned int uiSlice, double dShiftShortGRate)
{
    double df = 0.0;
    static DKMaille<double>      dX(g_uiNumFactors);
    static DKMaille<double>   dXNow(g_uiNumFactors);
    static DKMaille<long>   currentNodeStates(g_uiNumFactors);

    dX.resize(g_uiNumFactors);
    dXNow.resize(g_uiNumFactors);
    currentNodeStates.resize(g_uiNumFactors);

    double dt = g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice);

    switch(g_uiNumFactors)
    {
    case 1:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                currentNodeStates.at(0)=i1;
                unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                {
                    double dAD = g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uiIndex);
                    dX[0] = g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                    dXNow[0] = dX[0] ;
                    double dG = LatticeXToGRate(dXNow,uiSlice);
                    double dShortRate =  LatticeGToShortRate(dShiftShortGRate + dG,uiSlice);
                    double dShortRateDerivative = LatticeGToShortRateDerivative(dShiftShortGRate + dG,uiSlice);
                    double discountFactorForThePeriod = LatticePrimaryDiscount(uiSlice,dShortRate);
                    df += discountFactorForThePeriod * dAD * (- dShortRateDerivative * g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice));
                }
            } // while loop over nodes
        }
        break;
    case 2:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                    // Get the states from the next slice that we will need to use in connecting
                    if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                    {
                        double dAD = g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uiIndex);
                        dX[0] = g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                        dX[1] = g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex);
                        dXNow[0] = dX[0] ;
                        dXNow[1] = dX[1] ;
                        double dG = LatticeXToGRate(dXNow,uiSlice);
                        double dShortRate =  LatticeGToShortRate(dShiftShortGRate + dG,uiSlice);
                        double dShortRateDerivative = LatticeGToShortRateDerivative(dShiftShortGRate + dG,uiSlice);
                        double discountFactorForThePeriod = LatticePrimaryDiscount(uiSlice,dShortRate);
                        df += discountFactorForThePeriod * dAD * (- dShortRateDerivative * g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice));
                    }
                }
            } // while loop over nodes
        }
        break;
    case 3:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    for(int i3=g_LATTICE_MIN_STATE3.at(uiSlice);i3<=g_LATTICE_MAX_STATE3.at(uiSlice);i3++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                        // Get the states from the next slice that we will need to use in connecting
                        if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                        {
                            double dAD = g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uiIndex);
                            dX[0] = g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                            dX[1] = g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex);
                            dX[2] = g_LATTICE_X3_RATE.at(uiSlice).at(uiIndex);
                            dXNow[0] = dX[0] ;
                            dXNow[1] = dX[1] ;
                            dXNow[2] = dX[2] ;
                            double dG = LatticeXToGRate(dXNow,uiSlice);
                            double dShortRate =  LatticeGToShortRate(dShiftShortGRate + dG,uiSlice);
                            double dShortRateDerivative = LatticeGToShortRateDerivative(dShiftShortGRate + dG,uiSlice);
                            double discountFactorForThePeriod = LatticePrimaryDiscount(uiSlice,dShortRate);
                            df += discountFactorForThePeriod * dAD* (- dShortRateDerivative * g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice));
                        }
                    }
                }
            } // while loop over nodes
        }
        break;
    }

    return df;
}// LatticeCurrentDfDerivative(...)


double LatticeGuessTarget(unsigned uiSlice, double &target, double &errorBar)
{
    // this function calculates the variable target for the shifted G rate
    SetLatticeMarketDf(uiSlice+1);
    double marketDf = g_LATTICE_ZCB_PRICE.at(uiSlice+1);
    double df = LatticeCurrentDf(uiSlice, 0.0);
    double dt = g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice);
    target   = ((log(df)-log(marketDf))/dt)/g_dEquilibriumForwardBase.at(uiSlice);
    // target   = ((log(df)-log(marketDf))/dt);
    // Newton-Raphson to get the shifted G rate
    double dGuess=target;
    // double dGuess=0.;
    int iCounter=0;
    while (iCounter<100)
    {
        iCounter++;
        double dFunction = LatticeCurrentDf(uiSlice,dGuess)-marketDf;
        double dFunctionDerivative = LatticeCurrentDfDerivative(uiSlice,dGuess);
        double dRatio = dFunction/dFunctionDerivative;
        dGuess -= dRatio;
        if ( fabs(dRatio) < 1.e-10)
        {
            iCounter = 100;
        }
    }
    target=dGuess;
    errorBar = 0.2;
    return marketDf;
} // LatticeGuessTarget(...)

double LatticeSolveCalibration(unsigned uiSlice)
{
    double target;
    double errorBar;
    double marketDf = LatticeGuessTarget(uiSlice, target, errorBar);
    g_dShiftShortGRate.at(uiSlice)=target;
    return target;
}// LatticeSolveCalibration(...)

void LatticeInitialiseArrowDebreu(unsigned int uiNextSlice)
{
    unsigned int uiSize;

    if(g_uiNumFactors==3) uiSize=g_LATTICE_NUMBER_OF_STATES3.at(uiNextSlice)
                                     *g_LATTICE_NUMBER_OF_STATES2.at(uiNextSlice)*g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice);

    if(g_uiNumFactors==2) uiSize=g_LATTICE_NUMBER_OF_STATES2.at(uiNextSlice)
                                     *g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice);

    if(g_uiNumFactors==1) uiSize=g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice);

    for(unsigned int ui=0;ui<uiSize;ui++)
    {
        g_LATTICE_ARROW_DEBREU.at(uiNextSlice).at(ui)=0.;
        g_LATTICE_ARROW_DEBREU_PROB.at(uiNextSlice).at(ui)=0.;
    }
} // LatticeInitialiseArrowDebreu()

void LatticeInitialiseArrowDebreuBase(unsigned int uiNextSlice)
{
    unsigned int uiSize;

    if(g_uiNumFactors==3) uiSize=g_LATTICE_NUMBER_OF_STATES3.at(uiNextSlice)
                                     *g_LATTICE_NUMBER_OF_STATES2.at(uiNextSlice)*g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice);

    if(g_uiNumFactors==2) uiSize=g_LATTICE_NUMBER_OF_STATES2.at(uiNextSlice)
                                     *g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice);

    if(g_uiNumFactors==1) uiSize=g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice);

    for(unsigned int ui=0;ui<uiSize;ui++)
    {
        g_LATTICE_ARROW_DEBREU.at(uiNextSlice).at(ui)=0.;
        g_LATTICE_ARROW_DEBREU_PROB.at(uiNextSlice).at(ui)=0.;
    }
} // LatticeInitialiseArrowDebreuBase()

void LatticeInitialiseDiffusionProbability(unsigned int uiNextSlice)
{
    unsigned int uiSize;

    if(g_uiNumFactors==3) uiSize=g_LATTICE_NUMBER_OF_STATES3.at(uiNextSlice)
                                     *g_LATTICE_NUMBER_OF_STATES2.at(uiNextSlice)*g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice);

    if(g_uiNumFactors==2) uiSize=g_LATTICE_NUMBER_OF_STATES2.at(uiNextSlice)
                                     *g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice);

    for(unsigned int ui=0;ui<uiSize;ui++)
    {
        g_LATTICE_ARROW_DEBREU_PROB.at(uiNextSlice).at(ui)=0.;
        g_LATTICE_ARROW_DEBREU.at(uiNextSlice).at(ui)=0.;
    }
} // LatticeInitialiseDiffusionProbability()


void LatticeInitialiseArrowDebreuForeign(unsigned int uiNextSlice)
{
    unsigned int uiSize;

    if(g_uiNumFactors==3) uiSize=g_LATTICE_NUMBER_OF_STATES3.at(uiNextSlice)
                                     *g_LATTICE_NUMBER_OF_STATES2.at(uiNextSlice)*g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice);

    if(g_uiNumFactors==2) uiSize=g_LATTICE_NUMBER_OF_STATES2.at(uiNextSlice)
                                     *g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice);

    if(g_uiNumFactors==1) uiSize=g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice);

} // LatticeInitialiseArrowDebreuForeign()


double LatticeCalcConnectionProbability(unsigned int uiSlice,unsigned int uiIndex,DKMaille<long> &iBranch)
{
    // This is for trinomial
    double dProb;
    double dProb1;
    double dProb2;
    double dProb3;

    switch(g_uiNumFactors)
    {
    case 1:
        {
            if(iBranch.at(0)==-1) dProb1=g_LATTICE_DOWN_1.at(uiSlice).at(uiIndex);
            if(iBranch.at(0)==0)  dProb1=g_LATTICE_MID_1.at(uiSlice).at(uiIndex);
            if(iBranch.at(0)==1)  dProb1=g_LATTICE_UP_1.at(uiSlice).at(uiIndex);
            dProb=dProb1;

            if(dProb1<-1.e-10)
                throw("ERROR: Negative probabilities");

        }
        break;
    case 2:
        {

            if(iBranch.at(0)==-1) dProb1=g_LATTICE_DOWN_1.at(uiSlice).at(uiIndex);
            if(iBranch.at(0)==0)  dProb1=g_LATTICE_MID_1.at(uiSlice).at(uiIndex);
            if(iBranch.at(0)==1)  dProb1=g_LATTICE_UP_1.at(uiSlice).at(uiIndex);

            if(iBranch.at(1)==-1) dProb2=g_LATTICE_DOWN_2.at(uiSlice).at(uiIndex);
            if(iBranch.at(1)==0)  dProb2=g_LATTICE_MID_2.at(uiSlice).at(uiIndex);
            if(iBranch.at(1)==1)  dProb2=g_LATTICE_UP_2.at(uiSlice).at(uiIndex);

            if(dProb1<-1.e-10||dProb2<-1.e-10)
                throw("ERROR: Negative probabilities");

            dProb=dProb1*dProb2;
        }
        break;
    case 3:
        {
            if(iBranch.at(0)==-1) dProb1=g_LATTICE_DOWN_1.at(uiSlice).at(uiIndex);
            if(iBranch.at(0)==0)  dProb1=g_LATTICE_MID_1.at(uiSlice).at(uiIndex);
            if(iBranch.at(0)==1)  dProb1=g_LATTICE_UP_1.at(uiSlice).at(uiIndex);

            if(iBranch.at(1)==-1) dProb2=g_LATTICE_DOWN_2.at(uiSlice).at(uiIndex);
            if(iBranch.at(1)==0)  dProb2=g_LATTICE_MID_2.at(uiSlice).at(uiIndex);
            if(iBranch.at(1)==1)  dProb2=g_LATTICE_UP_2.at(uiSlice).at(uiIndex);

            if(iBranch.at(2)==-1) dProb3=g_LATTICE_DOWN_3.at(uiSlice).at(uiIndex);
            if(iBranch.at(2)==0)  dProb3=g_LATTICE_MID_3.at(uiSlice).at(uiIndex);
            if(iBranch.at(2)==1)  dProb3=g_LATTICE_UP_3.at(uiSlice).at(uiIndex);

            if(dProb1<-1.e-10||dProb2<-1.e-10||dProb3<-1.e-10)
                throw("ERROR: Negative probabilities");

            dProb=dProb1*dProb2*dProb3;
        }
        break;
    }
    return dProb;
} // LatticeCalcConnectionProbability()

void LatticeCalculateArrowDebreu(unsigned int uiSlice)
{
    unsigned int uiPreviousSlice=uiSlice-1;
    double dt = g_LATTICE_PERIOD_TIME_SPAN.at(uiPreviousSlice);
    DKMaille<long> currentNodeStates(g_uiNumFactors);
    DKMaille<long> targetNodeStates(g_uiNumFactors);
    // Do it for trinomial only first then move to mixed quadri- and tri-
    // Run loops over previous slice states
    DKMaille<long> iBranch(g_uiNumFactors);
    switch(g_uiNumFactors)
    {
    case 1:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiPreviousSlice);i1<=g_LATTICE_MAX_STATE1.at(uiPreviousSlice);i1++)
            {
                currentNodeStates.at(0)=i1;
                unsigned int uiIndex=LatticeIndexConvert(uiPreviousSlice,currentNodeStates);
                // Create nodes on next slice based on connection
                if(g_LATTICE_NODE_EXISTENCE.at(uiPreviousSlice).at(uiIndex))
                {
                    double dDiscount=LatticeDiscount(uiPreviousSlice,uiIndex);
                    // Forward Nodes
                    for(int uiBranch=-1;uiBranch<2;uiBranch++)
                    {
                        targetNodeStates.at(0)=g_LATTICE_TARGET_1.at(uiPreviousSlice).at(uiIndex)+uiBranch;
                        iBranch.at(0)=uiBranch;
                        double dConnectionProbability=LatticeCalcConnectionProbability(uiPreviousSlice,uiIndex,iBranch);
                        unsigned int uiTargetNodeIndex=LatticeIndexConvert(uiSlice,targetNodeStates);
                        // ADD the Arrow Debreu contribution to it
                        // Addition is there because of the many potential connections to the SAME node
                        g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uiTargetNodeIndex)+=
                            g_LATTICE_ARROW_DEBREU.at(uiPreviousSlice).at(uiIndex)*dConnectionProbability*dDiscount;
                        g_LATTICE_ARROW_DEBREU_PROB.at(uiSlice).at(uiTargetNodeIndex)+=
                            g_LATTICE_ARROW_DEBREU_PROB.at(uiPreviousSlice).at(uiIndex)*dConnectionProbability;
                    }
                }
            } // while loop over nodes
        }
        break;
    case 2:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiPreviousSlice);i1<=g_LATTICE_MAX_STATE1.at(uiPreviousSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiPreviousSlice);i2<=g_LATTICE_MAX_STATE2.at(uiPreviousSlice);i2++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiIndex=LatticeIndexConvert(uiPreviousSlice,currentNodeStates);
                    // Get the states from the next slice that we will need to use in connecting
                    if(g_LATTICE_NODE_EXISTENCE.at(uiPreviousSlice).at(uiIndex))
                    {
                        double dDiscount=LatticeDiscount(uiPreviousSlice,uiIndex);
                        // Forward Nodes
                        for(int uiBranch=-1;uiBranch<2;uiBranch++)
                        {
                            for(int uiBranch2=-1;uiBranch2<2;uiBranch2++)
                            {
                                targetNodeStates.at(0)=g_LATTICE_TARGET_1.at(uiPreviousSlice).at(uiIndex)+uiBranch;
                                targetNodeStates.at(1)=g_LATTICE_TARGET_2.at(uiPreviousSlice).at(uiIndex)+uiBranch2;
                                iBranch.at(0)=uiBranch;
                                iBranch.at(1)=uiBranch2;
                                double dConnectionProbability=LatticeCalcConnectionProbability(uiPreviousSlice,uiIndex,iBranch);
                                unsigned int uiTargetNodeIndex=LatticeIndexConvert(uiSlice,targetNodeStates);
                                // ADD the Arrow Debreu contribution to it
                                // Addition is there because of the many potential connections to the SAME node
                                g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uiTargetNodeIndex)+=
                                    g_LATTICE_ARROW_DEBREU.at(uiPreviousSlice).at(uiIndex)*dConnectionProbability*dDiscount;
                                g_LATTICE_ARROW_DEBREU_PROB.at(uiSlice).at(uiTargetNodeIndex)+=
                                    g_LATTICE_ARROW_DEBREU_PROB.at(uiPreviousSlice).at(uiIndex)*dConnectionProbability;
                            }
                        } // loop over branches
                    }
                }
            } // while loop over nodes
        }
        break;
    case 3:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiPreviousSlice);i1<=g_LATTICE_MAX_STATE1.at(uiPreviousSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiPreviousSlice);i2<=g_LATTICE_MAX_STATE2.at(uiPreviousSlice);i2++)
                {
                    for(int i3=g_LATTICE_MIN_STATE3.at(uiPreviousSlice);i3<=g_LATTICE_MAX_STATE3.at(uiPreviousSlice);i3++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiIndex=LatticeIndexConvert(uiPreviousSlice,currentNodeStates);
                        // Get the states from the next slice that we will need to use in connecting
                        if(g_LATTICE_NODE_EXISTENCE.at(uiPreviousSlice).at(uiIndex))
                        {
                            double dDiscount=LatticeDiscount(uiPreviousSlice,uiIndex);
                            // Forward Nodes
                            for(int uiBranch=-1;uiBranch<2;uiBranch++)
                            {
                                for(int uiBranch2=-1;uiBranch2<2;uiBranch2++)
                                {
                                    for(int uiBranch3=-1;uiBranch3<2;uiBranch3++)
                                    {
                                        targetNodeStates.at(0)=g_LATTICE_TARGET_1.at(uiPreviousSlice).at(uiIndex)+uiBranch;
                                        targetNodeStates.at(1)=g_LATTICE_TARGET_2.at(uiPreviousSlice).at(uiIndex)+uiBranch2;
                                        targetNodeStates.at(2)=g_LATTICE_TARGET_3.at(uiPreviousSlice).at(uiIndex)+uiBranch3;
                                        iBranch.at(0)=uiBranch;
                                        iBranch.at(1)=uiBranch2;
                                        iBranch.at(2)=uiBranch3;
                                        double dConnectionProbability=LatticeCalcConnectionProbability(uiPreviousSlice,uiIndex,iBranch);
                                        unsigned int uiTargetNodeIndex=LatticeIndexConvert(uiSlice,targetNodeStates);
                                        // ADD the Arrow Debreu contribution to it
                                        // Addition is there because of the many potential connections to the SAME node
                                        g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uiTargetNodeIndex)+=
                                            g_LATTICE_ARROW_DEBREU.at(uiPreviousSlice).at(uiIndex)*dConnectionProbability*dDiscount;
                                        g_LATTICE_ARROW_DEBREU_PROB.at(uiSlice).at(uiTargetNodeIndex)+=
                                            g_LATTICE_ARROW_DEBREU_PROB.at(uiPreviousSlice).at(uiIndex)*dConnectionProbability;
                                    }
                                }
                            } // loop over branches
                        }
                    }
                }
            } // while loop over nodes
        }
        break;
    }
} // LatticeCalculateArrowDebreu()

void LatticeCalculateArrowDebreuBase(unsigned int uiSlice)
{
    unsigned int uiPreviousSlice=uiSlice-1;
    double dt = g_LATTICE_PERIOD_TIME_SPAN.at(uiPreviousSlice);
    DKMaille<long> currentNodeStates(g_uiNumFactors);
    DKMaille<long> targetNodeStates(g_uiNumFactors);
    // Do it for trinomial only first then move to mixed quadri- and tri-
    // Run loops over previous slice states
    DKMaille<long> iBranch(g_uiNumFactors);
    switch(g_uiNumFactors)
    {
    case 1:
        {
            throw("LatticeCalculateArrowDebreuBase:: You should never come here");
        }
        break;
    case 2:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiPreviousSlice);i1<=g_LATTICE_MAX_STATE1.at(uiPreviousSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiPreviousSlice);i2<=g_LATTICE_MAX_STATE2.at(uiPreviousSlice);i2++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiIndex=LatticeIndexConvert(uiPreviousSlice,currentNodeStates);
                    // Get the states from the next slice that we will need to use in connecting
                    if(g_LATTICE_NODE_EXISTENCE.at(uiPreviousSlice).at(uiIndex))
                    {
                        double dDiscount=LatticeDiscount(uiPreviousSlice,uiIndex);
                        // Forward Nodes
                        for(int uiBranch=-1;uiBranch<2;uiBranch++)
                        {
                            for(int uiBranch2=-1;uiBranch2<2;uiBranch2++)
                            {
                                targetNodeStates.at(0)=g_LATTICE_TARGET_1.at(uiPreviousSlice).at(uiIndex)+uiBranch;
                                targetNodeStates.at(1)=g_LATTICE_TARGET_2.at(uiPreviousSlice).at(uiIndex)+uiBranch2;
                                iBranch.at(0)=uiBranch;
                                iBranch.at(1)=uiBranch2;
                                double dConnectionProbability=LatticeCalcConnectionProbability(uiPreviousSlice,uiIndex,iBranch);
                                unsigned int uiTargetNodeIndex=LatticeIndexConvert(uiSlice,targetNodeStates);
                                // ADD the Arrow Debreu contribution to it
                                // Addition is there because of the many potential connections to the SAME node
                                g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uiTargetNodeIndex)+=
                                    g_LATTICE_ARROW_DEBREU.at(uiPreviousSlice).at(uiIndex)*dConnectionProbability*dDiscount;
                                g_LATTICE_ARROW_DEBREU_PROB.at(uiSlice).at(uiTargetNodeIndex)+=
                                    g_LATTICE_ARROW_DEBREU_PROB.at(uiPreviousSlice).at(uiIndex)*dConnectionProbability;
                            }
                        } // loop over branches
                    }
                }
            } // while loop over nodes
        }
        break;
    case 3:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiPreviousSlice);i1<=g_LATTICE_MAX_STATE1.at(uiPreviousSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiPreviousSlice);i2<=g_LATTICE_MAX_STATE2.at(uiPreviousSlice);i2++)
                {
                    for(int i3=g_LATTICE_MIN_STATE3.at(uiPreviousSlice);i3<=g_LATTICE_MAX_STATE3.at(uiPreviousSlice);i3++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiIndex=LatticeIndexConvert(uiPreviousSlice,currentNodeStates);
                        // Get the states from the next slice that we will need to use in connecting
                        if(g_LATTICE_NODE_EXISTENCE.at(uiPreviousSlice).at(uiIndex))
                        {
                            double dDiscount=LatticeDiscount(uiPreviousSlice,uiIndex);
                            // Forward Nodes
                            for(int uiBranch=-1;uiBranch<2;uiBranch++)
                            {
                                for(int uiBranch2=-1;uiBranch2<2;uiBranch2++)
                                {
                                    for(int uiBranch3=-1;uiBranch3<2;uiBranch3++)
                                    {
                                        targetNodeStates.at(0)=g_LATTICE_TARGET_1.at(uiPreviousSlice).at(uiIndex)+uiBranch;
                                        targetNodeStates.at(1)=g_LATTICE_TARGET_2.at(uiPreviousSlice).at(uiIndex)+uiBranch2;
                                        targetNodeStates.at(2)=g_LATTICE_TARGET_3.at(uiPreviousSlice).at(uiIndex)+uiBranch3;
                                        iBranch.at(0)=uiBranch;
                                        iBranch.at(1)=uiBranch2;
                                        iBranch.at(2)=uiBranch3;
                                        double dConnectionProbability=LatticeCalcConnectionProbability(uiPreviousSlice,uiIndex,iBranch);
                                        unsigned int uiTargetNodeIndex=LatticeIndexConvert(uiSlice,targetNodeStates);
                                        // ADD the Arrow Debreu contribution to it
                                        // Addition is there because of the many potential connections to the SAME node
                                        g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uiTargetNodeIndex)+=
                                            g_LATTICE_ARROW_DEBREU.at(uiPreviousSlice).at(uiIndex)*dConnectionProbability*dDiscount;
                                        g_LATTICE_ARROW_DEBREU_PROB.at(uiSlice).at(uiTargetNodeIndex)+=
                                            g_LATTICE_ARROW_DEBREU_PROB.at(uiPreviousSlice).at(uiIndex)*dConnectionProbability;
                                    }
                                }
                            } // loop over branches
                        }
                    }
                }
            } // while loop over nodes
        }
        break;
    }
} // LatticeCalculateArrowDebreuBase()


void LatticeCalculateDiffusionProbability(unsigned int uiSlice)
{
    unsigned int uiPreviousSlice=uiSlice-1;
    double dt = g_LATTICE_PERIOD_TIME_SPAN.at(uiPreviousSlice);
    DKMaille<long> currentNodeStates(g_uiNumFactors);
    DKMaille<long> targetNodeStates(g_uiNumFactors);
    // This is a function that gets called only for hybrids
    // Do it for trinomial only first then move to mixed quadri- and tri-
    // Run loops over previous slice states
    DKMaille<long> iBranch(g_uiNumFactors);
    switch(g_uiNumFactors)
    {
    case 1:
        {
            throw("LatticeCalculateArrowDebreuBase:: You should never come here");
        }
        break;
    case 2:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiPreviousSlice);i1<=g_LATTICE_MAX_STATE1.at(uiPreviousSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiPreviousSlice);i2<=g_LATTICE_MAX_STATE2.at(uiPreviousSlice);i2++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiIndex=LatticeIndexConvert(uiPreviousSlice,currentNodeStates);
                    // Get the states from the next slice that we will need to use in connecting
                    if(g_LATTICE_NODE_EXISTENCE.at(uiPreviousSlice).at(uiIndex))
                    {
                        // double dDiscountConnection=LatticeDiscount(uiPreviousSlice,uiIndex);
                        // Foreign Green's function has to be calculated in foreign cash terms
                        // double dForeignDiscountConnection=LatticeDiscountForeign(uiPreviousSlice,uiIndex)*exp(g_dFromForeignDriftToCash.at(uiPreviousSlice)*dt);
                        for(unsigned int um=0;um<9;um++)
                        {
                            // g_LATTICE_ARROW_DEBREU_PROB.at(uiSlice).at(g_LATTICE_CONNECTION.at(uiPreviousSlice,um).at(uiIndex))+=
                            //  g_LATTICE_GREEN.at(uiPreviousSlice,um).at(uiIndex)*g_LATTICE_ARROW_DEBREU_PROB.at(uiPreviousSlice).at(uiIndex)
                            //  /dDiscountConnection;
                            g_LATTICE_ARROW_DEBREU.at(uiSlice).at(g_LATTICE_CONNECTION.at(uiPreviousSlice,um).at(uiIndex))+=
                                (1./(double)g_dConvert)*(double)g_LATTICE_GREEN.at(uiPreviousSlice,um).at(uiIndex)*g_LATTICE_ARROW_DEBREU.at(uiPreviousSlice).at(uiIndex);
                        }
                    }
                }
            } // while loop over nodes
        }
        break;
    case 3:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiPreviousSlice);i1<=g_LATTICE_MAX_STATE1.at(uiPreviousSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiPreviousSlice);i2<=g_LATTICE_MAX_STATE2.at(uiPreviousSlice);i2++)
                {
                    for(int i3=g_LATTICE_MIN_STATE3.at(uiPreviousSlice);i3<=g_LATTICE_MAX_STATE3.at(uiPreviousSlice);i3++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiIndex=LatticeIndexConvert(uiPreviousSlice,currentNodeStates);
                        // Get the states from the next slice that we will need to use in connecting
                        if(g_LATTICE_NODE_EXISTENCE.at(uiPreviousSlice).at(uiIndex))
                        {
                            // double dDiscountConnection=LatticeDiscount(uiPreviousSlice,uiIndex);
                            // Foreign Green's function has to be calculated in foreign cash terms
                            // double dForeignDiscountConnection=LatticeDiscountForeign(uiPreviousSlice,uiIndex)*exp(g_dFromForeignDriftToCash.at(uiPreviousSlice)*dt);
                            for(unsigned int um=0;um<27;um++)
                            {
                                // g_LATTICE_ARROW_DEBREU_PROB.at(uiSlice).at(g_LATTICE_CONNECTION.at(uiPreviousSlice,um).at(uiIndex))+=
                                //  g_LATTICE_GREEN.at(uiPreviousSlice,um).at(uiIndex)*g_LATTICE_ARROW_DEBREU_PROB.at(uiPreviousSlice).at(uiIndex)
                                //  /dDiscountConnection;
                                g_LATTICE_ARROW_DEBREU.at(uiSlice).at(g_LATTICE_CONNECTION.at(uiPreviousSlice,um).at(uiIndex))+=
                                    (1./(double)g_dConvert)*(double)g_LATTICE_GREEN.at(uiPreviousSlice,um).at(uiIndex)*g_LATTICE_ARROW_DEBREU.at(uiPreviousSlice).at(uiIndex);
                            }
                        }
                    }
                }
            } // while loop over nodes
        }
        break;
    }
} // LatticeCalculateDiffusionProbability()

void LatticeApplyCalibration(unsigned int uiSlice, double dShiftShortGRate)
{
    static DKMaille<double> xRates(g_uiNumFactors);
    xRates.resize(g_uiNumFactors);

    // Apply the shift
    double dt = g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice);
    static DKMaille<long> currentNodeStates(g_uiNumFactors);
    currentNodeStates.resize(g_uiNumFactors);

    switch(g_uiNumFactors)
    {
    case 1:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                currentNodeStates.at(0)=i1;
                unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                {
                    xRates[0] = g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                    double dG = LatticeXToGRate(xRates,uiSlice);
                    double dShortRate =  LatticeGToShortRate(dShiftShortGRate + dG,uiSlice);
                    LabelLatticeNodeWithShortRate(uiSlice,uiIndex,dShortRate);
                    LabelLatticeNodeWithGRate(uiSlice,uiIndex,dG);
                }
            } // while loop over nodes
        }
        break;
    case 2:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                    // Get the states from the next slice that we will need to use in connecting
                    if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                    {
                        xRates[0] = g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                        xRates[1] = g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex);
                        double dG = LatticeXToGRate(xRates,uiSlice);
                        double dShortRate =  LatticeGToShortRate(dShiftShortGRate + dG,uiSlice);
                        LabelLatticeNodeWithShortRate(uiSlice,uiIndex,dShortRate);
                        LabelLatticeNodeWithGRate(uiSlice,uiIndex,dG);
                    }
                }
            } // while loop over nodes
        }
        break;
    case 3:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    for(int i3=g_LATTICE_MIN_STATE3.at(uiSlice);i3<=g_LATTICE_MAX_STATE3.at(uiSlice);i3++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                        // Get the states from the next slice that we will need to use in connecting
                        if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                        {
                            xRates[0] = g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                            xRates[1] = g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex);
                            xRates[2] = g_LATTICE_X3_RATE.at(uiSlice).at(uiIndex);
                            double dG = LatticeXToGRate(xRates,uiSlice);
                            double dShortRate =  LatticeGToShortRate(dShiftShortGRate + dG,uiSlice);
                            LabelLatticeNodeWithShortRate(uiSlice,uiIndex,dShortRate);
                            LabelLatticeNodeWithGRate(uiSlice,uiIndex,dG);
                        }
                    }
                }
            } // while loop over nodes
        }
        break;
    }
} // LatticeApplyCalibration(...)

void LatticeApplyCalibrationForeign(unsigned int uiSlice, double dShiftShortGRateForeign)
{
    static DKMaille<double> xRates(g_uiNumFactors);
    xRates.resize(g_uiNumFactors);

    // Apply the shift
    double dt = g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice);
    static DKMaille<long> currentNodeStates(g_uiNumFactors);
    currentNodeStates.resize(g_uiNumFactors);

    switch(g_uiNumFactors)
    {
    case 1:
        {
            throw("LatticeApplyCalibrationBase:: You should not be here");
        }
        break;
    case 2:
        {
            DKMaille<double> dG(g_uiNumFactors);
            DKMaille<double> dShortRate(g_uiNumFactors);
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                    // Get the states from the next slice that we will need to use in connecting
                    if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                    {
                        xRates[0] = g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                        xRates[1] = g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex);
                        LatticeXToGRateHybrid(xRates, dG, uiSlice);
                        //  add the shift to the base lattice
                        dG.at(1)+=dShiftShortGRateForeign;
                        LatticeGToShortRateFXHybrid(dG, dShortRate,uiSlice);
                        LabelLatticeNodeWithShortRateFXHybrid(uiSlice,uiIndex,dShortRate,2);
                        dG.at(1)-=dShiftShortGRateForeign;
                        LabelLatticeNodeWithGRateHybrid(uiSlice,uiIndex,dG,2);
                    }
                }
            } // while loop over nodes
        }
        break;
    case 3:
        {
            DKMaille<double> dG(g_uiNumFactors);
            DKMaille<double> dShortRate(g_uiNumFactors);
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    for(int i3=g_LATTICE_MIN_STATE3.at(uiSlice);i3<=g_LATTICE_MAX_STATE3.at(uiSlice);i3++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                        // Get the states from the next slice that we will need to use in connecting
                        if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                        {
                            xRates[0] = g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                            xRates[1] = g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex);
                            xRates[2] = g_LATTICE_X3_RATE.at(uiSlice).at(uiIndex);
                            LatticeXToGRateHybrid(xRates, dG, uiSlice);
                            //  add the shift to the base lattice
                            dG.at(1)+=dShiftShortGRateForeign;
                            LatticeGToShortRateFXHybrid(dG, dShortRate,uiSlice);
                            LabelLatticeNodeWithShortRateFXHybrid(uiSlice,uiIndex,dShortRate,2);
                            dG.at(1)-=dShiftShortGRateForeign;
                            LabelLatticeNodeWithGRateHybrid(uiSlice,uiIndex,dG,2);
                            // Fill the G-FX holder but not the FX holder. At this stage both base and foreign have been calculated
                            LabelLatticeNodeWithGRateHybrid(uiSlice,uiIndex,dG,3);
                            // LabelLatticeNodeWithShortRateFXHybrid(uiSlice,uiIndex,dShortRate,3);
                        }
                    }
                }
            } // while loop over nodes
        }
        break;
    }
} // LatticeApplyCalibrationForeign(...)

void LatticeApplyCalibrationBase(unsigned int uiSlice, double dShiftShortGRateBase)
{
    static DKMaille<double> xRates(g_uiNumFactors);
    xRates.resize(g_uiNumFactors);

    // Apply the shift
    double dt = g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice);
    static DKMaille<long> currentNodeStates(g_uiNumFactors);
    currentNodeStates.resize(g_uiNumFactors);

    switch(g_uiNumFactors)
    {
    case 1:
        {
            throw("LatticeApplyCalibrationBase:: You should not be here");
        }
        break;
    case 2:
        {
            DKMaille<double> dG(g_uiNumFactors);
            DKMaille<double> dShortRate(g_uiNumFactors);
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                    // Get the states from the next slice that we will need to use in connecting
                    if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                    {
                        xRates[0] = g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                        xRates[1] = g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex);
                        LatticeXToGRateHybrid(xRates, dG, uiSlice);
                        //  add the shift to the base lattice
                        dG.at(0)+=dShiftShortGRateBase;
                        LatticeGToShortRateFXHybrid(dG, dShortRate,uiSlice);
                        LabelLatticeNodeWithShortRateFXHybrid(uiSlice,uiIndex,dShortRate,1);
                        dG.at(0)-=dShiftShortGRateBase;
                        LabelLatticeNodeWithGRateHybrid(uiSlice,uiIndex,dG,1);
                    }
                }
            } // while loop over nodes
        }
        break;
    case 3:
        {
            DKMaille<double> dG(g_uiNumFactors);
            DKMaille<double> dShortRate(g_uiNumFactors);
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    for(int i3=g_LATTICE_MIN_STATE3.at(uiSlice);i3<=g_LATTICE_MAX_STATE3.at(uiSlice);i3++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                        // Get the states from the next slice that we will need to use in connecting
                        if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                        {
                            xRates[0] = g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                            xRates[1] = g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex);
                            xRates[2] = g_LATTICE_X3_RATE.at(uiSlice).at(uiIndex);
                            LatticeXToGRateHybrid(xRates, dG, uiSlice);
                            //  add the shift to the base lattice
                            dG.at(0)+=dShiftShortGRateBase;
                            LatticeGToShortRateFXHybrid(dG, dShortRate,uiSlice);
                            LabelLatticeNodeWithShortRateFXHybrid(uiSlice,uiIndex,dShortRate,1);
                            dG.at(0)-=dShiftShortGRateBase;
                            LabelLatticeNodeWithGRateHybrid(uiSlice,uiIndex,dG,1);
                        }
                    }
                }
            } // while loop over nodes
        }
        break;
    }
} // LatticeApplyCalibrationBase(...)

void LatticeCalibrateSlice(unsigned int uiSlice)
{
    double dShiftShortGRate = LatticeSolveCalibration(uiSlice);
    LatticeApplyCalibration(uiSlice, dShiftShortGRate);
}// LatticeCalibrateSlice(...)

void LatticeCalibrateSliceBase(unsigned int uiSlice)
{
    double dShiftShortGRateBase=0.;
    // This is where we make the U(3)xU(3)xU(1) projection for computational efficiency
    if(g_dShiftShortGRateBase.entries()!=0)
    {
        dShiftShortGRateBase=g_dShiftShortGRateBase.at(uiSlice);
    }
    else
    {
        throw("You should not come here ... !");
    }
    LatticeApplyCalibrationBase(uiSlice, dShiftShortGRateBase);
}// LatticeCalibrateSliceBase(...)

void LatticeCalibrateSliceForeign(unsigned int uiSlice)
{
    double dShiftShortGRateForeign=0.;
    // This is where we make the U(3)xU(3)xU(1) projection for computational efficiency
    if(g_dShiftShortGRateForeign.entries()!=0)
    {
        dShiftShortGRateForeign=g_dShiftShortGRateForeign.at(uiSlice);
    }
    else
    {
        throw("You should not come here ... !");
    }
    LatticeApplyCalibrationForeign(uiSlice, dShiftShortGRateForeign);
}// LatticeCalibrateForeign(...)


// Iterate over the nodes in the slice setting the X rate so that we can use
// it subsequently in the smoothing analytics.
// Do this only in the case of single currency IR two or three-factor models
void LatticeSetXRatesInSlice(unsigned int uiSlice)
{
    // X=Y in one factor
    // if ((1 == g_uiNumFactors))
    //  return;

    switch(g_uiNumFactors)
    {
    case 1:
        {
            double dXX1;
            DKMaille<long> currentNodeStates(g_uiNumFactors);
            DKMaille<double> dYRate(1);
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                currentNodeStates.at(0)=i1;
                unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                // Get the states from the next slice that we will need to use in connecting
                if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                {
                    dYRate.at(0)=g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                    dXX1=dYRate.at(0)*g_S1_TD.at(uiSlice)/g_S1;
                    g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex)=dXX1;
                }
            } // while loop over nodes
        }
        break;
    case 2:
        {
            double dXX1, dXX2;
            DKMaille<long> currentNodeStates(g_uiNumFactors);
            DKMaille<double> dYRate(2);
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                    // Get the states from the next slice that we will need to use in connecting
                    if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                    {

                        dYRate.at(0)=g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                        dYRate.at(1)=g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex);
                        dXX1=LatticeFromYRateToXRate(dYRate,1)*g_S1_TD.at(uiSlice)/g_S1;
                        dXX2=LatticeFromYRateToXRate(dYRate,2)*g_S2_TD.at(uiSlice)/g_S2;
                        g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex)=dXX1;
                        g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex)=dXX2;
                    }
                }
            } // while loop over nodes
        }
        break;
    case 3:
        {
            double dXX1, dXX2, dXX3;
            DKMaille<long> currentNodeStates(g_uiNumFactors);
            DKMaille<double> dYRate(3);
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    for(int i3=g_LATTICE_MIN_STATE3.at(uiSlice);i3<=g_LATTICE_MAX_STATE3.at(uiSlice);i3++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                        // Get the states from the next slice that we will need to use in connecting
                        if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                        {
                            dYRate.at(0)=g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex);
                            dYRate.at(1)=g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex);
                            dYRate.at(2)=g_LATTICE_X3_RATE.at(uiSlice).at(uiIndex);
                            dXX1=LatticeFromYRateToXRate(dYRate,1)*g_S1_TD.at(uiSlice)/g_S1;
                            dXX2=LatticeFromYRateToXRate(dYRate,2)*g_S2_TD.at(uiSlice)/g_S2;
                            dXX3=LatticeFromYRateToXRate(dYRate,3)*g_S3_TD.at(uiSlice)/g_S3;
                            g_LATTICE_X1_RATE.at(uiSlice).at(uiIndex)=dXX1;
                            g_LATTICE_X2_RATE.at(uiSlice).at(uiIndex)=dXX2;
                            g_LATTICE_X3_RATE.at(uiSlice).at(uiIndex)=dXX3;
                        }
                    }
                }
            }// while loop over nodes
        }
        break;
    }

} // LatticeSetXRatesInSlice()

void LatticeSetGreen1F(unsigned int uiSlice,int uiBranch,
                       unsigned int uiBranchIndex,unsigned int uiIndex,
                       double dDiscount)
{
    double dProb1=0.;

    if(uiBranch==-1) dProb1=(g_LATTICE_DOWN_1.at(uiSlice).at(uiIndex));
    if(uiBranch==0) dProb1=(g_LATTICE_MID_1.at(uiSlice).at(uiIndex));
    if(uiBranch==1) dProb1=(g_LATTICE_UP_1.at(uiSlice).at(uiIndex));
    if(dProb1<-1.e-10)
        throw("ERROR: Negative probabilities");
    g_LATTICE_GREEN.at(uiSlice,uiBranchIndex).at(uiIndex)=dProb1*dDiscount*g_dConvert;

} // LatticeSetGreen1F

void LatticeSetGreen2F(unsigned int uiSlice,int uiBranch,int uiBranch2,
                       unsigned int uiBranchIndex,unsigned int uiIndex,
                       double dDiscount)
{
    double dProb1=0.;
    double dProb2=0.;

    if(uiBranch==-1) dProb1=(g_LATTICE_DOWN_1.at(uiSlice).at(uiIndex));
    if(uiBranch==0) dProb1=(g_LATTICE_MID_1.at(uiSlice).at(uiIndex));
    if(uiBranch==1) dProb1=(g_LATTICE_UP_1.at(uiSlice).at(uiIndex));

    if(uiBranch2==-1) dProb2=(g_LATTICE_DOWN_2.at(uiSlice).at(uiIndex));
    if(uiBranch2==0) dProb2=(g_LATTICE_MID_2.at(uiSlice).at(uiIndex));
    if(uiBranch2==1) dProb2=(g_LATTICE_UP_2.at(uiSlice).at(uiIndex));
    if(dProb1<-1.e-10||dProb2<-1.e-10)
        throw("ERROR: Negative probabilities");
    g_LATTICE_GREEN.at(uiSlice,uiBranchIndex).at(uiIndex)=dProb1*dProb2*dDiscount*g_dConvert;

} // LatticeSetGreen2F

void LatticeSetGreen3F(unsigned int uiSlice,int uiBranch,int uiBranch2,int uiBranch3,
                       unsigned int uiBranchIndex,unsigned int uiIndex,
                       double dDiscount)
{
    double dProb1=0.;
    double dProb2=0.;
    double dProb3=0.;

    if(uiBranch==-1) dProb1=(g_LATTICE_DOWN_1.at(uiSlice).at(uiIndex));
    if(uiBranch==0) dProb1=(g_LATTICE_MID_1.at(uiSlice).at(uiIndex));
    if(uiBranch==1) dProb1=(g_LATTICE_UP_1.at(uiSlice).at(uiIndex));

    if(uiBranch2==-1) dProb2=(g_LATTICE_DOWN_2.at(uiSlice).at(uiIndex));
    if(uiBranch2==0) dProb2=(g_LATTICE_MID_2.at(uiSlice).at(uiIndex));
    if(uiBranch2==1) dProb2=(g_LATTICE_UP_2.at(uiSlice).at(uiIndex));

    if(uiBranch3==-1) dProb3=(g_LATTICE_DOWN_3.at(uiSlice).at(uiIndex));
    if(uiBranch3==0) dProb3=(g_LATTICE_MID_3.at(uiSlice).at(uiIndex));
    if(uiBranch3==1) dProb3=(g_LATTICE_UP_3.at(uiSlice).at(uiIndex));
    if(dProb1<-1.e-10||dProb2<-1.e-10||dProb3<-1.e-10)
        throw("ERROR: Negative probabilities");
    g_LATTICE_GREEN.at(uiSlice,uiBranchIndex).at(uiIndex)=dProb1*dProb2*dProb3*dDiscount*g_dConvert;

} // LatticeSetGreen3F

// Vladimir's approximation
double dLatticeDecayFunctional(double dMeanReversion,unsigned int uiSlice,double dPaymentDate,unsigned int uiBaseOrForeign)
{
    double dFunctional=0.;
    double dDiscount1=0.;
    double dDiscount2=0.;
    double dForward=0.;
    double dObservationDate=g_LATTICE_DATE.at(uiSlice);
    if(g_uiIsDK==1)
    {
        unsigned int uiSlices=4;
        double dLocalRate=0.;
        if(uiBaseOrForeign==0) dLocalRate=g_dQParameter_Base*(g_dEquilibriumForwardShortBase.at(uiSlice)-g_dEquilibriumForwardBase.at(uiSlice))
                                              +g_dEquilibriumForwardBase.at(uiSlice);
        else dLocalRate=g_dQParameter_Foreign*(g_dEquilibriumForwardShortForeign.at(uiSlice)-g_dEquilibriumForwardForeign.at(uiSlice))
                            +g_dEquilibriumForwardForeign.at(uiSlice);

        for(unsigned int ui=0;ui<uiSlices;ui++)
        {
            double dt=(dPaymentDate-dObservationDate)/(double)uiSlices;
            double dDate1=dObservationDate+(double)ui*dt;
            double dDate2=dObservationDate+(double)(ui+1)*dt;
            if(uiBaseOrForeign==0)
            {
                dDiscount1=BaseDiscountInterpolate(dDate1);
                dDiscount2=BaseDiscountInterpolate(dDate2);
            }
            if(uiBaseOrForeign==1)
            {
                dDiscount1=ForeignDiscountInterpolate(dDate1);
                dDiscount2=ForeignDiscountInterpolate(dDate2);
            }
            dForward=(dDiscount1/dDiscount2-1.)/dt;
            double dK=0.;
            if(dForward<0.002) dK=0.002;
            else dK=dForward;

            double dQParameter;
            if(uiBaseOrForeign==0) dQParameter=g_dQParameter_Base;
            if(uiBaseOrForeign==1) dQParameter=g_dQParameter_Foreign;


            dFunctional+=(dQParameter*(dForward-dK)+dK)*(exp(-dMeanReversion*(dDate1-dObservationDate))-exp(-dMeanReversion*(dDate2-dObservationDate)));
        }
        dFunctional/=(dMeanReversion*dLocalRate);
    }
    else
    {
        dFunctional=(1.-exp(-dMeanReversion*(dPaymentDate-dObservationDate)))/dMeanReversion;
    }
    return dFunctional;
} // dLatticeDecayFunctional

double CalcGreensFunctionNormalisation(double dDF,
                                       unsigned int uiBaseForeign,
                                       double dObservationDate,
                                       double dBondMaturity,
                                       unsigned int uiSlice)
{

    double dNormalisation=0.;
    unsigned int uiiNumberOfStates;
    switch(g_uiNumFactors)
    {
    case 1:
        uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uiSlice);
        break;
    case 2:
        uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uiSlice)*g_LATTICE_NUMBER_OF_STATES2.at(uiSlice);
        break;
    case 3:
        uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uiSlice)*g_LATTICE_NUMBER_OF_STATES2.at(uiSlice)*g_LATTICE_NUMBER_OF_STATES3.at(uiSlice);
        break;
    }

    // Cases where the fixed income diffusion is restricted to one-factor
    // That incorporates hybrid products
    if(g_uiProductModelCode==1||g_uiProductModelCode>3)
    {
        double dDecay=0.;
        if(uiBaseForeign==0)
        {
            if(g_uiIsDK==0)
            {
                dDecay=(1.-exp(-g_MR1*(dBondMaturity-dObservationDate)))/g_MR1;
            }
            else
            {
                // Vladimir's Approximation
                dDecay=dLatticeDecayFunctional(g_MR1,uiSlice,dBondMaturity,0);
            }
        }
        if(uiBaseForeign==1)
        {
            if(g_uiIsDK==0)
            {
                dDecay=(1.-exp(-g_MR2*(dBondMaturity-dObservationDate)))/g_MR2;
            }
            else
            {
                // Vladimir's Approximation
                dDecay=dLatticeDecayFunctional(g_MR2,uiSlice,dBondMaturity,1);
            }
        }
        double dSpot=g_dSpotFXTilde;
        for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
        {
            if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uk)==true)
            {
                if(uiBaseForeign==0) dNormalisation+=g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uk)*exp(-dDecay*g_LATTICE_SHORT_RATE.at(uiSlice).at(uk));
                // if(uiBaseForeign==0) dNormalisation+=g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uk)*exp(-dDecay*g_LATTICE_X1_RATE.at(uiSlice).at(uk));
                if(uiBaseForeign==1) dNormalisation+=(g_LATTICE_SHORT_RATE_FX.at(uiSlice).at(uk)/dSpot)
                                                         *g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uk)*exp(-dDecay*g_LATTICE_SHORT_RATE_FOREIGN.at(uiSlice).at(uk));
            }
        }
        dNormalisation=dDF/dNormalisation;
    }

    // Cases where the fixed income diffusion is two- or three- factor in single currency
    if(g_uiProductModelCode==2||g_uiProductModelCode==3)
    {
        if(g_uiNumFactors==2)
        {
            double dDecay1=0.;
            double dDecay2=0.;
            dDecay1=(1.-exp(-g_MR1*(dBondMaturity-dObservationDate)))/g_MR1;
            dDecay2=(1.-exp(-g_MR2*(dBondMaturity-dObservationDate)))/g_MR2;
            for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
            {
                if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uk)==true)
                {
                    dNormalisation+=g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uk)*exp(-dDecay1*g_LATTICE_X1_RATE.at(uiSlice).at(uk)
                                    -dDecay2*g_LATTICE_X2_RATE.at(uiSlice).at(uk));
                }
            }
            dNormalisation=dDF/dNormalisation;
        }
        if(g_uiNumFactors==3)
        {
            double dDecay1=0.;
            double dDecay2=0.;
            double dDecay3=0.;
            dDecay1=(1.-exp(-g_MR1*(dBondMaturity-dObservationDate)))/g_MR1;
            dDecay2=(1.-exp(-g_MR2*(dBondMaturity-dObservationDate)))/g_MR2;
            dDecay3=(1.-exp(-g_MR3*(dBondMaturity-dObservationDate)))/g_MR3;
            for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
            {
                if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uk)==true)
                {
                    dNormalisation+=g_LATTICE_ARROW_DEBREU.at(uiSlice).at(uk)*exp(-dDecay1*g_LATTICE_X1_RATE.at(uiSlice).at(uk)
                                    -dDecay2*g_LATTICE_X2_RATE.at(uiSlice).at(uk)
                                    -dDecay3*g_LATTICE_X3_RATE.at(uiSlice).at(uk));
                }
            }
            dNormalisation=dDF/dNormalisation;
        }
    }
    return dNormalisation;
}

void CalcLatticeStringTension(unsigned int uiSlice)
{
    // Calculate Mean String Tension
    double dNorm=CalcGreensFunctionNormalisation(g_LATTICE_ZCB_PRICE.at(uiSlice+1),
                 0,
                 g_LATTICE_DATE.at(uiSlice),
                 g_LATTICE_DATE.at(uiSlice+1),uiSlice);
    g_LATTICE_NORM.at(uiSlice)=dNorm;
    g_LATTICE_DECAY.at(uiSlice)=(1.-exp(-g_MR1*(g_LATTICE_DATE.at(uiSlice+1)-g_LATTICE_DATE.at(uiSlice))))/g_MR1;
}


void CalcLatticeForeignStringTension(unsigned int uiSlice)
{
    // Calculate Mean String Tension
    double dNormForeign=CalcGreensFunctionNormalisation(g_LATTICE_ZCB_PRICE_FOREIGN.at(uiSlice+1),
                        1,
                        g_LATTICE_DATE.at(uiSlice),
                        g_LATTICE_DATE.at(uiSlice+1),uiSlice);
    g_LATTICE_NORM_FOREIGN.at(uiSlice)=dNormForeign;
    g_LATTICE_DECAY_FOREIGN.at(uiSlice)=(1.-exp(-g_MR2*(g_LATTICE_DATE.at(uiSlice+1)-g_LATTICE_DATE.at(uiSlice))))/g_MR2;
}



void LatticeForwardConstructionOfSlices()
{
    // Forward construction of the driftless lattice

    for (unsigned uiSlice = 0 ; uiSlice < g_uiNumSlices-1 ; uiSlice++)
    {
        // DO CALIBRATION HERE CONNECTIONS ARE NOT NEEDED FOR THIS
        SetLatticeMarketDf(uiSlice+1);
        if(g_uiProductModelCode>3) SetLatticeMarketDfForeign(uiSlice+1);

        // calibrate slice is not persisting in the string model.
        // Drop later
        if(g_uiProductModelCode<4)
            LatticeCalibrateSlice(uiSlice);
        if(g_uiProductModelCode>3)
        {
            LatticeCalibrateSliceBase(uiSlice);
            LatticeCalibrateSliceForeign(uiSlice);
        }

        // Get the dt from Now to Next
        double dt = g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice);

        if(g_bStringModel)
        {
            if(g_uiProductModelCode==1)
                CalcLatticeStringTension(uiSlice);
            if((g_uiProductModelCode==4||g_uiProductModelCode==5) )
            {
                CalcLatticeStringTension(uiSlice);
                CalcLatticeForeignStringTension(uiSlice);
            }
        }

        // Set the shortRate step for the NEXT slice
        DKMaille<double> dr(SetLatticeRateStep(uiSlice+1, dt));

        // Create limits for the lattice
        LatticeSetLimits(uiSlice);

        // Initialise the calc of the min and max for the NEXT slice
        LatticeInitMaxMinX(uiSlice+1);

        // Establish connections from each node on the current slice
        DKMaille<double> dUnconditionalProbability(g_uiBranching*g_uiNumFactors);
        DKMaille<long> currentNodeStates(g_uiNumFactors);

        switch(g_uiNumFactors)
        {
        case 1:
            {
                for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
                {
                    currentNodeStates.at(0)=i1;
                    unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                    // Get the states from the next slice that we will need to use in connecting
                    if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                    {
                        DKMaille<long> &centerStates=LatticeSearchNearestNodeOnNextSlice(uiSlice,
                                                     currentNodeStates,dt,dr,dUnconditionalProbability);
                        LatticeSetMaxMinX(centerStates,uiSlice+1);
                    }
                } // while loop over nodes
            }
            break;
        case 2:
            {
                for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
                {
                    for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                        // Get the states from the next slice that we will need to use in connecting
                        if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                        {
                            DKMaille<long> &centerStates=LatticeSearchNearestNodeOnNextSlice(uiSlice,
                                                         currentNodeStates,dt,dr,dUnconditionalProbability);
                            LatticeSetMaxMinX(centerStates,uiSlice+1);
                        }
                    }
                } // while loop over nodes
            }
            break;
        case 3:
            {
                for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
                {
                    for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                    {
                        for(int i3=g_LATTICE_MIN_STATE3.at(uiSlice);i3<=g_LATTICE_MAX_STATE3.at(uiSlice);i3++)
                        {
                            currentNodeStates.at(0)=i1;
                            currentNodeStates.at(1)=i2;
                            currentNodeStates.at(2)=i3;
                            unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                            // Get the states from the next slice that we will need to use in connecting
                            if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                            {
                                DKMaille<long> &centerStates=LatticeSearchNearestNodeOnNextSlice(uiSlice,
                                                             currentNodeStates,dt,dr,dUnconditionalProbability);
                                LatticeSetMaxMinX(centerStates,uiSlice+1);
                            }
                        }
                    }
                } // while loop over nodes
            }
            break;
        }

        //Now resize the next slice to the maximum possible number of nodes
        LatticeResizeSlice(uiSlice+1);
        // Based on the connection targets now establish the existence of nodes on the NEXT slice
        // We do it here for trinomial initially
        DKMaille<long> targetNodeStates(pow(3.,static_cast<double>(g_uiNumFactors)));
        switch(g_uiNumFactors)
        {
        case 1:
            {
                for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
                {
                    currentNodeStates.at(0)=i1;
                    unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                    double dDiscount=LatticeDiscount(uiSlice,uiIndex);
                    // Create nodes on next slice based on connection
                    if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                    {
                        // Forward Node Creation
                        for(int uiBranch=-1;uiBranch<2;uiBranch++)
                        {
                            targetNodeStates.at(0)=g_LATTICE_TARGET_1.at(uiSlice).at(uiIndex)+uiBranch;
                            unsigned int uiTargetNodeIndex=LatticeIndexConvert(uiSlice+1,targetNodeStates);
                            // Create connections and green functions
                            unsigned int uiBranchIndex=(unsigned int)((uiBranch+1));
                            g_LATTICE_CONNECTION.at(uiSlice,uiBranchIndex).at(uiIndex)=uiTargetNodeIndex;
                            LatticeSetGreen1F(uiSlice,uiBranch,uiBranchIndex,uiIndex,dDiscount);
                            // Create the Node if it does not exist already
                            if(g_LATTICE_NODE_EXISTENCE.at(uiSlice+1).at(uiTargetNodeIndex)==false)
                                PrepareTheLatticeNode(uiSlice+1, targetNodeStates);
                        }
                    }
                } // while loop over nodes
            }
            break;
        case 2:
            {
                for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
                {
                    for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                        double dDiscount=LatticeDiscount(uiSlice,uiIndex);
                        // Get the states from the next slice that we will need to use in connecting
                        if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                        {
                            double dProbabilityHandle=0.;
                            // Forward Node Creation
                            for(int uiBranch=-1;uiBranch<2;uiBranch++)
                            {
                                for(int uiBranch2=-1;uiBranch2<2;uiBranch2++)
                                {
                                    targetNodeStates.at(0)=g_LATTICE_TARGET_1.at(uiSlice).at(uiIndex)+uiBranch;
                                    targetNodeStates.at(1)=g_LATTICE_TARGET_2.at(uiSlice).at(uiIndex)+uiBranch2;
                                    unsigned int uiTargetNodeIndex=LatticeIndexConvert(uiSlice+1,targetNodeStates);
                                    // Create connections and green functions for cases of 3F model
                                    unsigned int uiBranchIndex=(unsigned int)(3*(uiBranch+1)+(uiBranch2+1));
                                    g_LATTICE_CONNECTION.at(uiSlice,uiBranchIndex).at(uiIndex)=uiTargetNodeIndex;
                                    LatticeSetGreen2F(uiSlice,uiBranch,uiBranch2,uiBranchIndex,uiIndex,dDiscount);
                                    // Create the Node if it does not exist already
                                    // Create the Node if and only if the probability is not zero
                                    if(g_uiProductModelCode==4) dProbabilityHandle=(1./(double)g_dConvert)*(double)g_LATTICE_GREEN.at(uiSlice,uiBranchIndex).at(uiIndex);
                                    if(g_uiProductModelCode==2) dProbabilityHandle=1.;
                                    if((g_LATTICE_NODE_EXISTENCE.at(uiSlice+1).at(uiTargetNodeIndex)==false)&&
                                            (dProbabilityHandle>1.e-9))
                                    {
                                        PrepareTheLatticeNode(uiSlice+1, targetNodeStates);
                                    }
                                }
                            } // loop over branches
                        }
                    }
                } // while loop over nodes
            }
            break;
        case 3:
            {
                for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
                {
                    for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                    {
                        for(int i3=g_LATTICE_MIN_STATE3.at(uiSlice);i3<=g_LATTICE_MAX_STATE3.at(uiSlice);i3++)
                        {
                            currentNodeStates.at(0)=i1;
                            currentNodeStates.at(1)=i2;
                            currentNodeStates.at(2)=i3;
                            unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                            double dDiscount=LatticeDiscount(uiSlice,uiIndex);
                            // Get the states from the next slice that we will need to use in connecting
                            if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                            {
                                double dProbabilityHandle=0.;
                                // Forward Node Creation
                                for(int uiBranch=-1;uiBranch<2;uiBranch++)
                                {
                                    for(int uiBranch2=-1;uiBranch2<2;uiBranch2++)
                                    {
                                        for(int uiBranch3=-1;uiBranch3<2;uiBranch3++)
                                        {
                                            targetNodeStates.at(0)=g_LATTICE_TARGET_1.at(uiSlice).at(uiIndex)+uiBranch;
                                            targetNodeStates.at(1)=g_LATTICE_TARGET_2.at(uiSlice).at(uiIndex)+uiBranch2;
                                            targetNodeStates.at(2)=g_LATTICE_TARGET_3.at(uiSlice).at(uiIndex)+uiBranch3;
                                            unsigned int uiTargetNodeIndex=LatticeIndexConvert(uiSlice+1,targetNodeStates);
                                            // Create connections and green functions
                                            unsigned int uiBranchIndex=(unsigned int)(3*3*(uiBranch+1)+3*(uiBranch2+1)+(uiBranch3+1));
                                            g_LATTICE_CONNECTION.at(uiSlice,uiBranchIndex).at(uiIndex)=uiTargetNodeIndex;
                                            LatticeSetGreen3F(uiSlice,uiBranch,uiBranch2,uiBranch3,uiBranchIndex,uiIndex,dDiscount);
                                            // Create the Node if and only if the probability is not zero
                                            if(g_uiProductModelCode==5) dProbabilityHandle=(1./(double)g_dConvert)*(double)g_LATTICE_GREEN.at(uiSlice,uiBranchIndex).at(uiIndex);
                                            if(g_uiProductModelCode==3) dProbabilityHandle=1.;
                                            if((g_LATTICE_NODE_EXISTENCE.at(uiSlice+1).at(uiTargetNodeIndex)==false)&&
                                                    (dProbabilityHandle>1.e-9))
                                            {
                                                PrepareTheLatticeNode(uiSlice+1, targetNodeStates);
                                            }
                                        }
                                    }
                                } // loop over branches
                            }
                        }
                    }
                } // while loop over nodes
            }
            break;
        }
        // Now Fill Arrow-Debreu in next slice
        // First initalise to zero
        unsigned int uiNextSlice=uiSlice+1;
        // Single currency fixed income
        // One cannot escape from building Green's function by multi-dimensional induction
        if(g_uiProductModelCode<4)
        {
            LatticeInitialiseArrowDebreu(uiNextSlice);
            LatticeCalculateArrowDebreu(uiNextSlice);
        }
        else // for dual ccy hybrids and three-factor fx we use a trick to project the diffusion so we can escape
        {
            if(g_dADPLimit!=1.0)
            {
                LatticeInitialiseDiffusionProbability(uiNextSlice);
                LatticeCalculateDiffusionProbability(uiNextSlice);
            }
        }

        // Resetting from Y to X is not necessary for smoothing filtering since Y dimension is equidistant as well as X dimension
        // However resetting X-rates is needed for single currency fixed income 2 and 3 factors
        if(g_uiProductModelCode==2||g_uiProductModelCode==3) LatticeSetXRatesInSlice(uiSlice);
        if(g_uiProductModelCode==1) LatticeSetXRatesInSlice(uiSlice);
        LatticeFreeSlice(uiSlice);



    } // for loop over slices


    // Calibrate Last Slice
    unsigned int uiLastSlice=g_uiNumSlices-1;
    SetLatticeMarketDf(uiLastSlice+1);
    if(g_uiProductModelCode>3) SetLatticeMarketDfForeign(uiLastSlice+1);

    if(g_uiProductModelCode<4)
        LatticeCalibrateSlice(uiLastSlice);

    if(g_uiProductModelCode>3)
    {
        LatticeCalibrateSliceBase(uiLastSlice);
        LatticeCalibrateSliceForeign(uiLastSlice);
    }

    if(g_uiProductModelCode==2||g_uiProductModelCode==3) LatticeSetXRatesInSlice(uiLastSlice);

} // LatticeForwardConstruction(...)



void LabelLatticeWithDerivatives()
{
    unsigned int ui=0;
    double dShockCutOff=0.;

    if(g_uiProductModelCode<4)
    {

        switch(g_uiNumFactors)
        {
        case 1:
            {
                g_S1_TD.resize(g_uiNumSlices);
                g_S1_DERIVATIVE_TD.resize(g_uiNumSlices);
                for(ui=0;ui<g_uiNumSlices;ui++)
                {
                    g_S1_TD.at(ui)=numerical_function(g_LATTICE_DATE.at(ui),g_dStdDevBaseX,g_dStdDevBaseZ);
                }
                for(ui=0;ui<g_uiNumSlices;ui++)
                {
                    if(g_dNoticeDates.at(0)>dShockCutOff) g_S1_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_LATTICE_DATE,g_S1_TD);
                    else g_S1_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_dStdDevBaseX,g_dStdDevBaseZ);
                }
            }
            break;
        case 2:
            {
                g_S1_TD.resize(g_uiNumSlices);
                g_S1_DERIVATIVE_TD.resize(g_uiNumSlices);
                g_S2_TD.resize(g_uiNumSlices);
                g_S2_DERIVATIVE_TD.resize(g_uiNumSlices);
                for(ui=0;ui<g_uiNumSlices;ui++)
                {
                    g_S1_TD.at(ui)=numerical_function(g_LATTICE_DATE.at(ui),g_dStdDevBaseX,g_dStdDevBaseZ);
                    g_S2_TD.at(ui)=(g_S2/g_S1)*g_S1_TD.at(ui);
                }
                for(ui=0;ui<g_uiNumSlices;ui++)
                {
                    if(g_dNoticeDates.at(0)>dShockCutOff) g_S1_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_LATTICE_DATE,g_S1_TD);
                    else g_S1_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_dStdDevBaseX,g_dStdDevBaseZ);
                    g_S2_DERIVATIVE_TD.at(ui)=(g_S2/g_S1)*g_S1_DERIVATIVE_TD.at(ui);
                }
            }
            break;
        case 3:
            {
                g_S1_TD.resize(g_uiNumSlices);
                g_S1_DERIVATIVE_TD.resize(g_uiNumSlices);
                g_S2_TD.resize(g_uiNumSlices);
                g_S2_DERIVATIVE_TD.resize(g_uiNumSlices);
                g_S3_TD.resize(g_uiNumSlices);
                g_S3_DERIVATIVE_TD.resize(g_uiNumSlices);
                for(ui=0;ui<g_uiNumSlices;ui++)
                {
                    g_S1_TD.at(ui)=numerical_function(g_LATTICE_DATE.at(ui),g_dStdDevBaseX,g_dStdDevBaseZ);
                    g_S2_TD.at(ui)=(g_S2/g_S1)*g_S1_TD.at(ui);
                    g_S3_TD.at(ui)=(g_S3/g_S1)*g_S1_TD.at(ui);
                }
                for(ui=0;ui<g_uiNumSlices;ui++)
                {
                    if(g_dNoticeDates.at(0)>dShockCutOff) g_S1_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_LATTICE_DATE,g_S1_TD);
                    else g_S1_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_dStdDevBaseX,g_dStdDevBaseZ);
                    g_S2_DERIVATIVE_TD.at(ui)=(g_S2/g_S1)*g_S1_DERIVATIVE_TD.at(ui);
                    g_S3_DERIVATIVE_TD.at(ui)=(g_S3/g_S1)*g_S1_DERIVATIVE_TD.at(ui);
                }
            }
            break;
        }
    }
    else
    {
        switch(g_uiNumFactors)
        {
        case 1:
            {
                throw("You should never come here");
            }
            break;
        case 2:
            {
                g_S1_TD.resize(g_uiNumSlices);
                g_S1_DERIVATIVE_TD.resize(g_uiNumSlices);
                g_S2_TD.resize(g_uiNumSlices);
                g_S2_DERIVATIVE_TD.resize(g_uiNumSlices);
                for(ui=0;ui<g_uiNumSlices;ui++)
                {
                    g_S1_TD.at(ui)=numerical_function(g_LATTICE_DATE.at(ui),g_dStdDevBaseX,g_dStdDevBaseZ);
                    g_S2_TD.at(ui)=numerical_function(g_LATTICE_DATE.at(ui),g_dStdDevForeignX,g_dStdDevForeignZ);
                }
                for(ui=0;ui<g_uiNumSlices;ui++)
                {
                    if(g_dNoticeDates.at(0)>dShockCutOff) g_S1_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_LATTICE_DATE,g_S1_TD);
                    else g_S1_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_dStdDevBaseX,g_dStdDevBaseZ);
                    if(g_dNoticeDates.at(0)>dShockCutOff) g_S2_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_LATTICE_DATE,g_S2_TD);
                    else g_S2_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_dStdDevForeignX,g_dStdDevForeignZ);
                }
            }
            break;
        case 3:
            {
                g_S1_TD.resize(g_uiNumSlices);
                g_S1_DERIVATIVE_TD.resize(g_uiNumSlices);
                g_S2_TD.resize(g_uiNumSlices);
                g_S2_DERIVATIVE_TD.resize(g_uiNumSlices);
                g_S3_TD.resize(g_uiNumSlices);
                g_S3_DERIVATIVE_TD.resize(g_uiNumSlices);
                LatticeVolsGlobalExtraction.Resize(g_uiNumSlices,7);
                for(ui=0;ui<g_uiNumSlices;ui++)
                {
                    g_S1_TD.at(ui)=numerical_function(g_LATTICE_DATE.at(ui),g_dStdDevBaseX,g_dStdDevBaseZ);
                    g_S2_TD.at(ui)=numerical_function(g_LATTICE_DATE.at(ui),g_dStdDevForeignX,g_dStdDevForeignZ);
                    g_S3_TD.at(ui)=numerical_function(g_LATTICE_DATE.at(ui),g_dSpotFXVolDatesTD,g_dSpotFXVolTD);
                }
                for(ui=0;ui<g_uiNumSlices;ui++)
                {
                    if(g_dNoticeDates.at(0)>dShockCutOff) g_S1_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_LATTICE_DATE,g_S1_TD);
                    else g_S1_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_dStdDevBaseX,g_dStdDevBaseZ);
                    if(g_dNoticeDates.at(0)>dShockCutOff) g_S2_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_LATTICE_DATE,g_S2_TD);
                    else g_S2_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_dStdDevForeignX,g_dStdDevForeignZ);
                    if(g_dNoticeDates.at(0)>dShockCutOff) g_S3_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_LATTICE_DATE,g_S3_TD);
                    else g_S3_DERIVATIVE_TD.at(ui)=numerical_derivative(g_LATTICE_DATE.at(ui),g_dSpotFXVolDatesTD,g_dSpotFXVolTD);
                }
                for(ui=0;ui<g_uiNumSlices;ui++)
                {
                     LatticeVolsGlobalExtraction.Elt(ui,0)=365.0*g_LATTICE_DATE.at(ui);
                     LatticeVolsGlobalExtraction.Elt(ui,1)=g_S1_TD.at(ui);
                     LatticeVolsGlobalExtraction.Elt(ui,2)=g_S1_DERIVATIVE_TD.at(ui);
                     LatticeVolsGlobalExtraction.Elt(ui,3)=g_S2_TD.at(ui);
                     LatticeVolsGlobalExtraction.Elt(ui,4)=g_S2_DERIVATIVE_TD.at(ui);
                     LatticeVolsGlobalExtraction.Elt(ui,5)=g_S3_TD.at(ui);
                     LatticeVolsGlobalExtraction.Elt(ui,6)=g_S3_DERIVATIVE_TD.at(ui);
                }
           }
            break;
        }
    }
}

void LatticeInitialiseOrigin()
{
    g_LATTICE_ZCB_PRICE.resize(g_uiNumSlices+1);
    if(g_uiProductModelCode>3) g_LATTICE_ZCB_PRICE_FOREIGN.resize(g_uiNumSlices+1);

    g_LATTICE_RATE_STEP1.resize(g_uiNumSlices);
    if(g_uiNumFactors>1) g_LATTICE_RATE_STEP2.resize(g_uiNumSlices);
    if(g_uiNumFactors>2) g_LATTICE_RATE_STEP3.resize(g_uiNumSlices);

    g_LATTICE_MAX_STATE1.resize(g_uiNumSlices);
    if(g_uiNumFactors>1) g_LATTICE_MAX_STATE2.resize(g_uiNumSlices);
    if(g_uiNumFactors>2) g_LATTICE_MAX_STATE3.resize(g_uiNumSlices);

    g_LATTICE_MIN_STATE1.resize(g_uiNumSlices);
    if(g_uiNumFactors>1) g_LATTICE_MIN_STATE2.resize(g_uiNumSlices);
    if(g_uiNumFactors>2) g_LATTICE_MIN_STATE3.resize(g_uiNumSlices);

    g_LATTICE_NUMBER_OF_STATES1.resize(g_uiNumSlices);
    if(g_uiNumFactors>1) g_LATTICE_NUMBER_OF_STATES2.resize(g_uiNumSlices);
    if(g_uiNumFactors>2) g_LATTICE_NUMBER_OF_STATES3.resize(g_uiNumSlices);

    g_LATTICE_MAX_STATE1.at(0)=0;
    if(g_uiNumFactors>1) g_LATTICE_MAX_STATE2.at(0)=0;
    if(g_uiNumFactors>2) g_LATTICE_MAX_STATE3.at(0)=0;

    g_LATTICE_MIN_STATE1.at(0)=0;
    if(g_uiNumFactors>1) g_LATTICE_MIN_STATE2.at(0)=0;
    if(g_uiNumFactors>2) g_LATTICE_MIN_STATE3.at(0)=0;

    g_LATTICE_NUMBER_OF_STATES1.at(0)=1;
    if(g_uiNumFactors>1) g_LATTICE_NUMBER_OF_STATES2.at(0)=1;
    if(g_uiNumFactors>2) g_LATTICE_NUMBER_OF_STATES3.at(0)=1;

    g_LATTICE_X1_RATE.at(0).resize(1);
    if(g_uiNumFactors>1) g_LATTICE_X2_RATE.at(0).resize(1);
    if(g_uiNumFactors>2) g_LATTICE_X3_RATE.at(0).resize(1);

    g_LATTICE_G_RATE.at(0).resize(1);
    g_LATTICE_SHORT_RATE.at(0).resize(1);
    g_LATTICE_ARROW_DEBREU.at(0).resize(1);
    g_LATTICE_ARROW_DEBREU_PROB.at(0).resize(1);
    g_LATTICE_NODE_EXISTENCE.at(0).resize(1);

    if(g_uiProductModelCode>3)
    {
        g_LATTICE_G_RATE_FOREIGN.at(0).resize(1);
        g_LATTICE_SHORT_RATE_FOREIGN.at(0).resize(1);
    }

    if(g_uiProductModelCode==5)
    {
        g_LATTICE_G_RATE_FX.at(0).resize(1);
        g_LATTICE_SHORT_RATE_FX.at(0).resize(1);
    }

    g_LATTICE_TARGET_1.at(0).resize(1);
    if(g_uiNumFactors>1) g_LATTICE_TARGET_2.at(0).resize(1);
    if(g_uiNumFactors>2) g_LATTICE_TARGET_3.at(0).resize(1);

    g_LATTICE_ARROW_DEBREU.at(0).at(0)=1.;
    g_LATTICE_ARROW_DEBREU_PROB.at(0).at(0)=1.;

    g_LATTICE_UP_1.at(0).resize(1);
    if(g_uiNumFactors>1) g_LATTICE_UP_2.at(0).resize(1);
    if(g_uiNumFactors>2) g_LATTICE_UP_3.at(0).resize(1);

    g_LATTICE_MID_1.at(0).resize(1);
    if(g_uiNumFactors>1) g_LATTICE_MID_2.at(0).resize(1);
    if(g_uiNumFactors>2) g_LATTICE_MID_3.at(0).resize(1);

    g_LATTICE_DOWN_1.at(0).resize(1);
    if(g_uiNumFactors>1) g_LATTICE_DOWN_2.at(0).resize(1);
    if(g_uiNumFactors>2) g_LATTICE_DOWN_3.at(0).resize(1);

    if(g_uiNumFactors==1)
    {
        for(unsigned int uk=0;uk<3;uk++)
        {
            g_LATTICE_CONNECTION.at(0,uk).resize(1);
            g_LATTICE_GREEN.at(0,uk).resize(1);
        }
    }

    if(g_uiNumFactors==2)
    {
        for(unsigned int uk=0;uk<9;uk++)
        {
            g_LATTICE_CONNECTION.at(0,uk).resize(1);
            g_LATTICE_GREEN.at(0,uk).resize(1);
        }
    }

    if(g_uiNumFactors==3)
    {
        for(unsigned int uk=0;uk<27;uk++)
        {
            g_LATTICE_CONNECTION.at(0,uk).resize(1);
            g_LATTICE_GREEN.at(0,uk).resize(1);
        }
    }

} // LatticeInitialiseOrigin()

void LatticeInitialise()
{
    g_dStdDevCutOff.resize(g_uiNumSlices-1);
    g_dUpper.resize(g_uiNumSlices-1);
    g_dLower.resize(g_uiNumSlices-1);
    g_dNextSliceLimitUp.resize(g_uiNumSlices-1);
    g_dCurrentSliceLimitUp.resize(g_uiNumSlices-1);
    g_dNextSliceLimitDown.resize(g_uiNumSlices-1);
    g_dCurrentSliceLimitDown.resize(g_uiNumSlices-1);
    g_dCurrentCenter.resize(g_uiNumSlices-1);
    g_dNextCenter.resize(g_uiNumSlices-1);
    g_iLimit.resize(g_uiNumSlices-1);
    g_iCurrentShift.resize(g_uiNumSlices-1);
    g_iNextShift.resize(g_uiNumSlices-1);

    g_dDiffusion_11.resize(g_uiNumSlices-1);
    if(g_uiNumFactors>1)
    {
        g_dDiffusion_21.resize(g_uiNumSlices-1);
        g_dDiffusion_12.resize(g_uiNumSlices-1);
        g_dDiffusion_22.resize(g_uiNumSlices-1);
    }
    if(g_uiNumFactors>2)
    {
        g_dDiffusion_23.resize(g_uiNumSlices-1);
        g_dDiffusion_32.resize(g_uiNumSlices-1);
        g_dDiffusion_31.resize(g_uiNumSlices-1);
        g_dDiffusion_33.resize(g_uiNumSlices-1);
        g_dDiffusion_13.resize(g_uiNumSlices-1);
    }

    g_LATTICE_X1_RATE.resize(g_uiNumSlices);
    if(g_uiNumFactors>1) g_LATTICE_X2_RATE.resize(g_uiNumSlices);
    if(g_uiNumFactors>2) g_LATTICE_X3_RATE.resize(g_uiNumSlices);

    g_LATTICE_G_RATE.resize(g_uiNumSlices);
    g_LATTICE_SHORT_RATE.resize(g_uiNumSlices);
    g_LATTICE_ARROW_DEBREU.resize(g_uiNumSlices);

    if(g_uiProductModelCode>3)
    {
        g_LATTICE_G_RATE_FOREIGN.resize(g_uiNumSlices);
        g_LATTICE_SHORT_RATE_FOREIGN.resize(g_uiNumSlices);
    }

    if(g_uiProductModelCode==5)
    {
        g_LATTICE_G_RATE_FX.resize(g_uiNumSlices);
        g_LATTICE_SHORT_RATE_FX.resize(g_uiNumSlices);
    }


    g_LATTICE_ARROW_DEBREU_PROB.resize(g_uiNumSlices);
    g_LATTICE_NODE_EXISTENCE.resize(g_uiNumSlices);

    g_LATTICE_TARGET_1.resize(g_uiNumSlices);
    if(g_uiNumFactors>1) g_LATTICE_TARGET_2.resize(g_uiNumSlices);
    if(g_uiNumFactors>2) g_LATTICE_TARGET_3.resize(g_uiNumSlices);

    g_LATTICE_UP_1.resize(g_uiNumSlices);
    if(g_uiNumFactors>1) g_LATTICE_UP_2.resize(g_uiNumSlices);
    if(g_uiNumFactors>2) g_LATTICE_UP_3.resize(g_uiNumSlices);

    g_LATTICE_MID_1.resize(g_uiNumSlices);
    if(g_uiNumFactors>1) g_LATTICE_MID_2.resize(g_uiNumSlices);
    if(g_uiNumFactors>2) g_LATTICE_MID_3.resize(g_uiNumSlices);

    g_LATTICE_DOWN_1.resize(g_uiNumSlices);
    if(g_uiNumFactors>1) g_LATTICE_DOWN_2.resize(g_uiNumSlices);
    if(g_uiNumFactors>2) g_LATTICE_DOWN_3.resize(g_uiNumSlices);

    if(g_uiNumFactors==1)
    {
        g_LATTICE_CONNECTION.resize(g_uiNumSlices,3);
        g_LATTICE_GREEN.resize(g_uiNumSlices,3);
    }

    if(g_uiNumFactors==2)
    {
        g_LATTICE_CONNECTION.resize(g_uiNumSlices,9);
        g_LATTICE_GREEN.resize(g_uiNumSlices,9);
    }

    if(g_uiNumFactors==3)
    {
        g_LATTICE_CONNECTION.resize(g_uiNumSlices,27);
        g_LATTICE_GREEN.resize(g_uiNumSlices,27);
    }

} // LatticeInitialise()

void U3U3U1_Projection_Base()
{

    // Do first the base currency
    // Move to one-factor base mode
    g_uiNumFactors=1;
    g_uiProductModelCode=1;
    // Create a virtual base currency 1D lattice
    // g_dShiftShortGRateBase.resize(g_uiNumSlices);
    LatticeInitialise();
    LatticeInitialiseOrigin();
    LabelLatticeWithDates();
    LabelLatticeWithDerivatives();
    InitLatticeModel();
    LabelLatticeLangevinDiffusion();
    SetLatticeMarketDf(0);
    SetLatticeRateStep(0, 0.0);
    DKMaille<long> indexStates(g_uiNumFactors,0L);
    PrepareTheLatticeNode(0, indexStates);
    LatticeForwardConstructionOfSlices();
    // Now get the shifts for each slice
    DKMaille<long> iMiddle(1); iMiddle.at(0)=0;
    // This works only for Gaussian rates modify appropriately for DK model
    g_d1DGreenBase.resize(g_uiNumSlices);
    g_d1DShortRateBase.resize(g_uiNumSlices);
    g_dShiftShortGRateBase.resize(g_uiNumSlices);
    for(unsigned ui=0;ui<g_uiNumSlices;ui++)
    {
        g_dShiftShortGRateBase.at(ui)=g_dShiftShortGRate.at(ui);
        int iFlag=iIsNoticeDate(ui,g_dNoticeDates);
        // Fill one-dimensional Green's function
        // if(iFlag>-1) // if notice date
        {
            unsigned int uiSize=g_LATTICE_NUMBER_OF_STATES1.at(ui);
            g_d1DGreenBase.at(ui).resize(uiSize);
            g_d1DShortRateBase.at(ui).resize(uiSize);
            for(unsigned int um=0;um<uiSize;um++)
            {
                g_d1DGreenBase.at(ui).at(um)=g_LATTICE_ARROW_DEBREU.at(ui).at(um);
                g_d1DShortRateBase.at(ui).at(um)=g_LATTICE_SHORT_RATE.at(ui).at(um);
            }
        }
    }
} // U3U3U1_Projection_Base()

void U3U3U1_Projection_Foreign(unsigned int uiForeignCurveSize)
{
    g_dQParameter_Base=g_dQParameter_Foreign;
    // Move to one-factor base mode
    g_uiNumFactors=1;
    g_uiProductModelCode=1;
    // Create a virtual foreign currency 1D lattice
    // g_dShiftShortGRateForeign.resize(g_uiNumSlices);
    // ///////////////////////////////////////////////////////////
    // Switch Temporarily Foreign for Base Including Vol Parameters
    // resize curves
    g_dDiscountCurveBase.resize(uiForeignCurveSize);
    g_dDiscountCurveBaseDates.resize(uiForeignCurveSize);
    g_dZCCurveBase.resize(uiForeignCurveSize);
    // assign curves
    g_dDiscountCurveBaseDates=g_dDiscountCurveForeignDates;
    g_dDiscountCurveBase=g_dDiscountCurveForeign;
    g_dZCCurveBase=g_dZCCurveForeign;
    /////////////////////////////////////////////////////////////////
    g_dEquilibriumForwardBase=g_dEquilibriumForwardForeign;
    g_dEquilibriumForwardShortBase=g_dEquilibriumForwardShortForeign;
    ////////////////////////////////////////////////////////////////
    // switch parameters
    g_S1=g_S2;
    g_MR1=g_MR2;
    g_dLimitXX1=g_dLimitXX2;
    g_dStdDevBaseX.resize(g_dStdDevForeignX.entries());
    g_dStdDevBaseZ.resize(g_dStdDevForeignZ.entries());
    g_dStdDevBaseX=g_dStdDevForeignX;
    g_dStdDevBaseZ=g_dStdDevForeignZ;
    //////////////////////////////////////////////////////////////
    LatticeInitialise();
    LatticeInitialiseOrigin();
    LabelLatticeWithDates();
    LabelLatticeWithDerivatives();
    InitLatticeModel();
    LabelLatticeLangevinDiffusion();
    SetLatticeMarketDf(0);
    SetLatticeRateStep(0, 0.0);
    DKMaille<long> indexStates(g_uiNumFactors,0L);
    PrepareTheLatticeNode(0, indexStates);
    LatticeForwardConstructionOfSlices();
    // Now get the shifts for each slice
    DKMaille<long> iMiddle(1); iMiddle.at(0)=0;
    g_dFromForeignDriftToCash.resize(g_uiNumSlices);
    // Calculate now the domestic currency adjustment so that drift is expressed in cash terms
    // The below is the result of an integration
    // Modify this for
    // 1. Time dependent case 2. DK model case
    double dFutureDate=0.;
    // Get the time-dependent FX Vol
    DKMaille<double> S3_local(g_uiNumSlices);
    // REMEMBER TO MODIFY THIS PART FOR NON-STANDARD DIFFUSIONS EVEN THOUGH IT IS NOT NECESSARY
    for(unsigned ui=0;ui<g_uiNumSlices;ui++)
    {
        S3_local.at(ui)=numerical_function(g_LATTICE_DATE.at(ui),g_dSpotFXVolDatesTD,g_dSpotFXVolTD);
    }
    for(ui=0;ui<g_uiNumSlices;ui++)
    {
        dFutureDate=g_LATTICE_DATE.at(ui);
        double dX1=-g_RHO2*g_S2*g_S3*(1.-exp(-g_MR2*dFutureDate))/g_MR2;
        double dX2=-g_RHO2*SpotFXVolIRStdDevQuantoIntegral(g_LATTICE_DATE,
                   S3_local, // Spot FX vol
                   g_S1_TD, // Foreign IR vol
                   g_MR2, // Foreign MR
                   0., // Start of integration
                   g_LATTICE_DATE.at(ui)); // end of integration
        g_dFromForeignDriftToCash.at(ui)=dX2;
    }
    // This works only for Gaussian rates modify appropriately for DK model
    g_d1DGreenForeign.resize(g_uiNumSlices);
    g_d1DShortRateForeign.resize(g_uiNumSlices);
    g_dShiftShortGRateForeign.resize(g_uiNumSlices);
    for(ui=0;ui<g_uiNumSlices;ui++)
    {
        g_dShiftShortGRateForeign.at(ui)=g_dShiftShortGRate.at(ui)+g_dFromForeignDriftToCash.at(ui);
        // g_dShiftShortGRateForeign.at(ui)=g_LATTICE_SHORT_RATE.at(ui).at(LatticeIndexConvert(ui,iMiddle))+g_dFromForeignDriftToCash.at(ui);
        int iFlag=iIsNoticeDate(ui,g_dNoticeDates);
        // Fill one-dimensional Green's function
        // if(iFlag>-1) // if notice date
        {
            unsigned int uiSize=g_LATTICE_NUMBER_OF_STATES1.at(ui);
            g_d1DGreenForeign.at(ui).resize(uiSize);
            g_d1DShortRateForeign.at(ui).resize(uiSize);
            for(unsigned int um=0;um<uiSize;um++)
            {
                g_d1DGreenForeign.at(ui).at(um)=g_LATTICE_ARROW_DEBREU.at(ui).at(um);
                g_d1DShortRateForeign.at(ui).at(um)=g_LATTICE_SHORT_RATE.at(ui).at(um);
            }
        }
    }
} // U3U3U1_Projection_Foreign()


void U3U3U3_Reinstate(double dSpotDate,
                      double dStdDevBase,
                      double dMeanReversionBase,
                      double dX1Limit,
                      DKMaille<double> &dBaseDates,
                      DKMaille<double> &dBaseRates,
                      DKMaille<double> &dStdDevBaseX,
                      DKMaille<double> &dStdDevBaseZ,
                      double dProductModelCode,
                      DKMaille<double> &dShortRateKOnSlicesBase,
                      DKMaille<double> &dShortRateOnSlicesBase,
                      double dSkewForeign)
{
    g_dQParameter_Base=dSkewForeign;
    g_S1=dStdDevBaseZ.at(0);
    g_MR1=dMeanReversionBase;
    g_dStdDevBaseX.resize(dStdDevBaseX.entries());
    g_dStdDevBaseZ.resize(dStdDevBaseZ.entries());
    g_dStdDevBaseX=dStdDevBaseX;
    g_dStdDevBaseZ=dStdDevBaseZ;
    g_dLimitXX1=dX1Limit;
    // Put the curves back
    g_dDiscountCurveBase.resize(dBaseDates.entries());
    g_dDiscountCurveBaseDates.resize(dBaseDates.entries());
    g_dZCCurveBase.resize(dBaseDates.entries());
    g_dEquilibriumForwardBase=dShortRateKOnSlicesBase;
    g_dEquilibriumForwardShortBase=dShortRateOnSlicesBase;
    for(unsigned int ui=0;ui<dBaseDates.entries();ui++)
    {
        g_dDiscountCurveBase.at(ui)=dBaseRates.at(ui);
        g_dDiscountCurveBaseDates.at(ui)=(dBaseDates.at(ui)-dSpotDate)/365.;
        if(g_dDiscountCurveBaseDates.at(ui)>0.) g_dZCCurveBase.at(ui)=-log(g_dDiscountCurveBase.at(ui))/g_dDiscountCurveBaseDates.at(ui);
        else g_dZCCurveBase.at(ui)=-log(g_dDiscountCurveBase.at(ui+1))/g_dDiscountCurveBaseDates.at(ui+1);
    }
    g_GridScaling1=2./3.;
    g_GridScaling2=2./3.;
    g_GridScaling3=2./3.;
    // Switch back to the desired configuration
    g_uiNumFactors=(unsigned int)dProductModelCode;
    g_uiProductModelCode=(unsigned int)dProductModelCode;
    // Dual currency fixed income only
    if(g_uiNumFactors==4) g_uiNumFactors=2;
    // Dual currency fixed income and foreign exchange
    if(g_uiNumFactors==5) g_uiNumFactors=3;
} // U3U3U3_Reinstate()



double LatticeFXOptionAnalytics(unsigned int uii,
                                unsigned int uk,
                                DKMaille<double> &dDecayDomestic,
                                DKMaille<double> &dDecayForeign,
                                double dOptionExpiry,
                                double dOptionPayment,
                                double dStrike,
                                double dCap,
                                DKMaille<double> &dNormBase,
                                DKMaille<double> &dNormForeign,
                                double dCallPut)
{
    // Introduction here of Green's Function Analytics for fast discounting of bonds and fx options
    // Expiry date of the fx option at the end of the lattice as a test
    double dDate=g_LATTICE_DATE.at(uii);
    double dExpiry=dOptionExpiry;
    // This is also the payment date of the coupon (to be changed)
    // Calculate today's zero-coupon domestic and foreign bonds
    double dBaseDiscount;
    double dForeignDiscount;
    double dForwardFX;
    double dVolFX;
    double dOptionFX;
    double dDiscount;
    dBaseDiscount=BondFromVFDKAnalytics_HW1F(dNormBase.at(0),dDecayDomestic.at(0),g_LATTICE_SHORT_RATE.at(uii).at(uk));
    dForeignDiscount=BondFromVFDKAnalytics_HW1F(dNormForeign.at(0),dDecayForeign.at(0),g_LATTICE_SHORT_RATE_FOREIGN.at(uii).at(uk));
    dForwardFX=g_LATTICE_SHORT_RATE_FX.at(uii).at(uk)*dForeignDiscount/dBaseDiscount;
    // Constant Spot FX Vol
    // 03.05.2003 Upgrade to time-dependent Spot FX Vol

    if(g_uiVolIsCalculated.at(uii)==0 || uii==0) /// JMP to allow analytical computation
    {
        if(g_uiIsDK==0)
        {
            g_dLocalFXVol=GetForwardFXVol(0.,dDate,dExpiry,dExpiry,g_dStrip,
                                          g_dStripDomesticStdDev,g_dStripForeignStdDev,
                                          g_dStripSpotFXVol,g_MR1,g_MR2,g_RHO1,g_RHO2,g_RHO);
        }
        else
        {
            g_dLocalFXVol=GetForwardFXVol_DK(0.,dDate,dExpiry,dExpiry,g_dStrip,
                                             g_dStripDomesticStdDev,g_dStripForeignStdDev,
                                             g_dStripSpotFXVol,g_MR1,g_MR2,g_RHO1,g_RHO2,g_RHO,
                                             g_LATTICE_DATE,g_dEquilibriumForwardShortBase,g_dEquilibriumForwardShortForeign);
        }

        if(uii>0)
            g_uiVolIsCalculated.at(uii)=1;
    }
    dVolFX=g_dLocalFXVol;
    if(dOptionExpiry!=dOptionPayment)
    {
        dDiscount=BondFromVFDKAnalytics_HW1F(dNormBase.at(1),dDecayDomestic.at(1),g_LATTICE_SHORT_RATE.at(uii).at(uk));
    }
    else
    {
        dDiscount=dBaseDiscount;
    }
    dOptionFX=dBlackScholes((dExpiry-dDate),dVolFX,dStrike,dForwardFX,0.)*dDiscount;
    double dOptionCap=dBlackScholes((dExpiry-dDate),dVolFX,dCap,dForwardFX,0.)*dDiscount;
    dOptionFX=dOptionFX-dOptionCap;
    if(dCallPut==1.0) dOptionFX=dOptionFX-dDiscount*(dForwardFX-dStrike);
    return dOptionFX;
}

double LatticeCashFlowAnalytics(unsigned int uii,
                                unsigned int uk,
                                DKMaille<double> &dDecayDomestic,
                                DKMaille<double> &dDecayForeign,
                                double dPaymentDate,
                                DKMaille<double> &dNormBase,
                                DKMaille<double> &dNormForeign)
{
    // Introduction here of Green's Function Analytics for fast discounting of bonds
    double dBaseDiscount;
    dBaseDiscount=BondFromVFDKAnalytics_HW1F(dNormBase.at(0),dDecayDomestic.at(0),g_LATTICE_SHORT_RATE.at(uii).at(uk));
    return dBaseDiscount;
}

double LatticeFXForwardAnalytics(unsigned int uii,
                                 unsigned int uk,
                                 DKMaille<double> &dDecayDomestic,
                                 DKMaille<double> &dDecayForeign,
                                 double dReset,
                                 double dPayment,
                                 double dStrike,
                                 DKMaille<double> &dNormBase,
                                 DKMaille<double> &dNormForeign,
                                 double dLongShort)
{
    // Introduction here of Green's Function Analytics for fast discounting of bonds and fx options
    // Expiry date of the fx option at the end of the lattice as a test
    double dDate=g_LATTICE_DATE.at(uii);
    // This is also the payment date of the coupon (to be changed)
    // Calculate today's zero-coupon domestic and foreign bonds
    double dBaseDiscount;
    double dForeignDiscount;
    double dForwardFX;
    double dFX;
    double dDiscount;
    dBaseDiscount=BondFromVFDKAnalytics_HW1F(dNormBase.at(0),dDecayDomestic.at(0),g_LATTICE_SHORT_RATE.at(uii).at(uk));
    dForeignDiscount=BondFromVFDKAnalytics_HW1F(dNormForeign.at(0),dDecayForeign.at(0),g_LATTICE_SHORT_RATE_FOREIGN.at(uii).at(uk));
    dForwardFX=g_LATTICE_SHORT_RATE_FX.at(uii).at(uk)*dForeignDiscount/dBaseDiscount;
    if(dReset!=dPayment)
    {
        dDiscount=BondFromVFDKAnalytics_HW1F(dNormBase.at(1),dDecayDomestic.at(1),g_LATTICE_SHORT_RATE.at(uii).at(uk));
    }
    else
    {
        dDiscount=dBaseDiscount;
    }
    dFX=(dForwardFX-dStrike)*dDiscount;
    if(dLongShort==1.0) dFX=(-dForwardFX+dStrike)*dDiscount;
    return dFX;
}

double LatticeFloatSwapIndexAnalytics(unsigned int uii,
                                      unsigned int uk,
                                      DKMaille2D<double> &dDecay,
                                      double dStartDate,
                                      double dEndDate,
                                      double dPaymentDate,
                                      double dAccrualBasis,
                                      double dFundingSpread,
                                      DKMaille<double> &dNorm)
{
    double dDiscount1=0.;
    double dDiscount2=0.;
    double dDiscount3=0.;
    double dDate=g_LATTICE_DATE.at(uii);
    if(g_uiProductModelCode>3)
    {
        if(dDate!=dStartDate)
        {
            dDiscount1=BondFromVFDKAnalytics_HW1F(dNorm.at(0),dDecay.at(0,0),g_LATTICE_SHORT_RATE.at(uii).at(uk));
        }
        else
            dDiscount1=1.;
        if(g_dIsFundingLegStochastic==0.) dDiscount1=BaseDiscountInterpolate(dStartDate)/BaseDiscountInterpolate(dDate);

        dDiscount2=BondFromVFDKAnalytics_HW1F(dNorm.at(1),dDecay.at(1,0),g_LATTICE_SHORT_RATE.at(uii).at(uk));
        if(g_dIsFundingLegStochastic==0.) dDiscount2=BaseDiscountInterpolate(dEndDate)/BaseDiscountInterpolate(dDate);

        if(dPaymentDate!=dDate)
        {
            if(dPaymentDate!=dEndDate)
            {
                dDiscount3=BondFromVFDKAnalytics_HW1F(dNorm.at(2),dDecay.at(2,0),g_LATTICE_SHORT_RATE.at(uii).at(uk));
            }
            else
                dDiscount3=dDiscount2;
        }
        else
            dDiscount3=1.;
        if(g_dIsFundingLegStochastic==0.) dDiscount3=BaseDiscountInterpolate(dPaymentDate)/BaseDiscountInterpolate(dDate);
    }
    if(g_uiProductModelCode<=3)
    {
        if(dDate!=dStartDate)
        {
            if(g_uiNumFactors==1) dDiscount1=BondFromVFDKAnalytics_HW1F(dNorm.at(0),
                                                 dDecay.at(0,0),
                                                 g_LATTICE_SHORT_RATE.at(uii).at(uk));
            if(g_uiNumFactors==2) dDiscount1=BondFromVFDKAnalytics_HW2F(dNorm.at(0),
                                                 dDecay.at(0,0),
                                                 g_LATTICE_X1_RATE.at(uii).at(uk),
                                                 dDecay.at(0,1),
                                                 g_LATTICE_X2_RATE.at(uii).at(uk));
            if(g_uiNumFactors==3) dDiscount1=BondFromVFDKAnalytics_HW3F(dNorm.at(0),
                                                 dDecay.at(0,0),
                                                 g_LATTICE_X1_RATE.at(uii).at(uk),
                                                 dDecay.at(0,1),
                                                 g_LATTICE_X2_RATE.at(uii).at(uk),
                                                 dDecay.at(0,2),
                                                 g_LATTICE_X3_RATE.at(uii).at(uk));
        }
        else
            dDiscount1=1.;
        if(g_uiNumFactors==1) dDiscount2=BondFromVFDKAnalytics_HW1F(dNorm.at(1),
                                             dDecay.at(1,0),
                                             g_LATTICE_SHORT_RATE.at(uii).at(uk));
        if(g_uiNumFactors==2) dDiscount2=BondFromVFDKAnalytics_HW2F(dNorm.at(1),
                                             dDecay.at(1,0),
                                             g_LATTICE_X1_RATE.at(uii).at(uk),
                                             dDecay.at(1,1),
                                             g_LATTICE_X2_RATE.at(uii).at(uk));
        if(g_uiNumFactors==3) dDiscount2=BondFromVFDKAnalytics_HW3F(dNorm.at(1),
                                             dDecay.at(1,0),
                                             g_LATTICE_X1_RATE.at(uii).at(uk),
                                             dDecay.at(1,1),
                                             g_LATTICE_X2_RATE.at(uii).at(uk),
                                             dDecay.at(1,2),
                                             g_LATTICE_X3_RATE.at(uii).at(uk));
        if(dPaymentDate!=dDate)
        {
            if(dPaymentDate!=dEndDate)
            {
                if(g_uiNumFactors==1) dDiscount3=BondFromVFDKAnalytics_HW1F(dNorm.at(2),
                                                     dDecay.at(2,0),
                                                     g_LATTICE_SHORT_RATE.at(uii).at(uk));
                if(g_uiNumFactors==2) dDiscount3=BondFromVFDKAnalytics_HW2F(dNorm.at(2),
                                                     dDecay.at(2,0),
                                                     g_LATTICE_X1_RATE.at(uii).at(uk),
                                                     dDecay.at(2,1),
                                                     g_LATTICE_X2_RATE.at(uii).at(uk));
                if(g_uiNumFactors==3) dDiscount3=BondFromVFDKAnalytics_HW3F(dNorm.at(2),
                                                     dDecay.at(2,0),
                                                     g_LATTICE_X1_RATE.at(uii).at(uk),
                                                     dDecay.at(2,1),
                                                     g_LATTICE_X2_RATE.at(uii).at(uk),
                                                     dDecay.at(2,2),
                                                     g_LATTICE_X3_RATE.at(uii).at(uk));
            }
            else
                dDiscount3=dDiscount2;
        }
        else
            dDiscount3=1.;
    }
    return ((dDiscount1/dDiscount2)-1.)*dDiscount3+dFundingSpread*dAccrualBasis*dDiscount3;
}

double GetSwaptionVol(unsigned int uii,double dStartDate,double dEndDate,double dPaymentDate)
{
    double dSwaptionTenorInYears=dEndDate-dStartDate; //
    double dSwaptionExpiryInYears=dStartDate-g_LATTICE_DATE.at(uii); //
    double dNoticePeriodInDays=0.; // arbitrarily set to 0 in general 1-2 days lag depending on currency
    double dJulianObservationDate=g_dSpotDate+g_LATTICE_DATE.at(uii)*365.;
    double dAbsVol=0.;
    if(g_dIsSwaptionCalibrationWithBasis==1.0)
    {
        dAbsVol=SwaptionAbsVol_VFDK_HW1To3F(dSwaptionExpiryInYears,
                                            dSwaptionTenorInYears,
                                            dNoticePeriodInDays,
                                            g_dDiscountCurveBaseDates_FixedIncomeExotics,
                                            g_dDiscountCurveBaseNoBasis_FixedIncomeExotics,
                                            g_dDiscountCurveBase_FixedIncomeExotics,
                                            g_dStripJulian,
                                            g_dStripDomesticStdDev,
                                            g_dModelParameters_FixedIncomeExotics,
                                            dJulianObservationDate);
    }
    else
    {
        dAbsVol=SwaptionAbsVol_VFDK_HW1To3F(dSwaptionExpiryInYears,
                                            dSwaptionTenorInYears,
                                            dNoticePeriodInDays,
                                            g_dDiscountCurveBaseDates_FixedIncomeExotics,
                                            g_dDiscountCurveBaseNoBasis_FixedIncomeExotics,
                                            g_dDiscountCurveBaseNoBasis_FixedIncomeExotics,
                                            g_dStripJulian,
                                            g_dStripDomesticStdDev,
                                            g_dModelParameters_FixedIncomeExotics,
                                            dJulianObservationDate);
    }
    return dAbsVol;
}


double LatticeInverseFRNAnalytics(unsigned int uii,
                                  unsigned int uk,
                                  DKMaille2D<double> &dDecay,
                                  double dStartDate,
                                  double dEndDate,
                                  double dPaymentDate,
                                  double dAccrualBasis,
                                  double dFundingSpread,
                                  double dOptionSpread,
                                  double dStrike,
                                  double dFactor,
                                  DKMaille<double> &dNorm)
{
    double dDiscount1=0.;
    double dDiscount2=0.;
    double dDiscount3=0.;
    double dDate=g_LATTICE_DATE.at(uii);

    if(g_uiProductModelCode>3)
        throw("Callable Inverse Floater Priced with Single Currency Models only");

    if(g_uiProductModelCode<=3)
    {
        if(dDate!=dStartDate)
        {
            if(g_uiNumFactors==1) dDiscount1=BondFromVFDKAnalytics_HW1F(dNorm.at(0),
                                                 dDecay.at(0,0),
                                                 g_LATTICE_SHORT_RATE.at(uii).at(uk));
            if(g_uiNumFactors==2) dDiscount1=BondFromVFDKAnalytics_HW2F(dNorm.at(0),
                                                 dDecay.at(0,0),
                                                 g_LATTICE_X1_RATE.at(uii).at(uk),
                                                 dDecay.at(0,1),
                                                 g_LATTICE_X2_RATE.at(uii).at(uk));
            if(g_uiNumFactors==3) dDiscount1=BondFromVFDKAnalytics_HW3F(dNorm.at(0),
                                                 dDecay.at(0,0),
                                                 g_LATTICE_X1_RATE.at(uii).at(uk),
                                                 dDecay.at(0,1),
                                                 g_LATTICE_X2_RATE.at(uii).at(uk),
                                                 dDecay.at(0,2),
                                                 g_LATTICE_X3_RATE.at(uii).at(uk));
        }
        else
            dDiscount1=1.;
        if(g_uiNumFactors==1) dDiscount2=BondFromVFDKAnalytics_HW1F(dNorm.at(1),
                                             dDecay.at(1,0),
                                             g_LATTICE_SHORT_RATE.at(uii).at(uk));
        if(g_uiNumFactors==2) dDiscount2=BondFromVFDKAnalytics_HW2F(dNorm.at(1),
                                             dDecay.at(1,0),
                                             g_LATTICE_X1_RATE.at(uii).at(uk),
                                             dDecay.at(1,1),
                                             g_LATTICE_X2_RATE.at(uii).at(uk));
        if(g_uiNumFactors==3) dDiscount2=BondFromVFDKAnalytics_HW3F(dNorm.at(1),
                                             dDecay.at(1,0),
                                             g_LATTICE_X1_RATE.at(uii).at(uk),
                                             dDecay.at(1,1),
                                             g_LATTICE_X2_RATE.at(uii).at(uk),
                                             dDecay.at(1,2),
                                             g_LATTICE_X3_RATE.at(uii).at(uk));
        if(dPaymentDate!=dDate)
        {
            if(dPaymentDate!=dEndDate)
            {
                if(g_uiNumFactors==1) dDiscount3=BondFromVFDKAnalytics_HW1F(dNorm.at(2),
                                                     dDecay.at(2,0),
                                                     g_LATTICE_SHORT_RATE.at(uii).at(uk));
                if(g_uiNumFactors==2) dDiscount3=BondFromVFDKAnalytics_HW2F(dNorm.at(2),
                                                     dDecay.at(2,0),
                                                     g_LATTICE_X1_RATE.at(uii).at(uk),
                                                     dDecay.at(2,1),
                                                     g_LATTICE_X2_RATE.at(uii).at(uk));
                if(g_uiNumFactors==3) dDiscount3=BondFromVFDKAnalytics_HW3F(dNorm.at(2),
                                                     dDecay.at(2,0),
                                                     g_LATTICE_X1_RATE.at(uii).at(uk),
                                                     dDecay.at(2,1),
                                                     g_LATTICE_X2_RATE.at(uii).at(uk),
                                                     dDecay.at(2,2),
                                                     g_LATTICE_X3_RATE.at(uii).at(uk));
            }
            else
                dDiscount3=dDiscount2;
        }
        else
            dDiscount3=1.;
    }
    double dRate=((dDiscount1/dDiscount2)-1.)/dAccrualBasis+dFundingSpread;
    double dCoupon=0;
    // Determine whether we calculate or not the option
    // This is based on whether there is a notice period greater than 2 days
    // If yes we have to calculate the option price using analytics
    // we also have to calculate any convexity adjustment
    if(g_LATTICE_DATE.at(uii)-dStartDate<(-3./365.))
    {
        // Calculate Absolute Vol in forward observation mode
        double dSwaptionVol=0.;
        if(g_uiVolIsCalculated.at(uii)==0)
        {
            g_dLocalSwaptionVol=GetSwaptionVol(uii,dStartDate,dEndDate,dPaymentDate);
            g_uiVolIsCalculated.at(uii)=1;
        }
        // Calculate Convexity Adjustment
        double dConvexityAdjustment=0.;
        dRate+=dConvexityAdjustment;
        dRate*=dFactor;
        double dAnnuity=dAccrualBasis*dDiscount3;
        dCoupon=GetOptionPrice(1.,dRate,dAnnuity,dStrike,g_dLocalSwaptionVol,-g_LATTICE_DATE.at(uii)+dStartDate)/10000.+dOptionSpread*dAccrualBasis*dDiscount3;
    }
    else // Normal Calculation in the lattice
    {
        if(dRate<=dStrike)
        {
            dCoupon=(dStrike-dRate*dFactor)*dAccrualBasis*dDiscount3+dOptionSpread*dAccrualBasis*dDiscount3;
        }
        else
            dCoupon=dOptionSpread*dAccrualBasis*dDiscount3;
    }
    return dCoupon;
}

double LatticeFixedSwapIndexAnalytics(unsigned int uii,
                                      unsigned int uk,
                                      double dNorm,
                                      DKMaille<double> &dDecay,
                                      double dCouponPaymentDate,
                                      double dAccrualBasis,
                                      double dCashFlow)
{
    if(g_uiProductModelCode>3)
        throw("You cannot use LatticeFixedSwapIndexAnalytics for hybrid products");
    double dDiscount3=0.;
    switch(g_uiNumFactors)
    {
    case 1:
        {
            if(dCouponPaymentDate!=g_LATTICE_DATE.at(uii))
                dDiscount3=BondFromVFDKAnalytics_HW1F(dNorm,
                                                      dDecay.at(0),
                                                      g_LATTICE_SHORT_RATE.at(uii).at(uk));
            else
                dDiscount3=1.;
        }
        break;
    case 2:
        {
            if(dCouponPaymentDate!=g_LATTICE_DATE.at(uii))
                dDiscount3=BondFromVFDKAnalytics_HW2F(dNorm,
                                                      dDecay.at(0),
                                                      g_LATTICE_X1_RATE.at(uii).at(uk),
                                                      dDecay.at(1),
                                                      g_LATTICE_X2_RATE.at(uii).at(uk));
            else
                dDiscount3=1.;
        }
        break;
    case 3:
        {
            if(dCouponPaymentDate!=g_LATTICE_DATE.at(uii))
                dDiscount3=BondFromVFDKAnalytics_HW3F(dNorm,
                                                      dDecay.at(0),
                                                      g_LATTICE_X1_RATE.at(uii).at(uk),
                                                      dDecay.at(1),
                                                      g_LATTICE_X2_RATE.at(uii).at(uk),
                                                      dDecay.at(2),
                                                      g_LATTICE_X3_RATE.at(uii).at(uk));
            else
                dDiscount3=1.;
        }
        break;
    }
    return dCashFlow*dDiscount3*dAccrualBasis;
}
double LatticeFixedSwapIndexAnalyticsDomestic(unsigned int uii,
        unsigned int uk,
        double dNorm,
        double dDecay,
        double dCouponPaymentDate,
        double dAccrualBasis,
        double dCashFlow)
{
    double dDiscount3=0.;
    if(dCouponPaymentDate!=g_LATTICE_DATE.at(uii))
        dDiscount3=BondFromVFDKAnalytics_HW1F(dNorm,dDecay,g_LATTICE_SHORT_RATE.at(uii).at(uk));
    else
        dDiscount3=1.;

    return dCashFlow*dDiscount3*dAccrualBasis;
}
double LatticeFixedSwapIndexAnalyticsForeign(unsigned int uii,
        unsigned int uk,
        double dNorm,
        double dDecay,
        double dNormForeign,
        double dDecayForeign,
        double dCouponPaymentDate,
        double dAccrualBasis,
        double dCashFlow)
{
    double dDiscount3Base=0.;
    double dDiscount3Foreign=0.;
    //  Calculate domestic discount
    if(dCouponPaymentDate!=g_LATTICE_DATE.at(uii))
        dDiscount3Base=BondFromVFDKAnalytics_HW1F(dNorm,dDecay,g_LATTICE_SHORT_RATE.at(uii).at(uk));
    else
        dDiscount3Base=1.;
    //  Calculate foreign discount
    if(dCouponPaymentDate!=g_LATTICE_DATE.at(uii))
        dDiscount3Foreign=BondFromVFDKAnalytics_HW1F(dNormForeign,dDecayForeign,g_LATTICE_SHORT_RATE_FOREIGN.at(uii).at(uk));
    else
        dDiscount3Foreign=1.;
    // Calculate forward FX
    double dForwardFX=g_LATTICE_SHORT_RATE_FX.at(uii).at(uk)*dDiscount3Foreign/dDiscount3Base;
    return dCashFlow*dDiscount3Base*dAccrualBasis*dForwardFX;
}
int uiFirstFXCouponAfterNoticeDate(double dNoticeDate,DKMaille<double> &dFXResetDates,DKMaille<double> &dFXPaymentDates)
{
    int iCash=-1;
    unsigned int ui=0;
    while(ui<dFXResetDates.entries()&&(iCash==-1))
    {
        if((dNoticeDate<=dFXResetDates.at(ui)&&dFXPaymentDates.at(ui)>dFXResetDates.at(ui)+2./12.) // in advance
                ||(dNoticeDate+2./12.<dFXResetDates.at(ui))) // in arrears
        {
            iCash=(int)ui;
        }
        ui++;
    }
    return iCash;
}
int uiFirstIndexCouponAfterNoticeDate(double dNoticeDate,DKMaille<double> &dIndexResetDates)
{
    int iCash=-1;
    unsigned int ui=0;
    while(ui<dIndexResetDates.entries()&&(iCash==-1))
    {
        // Reset is 2 days before start as all funding is converted in JPY
        if(dNoticeDate<=dIndexResetDates.at(ui))
        {
            iCash=(int)ui;
        }
        ui++;
    }
    return iCash;
}
/*
// old functionality
int uiFirstCouponAfterNoticeDate(double dNoticeDate,DKMaille<double> &dIndexPaymentDates)
{
  int iCash=-1;
	unsigned int ui=0;
	while(ui<dIndexPaymentDates.entries()&&(iCash==-1))
	{
			// write this properly in the future
			if(dNoticeDate+30./365.<dIndexPaymentDates.at(ui)) 
			{
				iCash=(int)ui;
			}			
			ui++;
	}
	return iCash;
}
*/

/*
int uiFirstFXCouponAfterNoticeDate(double dNoticeDate, 
																   DKMaille<double> &dFXResetDates, 
																	 DKMaille<double> &dFXPaymentDates)
{
	int i, n = dFXResetDates.entries(), ResetInAdvance;
	double CouponTenor;

	// Evaluate once and for all whether the set is in advance:
	if (n > 1) // At least two FX coupons in the set.
	{
		// All periods being contiguous, if period 1's set is
		// before period 0's payment, then the set is in advance.
		ResetInAdvance = (dFXResetDates.at(1) < dFXPaymentDates.at(0));
		// Round the coupon tenor to a multiple of 12:
		CouponTenor = \
			(double) \
			((int) (((dFXPaymentDates.at(1) - dFXPaymentDates.at(0)) * 12.0) + 0.5) / 12);
	}
	else if (n) // A single coupon in the set.
	{
		// Since there is only one FX Coupon, there is no way to calculate the
		// tenor with the information we have. Suppose then that the tenor is 3 months:
		CouponTenor = 0.25;
		ResetInAdvance = (dFXResetDates.at(0) < dFXPaymentDates.at(0) - CouponTenor);
		// Note: this is an approximation as the set may be 3 months before the end of
		// a 1-year coupon.
	}
	else throw("You should never come here"); 
	// as the function was called on an empty FX Coupon Set.

	// Find the first FX Coupon to value:
	// For each FX Coupon in the set:
	for (i = 0; (i < n); ++i)
	{
		// If reset is in advance:
		if (ResetInAdvance)
		{
			// If the notice is before the reset, then start valuation
			// at the current FX coupon:
			if (dNoticeDate <= dFXResetDates.at(i)) break;
		}
		// Otherwise, reset is in arrears:
		else
		{
			// If the notice is before what the "in advance" reset would be,
			// then start valuation at the current FX Coupon:
			if (dNoticeDate <= dFXResetDates.at(i) - CouponTenor) break;
		}
	}
	return i;
}
*/
double dAnalyticsTest(double dMarketDF0,double dMarketDF1,double dMeanReversion,double dStandardDeviation,double dDecay,double dDate,double dPaymentDate)
{
    double dR=-log(dMarketDF0)/dDate;
    double dTerm1=exp(-dMeanReversion*dPaymentDate)-exp(-dMeanReversion*dDate);
    double dTerm2=exp(2.*dMeanReversion*dDate)-1.;
    return exp(log(dMarketDF1/dMarketDF0)-dDecay*(-dR)-0.25*(pow(dStandardDeviation,2.)/pow(dMeanReversion,3.))*dTerm1*dTerm1*dTerm2);
}
// This function cannot be used for hybrid products
void LatticeSliceGreenCashCalcs(unsigned int uii, double dPaymentDate, double *dNorm,DKMaille<double> &dDecay)
{
    if(g_uiProductModelCode>3) throw("You cannot use LatticeSliceGreenCashCalcs for hybrid products");
    // Current lattice date
    double dDate=g_LATTICE_DATE.at(uii);
    // Relevant date for the fixed payment
    double dMarketDF1 = BaseDiscountInterpolate(dPaymentDate);
    // This will not work for time-dependent mean reversion.
    // Note from dk: time-dependent mean reversion models for fixed income have been seen to be unstable.
    switch(g_uiNumFactors)
    {
    case 1:
        {
            if(dDate!=dPaymentDate)
            {
                *dNorm=CalcGreensFunctionNormalisation(dMarketDF1,0,dDate,dPaymentDate,uii);
            }
            else
            {
                *dNorm=1.;
            }
            if(g_uiIsDK==0)
            {
                dDecay.at(0)=(1.-exp(-g_MR1*(dPaymentDate-dDate)))/g_MR1;
            }
            else
            {
                // Vladimir
                dDecay.at(0)=dLatticeDecayFunctional(g_MR1,uii,dPaymentDate,0);
            }
        }
        break;
    case 2:
        {
            if(dDate!=dPaymentDate)
            {
                *dNorm=CalcGreensFunctionNormalisation(dMarketDF1,0,dDate,dPaymentDate,uii);
            }
            else
            {
                *dNorm=1.;
            }
            dDecay.at(0)=(1.-exp(-g_MR1*(dPaymentDate-dDate)))/g_MR1;
            dDecay.at(1)=(1.-exp(-g_MR2*(dPaymentDate-dDate)))/g_MR2;
        }
        break;
    case 3:
        {
            if(dDate!=dPaymentDate)
            {
                *dNorm=CalcGreensFunctionNormalisation(dMarketDF1,0,dDate,dPaymentDate,uii);
            }
            else
            {
                *dNorm=1.;
            }
            dDecay.at(0)=(1.-exp(-g_MR1*(dPaymentDate-dDate)))/g_MR1;
            dDecay.at(1)=(1.-exp(-g_MR2*(dPaymentDate-dDate)))/g_MR2;
            dDecay.at(2)=(1.-exp(-g_MR3*(dPaymentDate-dDate)))/g_MR3;
        }
        break;
    }
}
void LatticeSliceGreenLIBORCalcs(unsigned int uii,double dStartDate,double dEndDate, double dPaymentDate, DKMaille<double> &dNorm,DKMaille2D<double> &dDecay)
{
    // Function used for single currency fixed income as well as hybrid products
    // Current lattice date
    double dDate=g_LATTICE_DATE.at(uii);
    // Relevant dates for the float payment
    // This incorporates flexile dates for payment, ie arrears/delay
    double dMarketDF1 = BaseDiscountInterpolate(dStartDate);
    double dMarketDF2 = BaseDiscountInterpolate(dEndDate);
    double dMarketDF3 = BaseDiscountInterpolate(dPaymentDate);
    // This will not work for time-dependent mean reversion.
    // Note from dk: time-dependent mean reversion models for fixed income are found to be non-stable.
    if(dDate!=dStartDate)
    {
        dNorm.at(0)=CalcGreensFunctionNormalisation(dMarketDF1,0,dDate,dStartDate,uii);
    }
    else
    {
        dNorm.at(0)=1.;
    }
    dNorm.at(1)=CalcGreensFunctionNormalisation(dMarketDF2,0,dDate,dEndDate,uii);
    if(dPaymentDate!=dDate)
    {
        if(dPaymentDate!=dEndDate)
        {
            dNorm.at(2)=CalcGreensFunctionNormalisation(dMarketDF3,0,dDate,dPaymentDate,uii);
        }
        else
        {
            dNorm.at(2)=dNorm.at(1);
        }
    }
    else
        dNorm.at(2)=1.;

    // hybrid product discounted by 1F- domestic q-model
    if(g_uiProductModelCode>3)
    {
        if(g_uiIsDK==0)
        {
            dDecay.at(0,0)=(1.-exp(-g_MR1*(dStartDate-dDate)))/g_MR1;
            dDecay.at(1,0)=(1.-exp(-g_MR1*(dEndDate-dDate)))/g_MR1;
            dDecay.at(2,0)=(1.-exp(-g_MR1*(dPaymentDate-dDate)))/g_MR1;
        }
        else
        {
            // Vladimir
            dDecay.at(0,0)=dLatticeDecayFunctional(g_MR1,uii,dStartDate,0);
            dDecay.at(1,0)=dLatticeDecayFunctional(g_MR1,uii,dEndDate,0);
            dDecay.at(2,0)=dLatticeDecayFunctional(g_MR1,uii,dPaymentDate,0);
        }
    }
    if(g_uiProductModelCode<=3)
    {
        switch(g_uiNumFactors)
        {
        case 3:
            {
                dDecay.at(0,2)=(1.-exp(-g_MR3*(dStartDate-dDate)))/g_MR3;
                dDecay.at(1,2)=(1.-exp(-g_MR3*(dEndDate-dDate)))/g_MR3;
                dDecay.at(2,2)=(1.-exp(-g_MR3*(dPaymentDate-dDate)))/g_MR3;
            }
        case 2:
            {
                dDecay.at(0,1)=(1.-exp(-g_MR2*(dStartDate-dDate)))/g_MR2;
                dDecay.at(1,1)=(1.-exp(-g_MR2*(dEndDate-dDate)))/g_MR2;
                dDecay.at(2,1)=(1.-exp(-g_MR2*(dPaymentDate-dDate)))/g_MR2;
            }
        case 1:
            {
                if(g_uiIsDK==0)
                {
                    dDecay.at(0,0)=(1.-exp(-g_MR1*(dStartDate-dDate)))/g_MR1;
                    dDecay.at(1,0)=(1.-exp(-g_MR1*(dEndDate-dDate)))/g_MR1;
                    dDecay.at(2,0)=(1.-exp(-g_MR1*(dPaymentDate-dDate)))/g_MR1;
                }
                else
                {
                    // Vladimir
                    dDecay.at(0,0)=dLatticeDecayFunctional(g_MR1,uii,dStartDate,0);
                    dDecay.at(1,0)=dLatticeDecayFunctional(g_MR1,uii,dEndDate,0);
                    dDecay.at(2,0)=dLatticeDecayFunctional(g_MR1,uii,dPaymentDate,0);
                }
            }
        }
    }
}

void LatticeSliceGreenFXCalcs(unsigned int uii,double dDate,double *dNormBase,double *dNormForeign,double *dDecayDomestic,double *dDecayForeign)
{
    double dMarketDF = BaseDiscountInterpolate(dDate);
    double dForeignDF = ForeignDiscountInterpolate(dDate);
    // This will not work for time-dependent mean reversion.
    // Note from dk: time-dependent mean reversion models for fixed income are found to be non-stable.
    *dNormBase=CalcGreensFunctionNormalisation(dMarketDF,0,g_LATTICE_DATE.at(uii),dDate,uii);
    *dNormForeign=CalcGreensFunctionNormalisation(dForeignDF,1,g_LATTICE_DATE.at(uii),dDate,uii);

    if(g_uiIsDK==0)
    {
        *dDecayDomestic=(1.-exp(-g_MR1*(dDate-g_LATTICE_DATE.at(uii))))/g_MR1;
        *dDecayForeign=(1.-exp(-g_MR2*(dDate-g_LATTICE_DATE.at(uii))))/g_MR2;
    }
    else
    {
        *dDecayDomestic=dLatticeDecayFunctional(g_MR1,uii,dDate,0);
        *dDecayForeign=dLatticeDecayFunctional(g_MR2,uii,dDate,0);
    }

}
void LatticeFillFX(double dSpotFX)
{
    // String Model for cash diffusion ONLY
    DKMaille<double> dGRateCorrection(g_uiNumSlices);
    for(unsigned ui=0;ui<g_uiNumSlices;ui++)
    {
        int iFlag=iIsNoticeDate(ui,g_dNoticeDates);
        dGRateCorrection.at(ui)=0.;
        // only do this for notice dates

        if(iFlag>-1&&g_dADPLimit!=1.0 || ui==0) /// JMP to allow analytical fx option prices at 0
        {
            // Size the slice
            unsigned int uiNumberOfStates;
            uiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(ui)
                             *g_LATTICE_NUMBER_OF_STATES2.at(ui)*g_LATTICE_NUMBER_OF_STATES3.at(ui);
            // Calc the string tension
            dGRateCorrection.at(ui)=0.;
            double dNorm=0.;
            for(unsigned uk=0;uk<uiNumberOfStates;uk++)
            {
                if(g_LATTICE_NODE_EXISTENCE.at(ui).at(uk)==true)
                {
                    dNorm+=g_LATTICE_ARROW_DEBREU.at(ui).at(uk)*exp(g_LATTICE_G_RATE_FX.at(ui).at(uk));
                }
            }
            // if normal
            dGRateCorrection.at(ui)=log(g_LATTICE_ZCB_PRICE_FOREIGN.at(ui)/dNorm);
        }
    }
    for(ui=0;ui<g_uiNumSlices;ui++)
    {
        int iFlag=iIsNoticeDate(ui,g_dNoticeDates);
        // only do this for notice dates

        if(iFlag>-1 || ui == 0) /// JMP to allow analytical fx option prices at 0
        {
            // Size the slice
            unsigned int uiNumberOfStates;
            uiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(ui)
                             *g_LATTICE_NUMBER_OF_STATES2.at(ui)*g_LATTICE_NUMBER_OF_STATES3.at(ui);
            for(unsigned uk=0;uk<uiNumberOfStates;uk++)
            {
                if(g_LATTICE_NODE_EXISTENCE.at(ui).at(uk)==true)
                {
                    g_LATTICE_SHORT_RATE_FX.at(ui).at(uk)=dSpotFX*exp(g_LATTICE_G_RATE_FX.at(ui).at(uk)+dGRateCorrection.at(ui));
                }
            }
        }
    }
} // LatticeFillFX

void LatticeApplySmoothing(unsigned int uii,
                           unsigned int uiSmoothingOrder,
                           DKMaille< DKMaille<double> > &dOption,
                           DKMaille< DKMaille<double> > &dFunctionToSmooth)
{
    switch(g_uiNumFactors)
    {
    case 3:
        {
            DKMaille<long> currentNodeStates(g_uiNumFactors);
            int iAdjustmentIndex=0;
            double Adjustment=0.;

            unsigned int uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uii)*g_LATTICE_NUMBER_OF_STATES2.at(uii)*g_LATTICE_NUMBER_OF_STATES3.at(uii);
            DKMaille<double> dRawOption(uiiNumberOfStates, 0.0);

            for(int i1=g_LATTICE_MIN_STATE1.at(uii);i1<=g_LATTICE_MAX_STATE1.at(uii);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uii);i2<=g_LATTICE_MAX_STATE2.at(uii);i2++)
                {
                    for(int i3=g_LATTICE_MIN_STATE3.at(uii);i3<=g_LATTICE_MAX_STATE3.at(uii);i3++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiStateIndex=LatticeIndexConvert(uii,currentNodeStates);
                        //if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uiStateIndex))
                        //    dRawOption.at(uiStateIndex)=dOption.at(uii).at(uiStateIndex);
						if( (uiStateIndex < g_LATTICE_NODE_EXISTENCE.at(uii).entries()) && g_LATTICE_NODE_EXISTENCE.at(uii).at(uiStateIndex) )
						{
							dRawOption.at(uiStateIndex)=dOption.at(uii).at(uiStateIndex);
						}
                    }
                }
            }

            unsigned int uiN1=g_LATTICE_NUMBER_OF_STATES1.at(uii);
            DKMaille<double> dOption1F(uiN1, 0.0);
            DKMaille<double> dFunctionToSmooth1F(uiN1, 0.0);

            for(int i3=g_LATTICE_MIN_STATE3.at(uii);i3<=g_LATTICE_MAX_STATE3.at(uii);i3++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uii);i2<=g_LATTICE_MAX_STATE2.at(uii);i2++)
                {
                    for(int i1=g_LATTICE_MIN_STATE1.at(uii);i1<=g_LATTICE_MAX_STATE1.at(uii);i1++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiStateIndex=LatticeIndexConvert(uii,currentNodeStates);
						dOption1F.at(i1-g_LATTICE_MIN_STATE1.at(uii))=dRawOption.at(uiStateIndex);
						dFunctionToSmooth1F.at(i1-g_LATTICE_MIN_STATE1.at(uii))=dFunctionToSmooth.at(uii).at(uiStateIndex);
                    }
                    SmoothingAlgorithm(iAdjustmentIndex,
                                       Adjustment,
                                       dFunctionToSmooth1F,
                                       dOption1F);
                    currentNodeStates.at(0)=iAdjustmentIndex+g_LATTICE_MIN_STATE1.at(uii);
                    unsigned int uiStateIndex=LatticeIndexConvert(uii,currentNodeStates);
					if ( (uiStateIndex < g_LATTICE_NODE_EXISTENCE.at(uii).entries()) && g_LATTICE_NODE_EXISTENCE.at(uii).at(uiStateIndex) )
					{
						dOption.at(uii).at(uiStateIndex)+=Adjustment/3.;
					}
                }
            }

            unsigned int uiN2=g_LATTICE_NUMBER_OF_STATES2.at(uii);
            DKMaille<double> dOption2F(uiN2, 0.0);
            DKMaille<double> dFunctionToSmooth2F(uiN2, 0.0);

            for(i3=g_LATTICE_MIN_STATE3.at(uii);i3<=g_LATTICE_MAX_STATE3.at(uii);i3++)
            {
                for(i1=g_LATTICE_MIN_STATE1.at(uii);i1<=g_LATTICE_MAX_STATE1.at(uii);i1++)
                {
                    unsigned int uiCountStates=0;
                    for(int i2=g_LATTICE_MIN_STATE2.at(uii);i2<=g_LATTICE_MAX_STATE2.at(uii);i2++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiStateIndex=LatticeIndexConvert(uii,currentNodeStates);
	                    dOption2F.at(i2-g_LATTICE_MIN_STATE2.at(uii))=dRawOption.at(uiStateIndex);
		                dFunctionToSmooth2F.at(i2-g_LATTICE_MIN_STATE2.at(uii))=dFunctionToSmooth.at(uii).at(uiStateIndex);
                    }
                    SmoothingAlgorithm(iAdjustmentIndex,
                                       Adjustment,
                                       dFunctionToSmooth2F,
                                       dOption2F);
                    currentNodeStates.at(1)=iAdjustmentIndex+g_LATTICE_MIN_STATE2.at(uii);
                    unsigned int uiStateIndex=LatticeIndexConvert(uii,currentNodeStates);
                    //dOption.at(uii).at(LatticeIndexConvert(uii,currentNodeStates))+=Adjustment/3.;
					if ( (uiStateIndex < g_LATTICE_NODE_EXISTENCE.at(uii).entries()) && g_LATTICE_NODE_EXISTENCE.at(uii).at(uiStateIndex) )
					{
						dOption.at(uii).at(uiStateIndex)+=Adjustment/3.;
					}
                }
            }

            unsigned int uiN3=g_LATTICE_NUMBER_OF_STATES3.at(uii);
            DKMaille<double> dOption3F(uiN3, 0.0);
            DKMaille<double> dFunctionToSmooth3F(uiN3, 0.0);

            for(i1=g_LATTICE_MIN_STATE1.at(uii);i1<=g_LATTICE_MAX_STATE1.at(uii);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uii);i2<=g_LATTICE_MAX_STATE2.at(uii);i2++)
                {
                    unsigned int uiCountStates=0;
                    for(int i3=g_LATTICE_MIN_STATE3.at(uii);i3<=g_LATTICE_MAX_STATE3.at(uii);i3++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiStateIndex=LatticeIndexConvert(uii,currentNodeStates);
						dOption3F.at(i3-g_LATTICE_MIN_STATE3.at(uii))=dRawOption.at(uiStateIndex);
						dFunctionToSmooth3F.at(i3-g_LATTICE_MIN_STATE3.at(uii))=dFunctionToSmooth.at(uii).at(uiStateIndex);
                    }
                    SmoothingAlgorithm(iAdjustmentIndex,
                                       Adjustment,
                                       dFunctionToSmooth3F,
                                       dOption3F);
                    currentNodeStates.at(2)=iAdjustmentIndex+g_LATTICE_MIN_STATE3.at(uii);
                    unsigned int uiStateIndex=LatticeIndexConvert(uii,currentNodeStates);
                    //dOption.at(uii).at(LatticeIndexConvert(uii,currentNodeStates))+=Adjustment/3.;
					if ( (uiStateIndex < g_LATTICE_NODE_EXISTENCE.at(uii).entries()) && g_LATTICE_NODE_EXISTENCE.at(uii).at(uiStateIndex) )
					{
						dOption.at(uii).at(uiStateIndex)+=Adjustment/3.;
					}
                }
            }
        }
        break;
    case 2:
        {

            unsigned int uiN1=g_LATTICE_NUMBER_OF_STATES1.at(uii);
            DKMaille<long> currentNodeStates(g_uiNumFactors);
            DKMaille<double> dOption1F(uiN1);
            DKMaille<double> dFunctionToSmooth1F(uiN1);
            int iAdjustmentIndex=0;
            double Adjustment=0.;

            unsigned int uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uii)*g_LATTICE_NUMBER_OF_STATES2.at(uii);
            DKMaille<double> dRawOption(uiiNumberOfStates);

            for(int i1=g_LATTICE_MIN_STATE1.at(uii);i1<=g_LATTICE_MAX_STATE1.at(uii);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uii);i2<=g_LATTICE_MAX_STATE2.at(uii);i2++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiStateIndex=LatticeIndexConvert(uii,currentNodeStates);
                    if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uiStateIndex))
                        dRawOption.at(uiStateIndex)=dOption.at(uii).at(uiStateIndex);
                }
            }
            for(int i2=g_LATTICE_MIN_STATE2.at(uii);i2<=g_LATTICE_MAX_STATE2.at(uii);i2++)
            {
                for(int i1=g_LATTICE_MIN_STATE1.at(uii);i1<=g_LATTICE_MAX_STATE1.at(uii);i1++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiStateIndex=LatticeIndexConvert(uii,currentNodeStates);
					dOption1F.at(i1-g_LATTICE_MIN_STATE1.at(uii))=dRawOption.at(uiStateIndex);
					dFunctionToSmooth1F.at(i1-g_LATTICE_MIN_STATE1.at(uii))=dFunctionToSmooth.at(uii).at(uiStateIndex);
                }
                SmoothingAlgorithm(iAdjustmentIndex,
                                   Adjustment,
                                   dFunctionToSmooth1F,
                                   dOption1F);
                currentNodeStates.at(0)=iAdjustmentIndex+g_LATTICE_MIN_STATE1.at(uii);
                unsigned int uiStateIndex=LatticeIndexConvert(uii,currentNodeStates);
				dOption.at(uii).at(uiStateIndex)+=Adjustment/2.;
            }

            unsigned int uiN2=g_LATTICE_NUMBER_OF_STATES2.at(uii);
            DKMaille<double> dOption2F(uiN2);
            DKMaille<double> dFunctionToSmooth2F(uiN2);

            for(i1=g_LATTICE_MIN_STATE1.at(uii);i1<=g_LATTICE_MAX_STATE1.at(uii);i1++)
            {
                unsigned int uiCountStates=0;
                for(int i2=g_LATTICE_MIN_STATE2.at(uii);i2<=g_LATTICE_MAX_STATE2.at(uii);i2++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiStateIndex=LatticeIndexConvert(uii,currentNodeStates);
					dOption2F.at(i2-g_LATTICE_MIN_STATE2.at(uii))=dRawOption.at(uiStateIndex);
					dFunctionToSmooth2F.at(i2-g_LATTICE_MIN_STATE2.at(uii))=dFunctionToSmooth.at(uii).at(uiStateIndex);
                }
                SmoothingAlgorithm(iAdjustmentIndex,
                                   Adjustment,
                                   dFunctionToSmooth2F,
                                   dOption2F);
                currentNodeStates.at(1)=iAdjustmentIndex+g_LATTICE_MIN_STATE2.at(uii);
                unsigned int uiStateIndex=LatticeIndexConvert(uii,currentNodeStates);
				dOption.at(uii).at(LatticeIndexConvert(uii,currentNodeStates))+=Adjustment/2.;
            }
        }
        break;
    case 1:
        {
            unsigned int uiN1=g_LATTICE_NUMBER_OF_STATES1.at(uii);
            DKMaille<long> currentNodeStates(g_uiNumFactors);
            DKMaille<double> dOption1F(uiN1);
            DKMaille<double> dFunctionToSmooth1F(uiN1);
            int iAdjustmentIndex=0;
            double Adjustment=0.;

            unsigned int uiCountStates=0;
            for(int i1=g_LATTICE_MIN_STATE1.at(uii);i1<=g_LATTICE_MAX_STATE1.at(uii);i1++)
            {
                currentNodeStates.at(0)=i1;
                if(g_LATTICE_NODE_EXISTENCE.at(uii).at(LatticeIndexConvert(uii,currentNodeStates)))
                    uiCountStates++;
                dOption1F.at(i1-g_LATTICE_MIN_STATE1.at(uii))=dOption.at(uii).at(LatticeIndexConvert(uii,currentNodeStates));
                dFunctionToSmooth1F.at(i1-g_LATTICE_MIN_STATE1.at(uii))=dFunctionToSmooth.at(uii).at(LatticeIndexConvert(uii,currentNodeStates));
            }
            SmoothingAlgorithm(iAdjustmentIndex,
                               Adjustment,
                               dFunctionToSmooth1F,
                               dOption1F);
            currentNodeStates.at(0)=iAdjustmentIndex+g_LATTICE_MIN_STATE1.at(uii);
            dOption.at(uii).at(LatticeIndexConvert(uii,currentNodeStates))+=Adjustment;
        }
        break;
    } // switch on num factors
} // LatticeApplySmoothing()

void LatticeGetDiscountedSecurity(unsigned int uiSlice,DKMaille< DKMaille<double> > &dSecurity)
{
    unsigned int uiiNumberOfStates;
    unsigned int uiNumberOfBranches;

    switch(g_uiNumFactors)
    {
    case 3:
        {
            uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uiSlice)*g_LATTICE_NUMBER_OF_STATES2.at(uiSlice)*g_LATTICE_NUMBER_OF_STATES3.at(uiSlice);
            uiNumberOfBranches=27;
        }
        break;
    case 2:
        {
            uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uiSlice)*g_LATTICE_NUMBER_OF_STATES2.at(uiSlice);
            uiNumberOfBranches=9;
        }
        break;
    case 1:
        {
            uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uiSlice);
            uiNumberOfBranches=3;
        }
    }
    for(unsigned int uk=0;uk<uiiNumberOfStates;uk++)
    {
        if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uk)==true)
        {
            for(unsigned int um=0;um<uiNumberOfBranches;um++)
            {
                dSecurity.at(uiSlice).at(uk)+=(1./(double)g_dConvert)*(double)g_LATTICE_GREEN.at(uiSlice,um).at(uk)
                                              *dSecurity.at(uiSlice+1).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uk));
            }
        } // if node exists
    } // loop on states
}

DKMaille<double> &CreateNewStripDate(void)
{
    static DKMaille<double> dNewStripDates;
    dNewStripDates.resize(g_LATTICE_DATE.entries()+g_ppy);

    // les dates avant la derniere notice date sont les memes que dans le lattice
    for(unsigned int ui=0;ui<g_LATTICE_DATE.entries();ui++)
    {
        dNewStripDates.at(ui)=g_LATTICE_DATE.at(ui);
    }

    // la derniere date est l'expiry date
    dNewStripDates.at(dNewStripDates.entries()-1)=g_dLastDate;

    // Si le nombre de slice entre la derniere notice et l'expiration est > 1, on calcule les dates correspondantes!
    for(ui=0;ui<g_ppy-1;ui++)
    {
        dNewStripDates.at(g_LATTICE_DATE.entries()+ui)=dNewStripDates.at(g_LATTICE_DATE.entries()-1)
                +(ui+1)*(g_dLastDate-dNewStripDates.at(g_LATTICE_DATE.entries()-1))/g_ppy;
    }
    return dNewStripDates;
}


void PropagateTotalCoupon(unsigned int uiSlice,DKMaille< DKMaille<double> > &dTotalCoupon)
{
    unsigned int uiNextSlice=uiSlice+1;
    double dt = g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice);
    DKMaille<long> currentNodeStates(g_uiNumFactors);
    DKMaille<long> targetNodeStates(g_uiNumFactors);
    unsigned int uiiNumberOfStates=0;
    if(g_uiNumFactors==3)
    {
        uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice)*
                          g_LATTICE_NUMBER_OF_STATES2.at(uiNextSlice)*
                          g_LATTICE_NUMBER_OF_STATES3.at(uiNextSlice);
    }
    if(g_uiNumFactors==2)
    {
        uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice)*
                          g_LATTICE_NUMBER_OF_STATES2.at(uiNextSlice);
    }
    if(g_uiNumFactors==1)
    {
        uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uiNextSlice);
    }
    DKMaille<long> iBranch(g_uiNumFactors);
    DKMaille<double> dConditionalGreen;
    dTotalCoupon.at(uiNextSlice).resize(uiiNumberOfStates);
    dConditionalGreen.resize(uiiNumberOfStates);
    double dProbability=0.;
    for(unsigned int ug=0;ug<uiiNumberOfStates;ug++)
        dConditionalGreen.at(ug)=0.;
    for(ug=0;ug<uiiNumberOfStates;ug++)
        dTotalCoupon.at(uiNextSlice).at(ug)=0.;
    switch(g_uiNumFactors)
    {
    case 1:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                currentNodeStates.at(0)=i1;
                unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                // Get the states from the next slice that we will need to use in connecting
                if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                {
                    double dDiscountConnection=LatticeDiscount(uiSlice,uiIndex);
                    for(unsigned int um=0;um<3;um++)
                    {
                        if(g_LATTICE_NODE_EXISTENCE.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex)))
                        {
                            dProbability=(1./(double)g_dConvert)*(double)g_LATTICE_GREEN.at(uiSlice,um).at(uiIndex)/dDiscountConnection;
                            dConditionalGreen.at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))+=dProbability;
                            dTotalCoupon.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))+=dProbability*dTotalCoupon.at(uiSlice).at(uiIndex);
                        }
                    }
                }
            } // while loop over nodes

            // Renormalisation of probability
            for(i1=g_LATTICE_MIN_STATE1.at(uiNextSlice);i1<=g_LATTICE_MAX_STATE1.at(uiNextSlice);i1++)
            {
                currentNodeStates.at(0)=i1;
                unsigned int uiIndex=LatticeIndexConvert(uiNextSlice,currentNodeStates);
                // Get the states from the next slice that we will need to use in connecting
                if(g_LATTICE_NODE_EXISTENCE.at(uiNextSlice).at(uiIndex))
                {
                    if(dConditionalGreen.at(uiIndex)!=0.) dTotalCoupon.at(uiNextSlice).at(uiIndex)/=dConditionalGreen.at(uiIndex);
                }
            } // while loop over nodes
        }
        break;
    case 2:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                    // Get the states from the next slice that we will need to use in connecting
                    if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                    {
                        double dDiscountConnection=LatticeDiscount(uiSlice,uiIndex);
                        for(unsigned int um=0;um<9;um++)
                        {
                            if(g_LATTICE_NODE_EXISTENCE.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex)))
                            {
                                dProbability=(1./(double)g_dConvert)*(double)g_LATTICE_GREEN.at(uiSlice,um).at(uiIndex)/dDiscountConnection;
                                dConditionalGreen.at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))+=dProbability;
                                dTotalCoupon.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))+=dProbability*dTotalCoupon.at(uiSlice).at(uiIndex);
                            }
                        }
                    }
                }
            } // while loop over nodes

            // Renormalisation of probability
            for(i1=g_LATTICE_MIN_STATE1.at(uiNextSlice);i1<=g_LATTICE_MAX_STATE1.at(uiNextSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiNextSlice);i2<=g_LATTICE_MAX_STATE2.at(uiNextSlice);i2++)
                {
                    currentNodeStates.at(0)=i1;
                    currentNodeStates.at(1)=i2;
                    unsigned int uiIndex=LatticeIndexConvert(uiNextSlice,currentNodeStates);
                    // Get the states from the next slice that we will need to use in connecting
                    if(g_LATTICE_NODE_EXISTENCE.at(uiNextSlice).at(uiIndex))
                    {
                        if(dConditionalGreen.at(uiIndex)!=0.) dTotalCoupon.at(uiNextSlice).at(uiIndex)/=dConditionalGreen.at(uiIndex);
                    }
                }
            } // while loop over nodes
        }
        break;
    case 3:
        {
            for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                {
                    for(int i3=g_LATTICE_MIN_STATE3.at(uiSlice);i3<=g_LATTICE_MAX_STATE3.at(uiSlice);i3++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                        // Get the states from the next slice that we will need to use in connecting
                        if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                        {
                            double dDiscountConnection=LatticeDiscount(uiSlice,uiIndex);
                            for(unsigned int um=0;um<27;um++)
                            {
                                if(g_LATTICE_NODE_EXISTENCE.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex)))
                                {
                                    dProbability=(1./(double)g_dConvert)*(double)g_LATTICE_GREEN.at(uiSlice,um).at(uiIndex)/dDiscountConnection;
                                    dConditionalGreen.at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))+=dProbability;
                                    dTotalCoupon.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))+=dProbability*dTotalCoupon.at(uiSlice).at(uiIndex);
                                }
                            }
                        }
                    }
                }
            } // while loop over nodes

            // Renormalisation of probability
            for(i1=g_LATTICE_MIN_STATE1.at(uiNextSlice);i1<=g_LATTICE_MAX_STATE1.at(uiNextSlice);i1++)
            {
                for(int i2=g_LATTICE_MIN_STATE2.at(uiNextSlice);i2<=g_LATTICE_MAX_STATE2.at(uiNextSlice);i2++)
                {
                    for(int i3=g_LATTICE_MIN_STATE3.at(uiNextSlice);i3<=g_LATTICE_MAX_STATE3.at(uiNextSlice);i3++)
                    {
                        currentNodeStates.at(0)=i1;
                        currentNodeStates.at(1)=i2;
                        currentNodeStates.at(2)=i3;
                        unsigned int uiIndex=LatticeIndexConvert(uiNextSlice,currentNodeStates);
                        // Get the states from the next slice that we will need to use in connecting
                        if(g_LATTICE_NODE_EXISTENCE.at(uiNextSlice).at(uiIndex))
                        {
                            if(dConditionalGreen.at(uiIndex)!=0.) dTotalCoupon.at(uiNextSlice).at(uiIndex)/=dConditionalGreen.at(uiIndex);
                        }
                    }
                }
            } // while loop over nodes
        }
        break;
    }
}

DKMaille<double> Lattice_HWVFDK_3F_Calc(DKMaille2D<double> dBoosterData)
{
    if(g_dDeltaFlag==1.) g_uiDeltaTick+=1;
    // Tests on the resulting Green's function
    double dLatticeDF =0.;
    double dLatticeDFForeign =0.;
    double dMarketDF=0.;
    double dMarketDFForeign=0.;
    double dGreenOption=0.;
    double dBaseInterestRate=0;
    double dLatticeBaseInterestRate=0;
    double dForeignInterestRate=0;
    double dLatticeForeignInterestRate=0;

    NbNodesPerSliceGlobalExtraction.Resize(g_uiNumSlices);

    for(unsigned ui=0;ui<g_uiNumSlices;ui++)
    {
        NbNodesPerSliceGlobalExtraction[ui]=0;

        dLatticeDF = 0. ;
        if(g_uiProductModelCode>3) dLatticeDFForeign = 0.;
        // Size the slice
        unsigned int uiNumberOfStates;
        if(g_uiNumFactors==3) uiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(ui)
                                                   *g_LATTICE_NUMBER_OF_STATES2.at(ui)*g_LATTICE_NUMBER_OF_STATES3.at(ui);
        if(g_uiNumFactors==2) uiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(ui)
                                                   *g_LATTICE_NUMBER_OF_STATES2.at(ui);
        if(g_uiNumFactors==1) uiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(ui);
        int iFlag=iIsNoticeDate(ui,g_dNoticeDates);
        if(iFlag>-1)
        {
            unsigned int uiRealStates=0;
            dMarketDF = BaseDiscountInterpolate(g_LATTICE_DATE.at(ui));
            if(g_uiProductModelCode>3) dMarketDFForeign = ForeignDiscountInterpolate(g_LATTICE_DATE.at(ui));
            for(unsigned uk=0;uk<uiNumberOfStates;uk++)
            {
                if(g_LATTICE_NODE_EXISTENCE.at(ui).at(uk)==true)
                {
                    ++NbNodesPerSliceGlobalExtraction[ui];

                    dLatticeDF+=g_LATTICE_ARROW_DEBREU.at(ui).at(uk);
                    uiRealStates++;
                }
            }
            dBaseInterestRate=log(dMarketDF)/(-g_LATTICE_DATE.at(ui));
            dLatticeBaseInterestRate=log(dLatticeDF)/(-g_LATTICE_DATE.at(ui));
        }
        else
        {
            for(unsigned uk=0;uk<uiNumberOfStates;uk++)
            {
                if(g_LATTICE_NODE_EXISTENCE.at(ui).at(uk))
                    ++NbNodesPerSliceGlobalExtraction[ui];
            }
        }
    }
    if(g_dType<3.)
    {
        DKMaille< DKMaille<double> > dOption;
        dOption.resize(g_uiNumSlices);
        double dBermHolder;
        int uii=0;
        if(g_uiProductModelCode==5)
        {
            for(uii=g_uiNumSlices-1;uii>=0;uii--)
            {
                // Size the slice
                unsigned int uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uii)*g_LATTICE_NUMBER_OF_STATES2.at(uii)
                                               *g_LATTICE_NUMBER_OF_STATES3.at(uii);
                // initialise
                dOption.at(uii).resize(uiiNumberOfStates);
                for(unsigned int uf=0;uf<uiiNumberOfStates;uf++) dOption.at(uii).at(uf)=0.;
                // European Option
                // Case where lattice extends to  european option expiry
                if((uii==g_uiNumSlices-1)&&(g_LATTICE_DATE.at(g_uiNumSlices-1)==g_dOptionExpiry))
                {
                    for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                    {
                        if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                        {
                            if(g_dType==1.0||g_dType==2.0)
                            {
                                // Rolling Strike for American option
                                if(g_dType==2.0)
                                {
                                    // g_dStrike=g_dSpotFX*g_LATTICE_ZCB_PRICE_FOREIGN.at(uii)/g_LATTICE_ZCB_PRICE.at(uii);
                                }
                                if(g_LATTICE_SHORT_RATE_FX.at(uii).at(uk)>g_dStrike)
                                {
                                    dOption.at(uii).at(uk)=(g_LATTICE_SHORT_RATE_FX.at(uii).at(uk)-g_dStrike);
                                }
                                else
                                {
                                    dOption.at(uii).at(uk)=0.;
                                }
                            }
                            if(g_dType==1.5||g_dType==2.5)
                            {
                                // Rolling Strike for American option
                                if(g_dType==2.5)
                                {
                                    // g_dStrike=g_dSpotFX*g_LATTICE_ZCB_PRICE_FOREIGN.at(uii)/g_LATTICE_ZCB_PRICE.at(uii);
                                }
                                if(g_LATTICE_SHORT_RATE_FX.at(uii).at(uk)<g_dStrike)
                                {
                                    dOption.at(uii).at(uk)=-g_LATTICE_SHORT_RATE_FX.at(uii).at(uk)+g_dStrike;
                                }
                                else
                                {
                                    dOption.at(uii).at(uk)=0.;
                                }
                            }
                        }
                    }
                }
                else if((uii==g_uiNumSlices-1)&&(g_LATTICE_DATE.at(g_uiNumSlices-1)!=g_dOptionExpiry))
                {
                    // Introduction here of Green's Function Analytics for fast discounting of bonds and fx options
                    if(g_dADPLimit!=1.0)
                    {
                        // Preparation for node calculation
                        double dMarketDF = BaseDiscountInterpolate(g_dOptionExpiry);
                        double dMarketDFForeign = ForeignDiscountInterpolate(g_dOptionExpiry);
                        // This will not work for time-dependent mean reversion.
                        // Note from dk: time-dependent mean reversion models for fixed income are found to be non-stable.
                        DKMaille<double> dNormBase(2); DKMaille<double> dNormForeign(2); DKMaille<double> dDecayDomestic(2); DKMaille<double> dDecayForeign(2);
                        dNormBase.at(0)=CalcGreensFunctionNormalisation(dMarketDF,0,g_LATTICE_DATE.at(uii),g_dOptionExpiry,uii);
                        dNormBase.at(1)=dNormBase.at(0);
                        dNormForeign.at(0)=CalcGreensFunctionNormalisation(dMarketDFForeign,1,g_LATTICE_DATE.at(uii),g_dOptionExpiry,uii);
                        dNormForeign.at(1)=dNormForeign.at(0);
                        if(g_uiIsDK==0)
                        {
                            dDecayDomestic.at(0)=(1.-exp(-g_MR1*(g_dOptionExpiry-g_LATTICE_DATE.at(uii))))/g_MR1;
                            dDecayDomestic.at(1)=dDecayDomestic.at(0);
                            dDecayForeign.at(0)=(1.-exp(-g_MR2*(g_dOptionExpiry-g_LATTICE_DATE.at(uii))))/g_MR2;
                            dDecayForeign.at(1)=dDecayForeign.at(0);
                        }
                        else
                        {
                            // Vladimir's Approximation
                            dDecayDomestic.at(0)=dLatticeDecayFunctional(g_MR1,uii,g_dOptionExpiry,0);
                            dDecayDomestic.at(1)=dDecayDomestic.at(0);
                            dDecayForeign.at(0)=dLatticeDecayFunctional(g_MR2,uii,g_dOptionExpiry,1);
                            dDecayForeign.at(1)=dDecayForeign.at(0);
                        }

                        // Set the payment date to the option expiry date
                        double dPaymentDate=g_dOptionExpiry;
                        // Initialise local FX volatility calculation
                        g_uiVolIsCalculated.at(uii)=0;
                        for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                        {
                            if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                            {
                                // The function returns the value of the underlying.
                                // In that case it is a simple FX Option cashflow
                                // Rolling Strike for American option
                                if(g_dType==2.0||g_dType==2.5)
                                {
                                    // g_dStrike=g_dSpotFX*g_LATTICE_ZCB_PRICE_FOREIGN.at(uii)/g_LATTICE_ZCB_PRICE.at(uii);
                                }
                                if(g_dType==1.0||g_dType==2.0)
                                {
                                    dOption.at(uii).at(uk)=LatticeFXOptionAnalytics(uii,
                                                           uk,
                                                           dDecayDomestic,
                                                           dDecayForeign,
                                                           g_dOptionExpiry,
                                                           dPaymentDate,
                                                           g_dStrike,
                                                           1000000.,
                                                           dNormBase,
                                                           dNormForeign,
                                                           0.0);
                                }
                                if(g_dType==1.5||g_dType==2.5)
                                {
                                    dOption.at(uii).at(uk)=LatticeFXOptionAnalytics(uii,
                                                           uk,
                                                           dDecayDomestic,
                                                           dDecayForeign,
                                                           g_dOptionExpiry,
                                                           dPaymentDate,
                                                           g_dStrike,
                                                           1000000.,
                                                           dNormBase,
                                                           dNormForeign,
                                                           1.0);
                                }
                            }
                        }
                    }
                    else
                        throw("You are missing the Green's function for such calculation");
                }
                else
                {
                    for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                    {
                        if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                        {
                            for(unsigned int um=0;um<27;um++)
                            {
                                dOption.at(uii).at(uk)+=(1./(double)g_dConvert)*(double)g_LATTICE_GREEN.at(uii,um).at(uk)
                                                        *dOption.at(uii+1).at(g_LATTICE_CONNECTION.at(uii,um).at(uk));
                            }
                            // Calc the bermudan holder
                            if(g_dType==2.)
                            {
                                // g_dStrike=g_dSpotFX*g_LATTICE_ZCB_PRICE_FOREIGN.at(uii)/g_LATTICE_ZCB_PRICE.at(uii);
                                if(g_LATTICE_SHORT_RATE_FX.at(uii).at(uk)>g_dStrike)
                                {
                                    dBermHolder=(g_LATTICE_SHORT_RATE_FX.at(uii).at(uk)-g_dStrike);
                                }
                                else
                                {
                                    dBermHolder=0.;
                                }
                                // Calc the bermudan option
                                if(dBermHolder>=dOption.at(uii).at(uk))
                                    dOption.at(uii).at(uk)=dBermHolder;
                            } // if american option call
                            if(g_dType==2.5)
                            {
                                // g_dStrike=g_dSpotFX*g_LATTICE_ZCB_PRICE_FOREIGN.at(uii)/g_LATTICE_ZCB_PRICE.at(uii);
                                if(g_LATTICE_SHORT_RATE_FX.at(uii).at(uk)<g_dStrike)
                                {
                                    dBermHolder=(-g_LATTICE_SHORT_RATE_FX.at(uii).at(uk)+g_dStrike);
                                }
                                else
                                {
                                    dBermHolder=0.;
                                }
                                // Calc the bermudan option
                                if(dBermHolder>=dOption.at(uii).at(uk))
                                    dOption.at(uii).at(uk)=dBermHolder;
                            } // if american option put
                        } // if node exists
                    } // loop on states
                } // if on backwardation
                if(g_dDeltaFlag==1.&&g_uiDeltaTick==1) LatticeFreeGreen(uii);
                if(g_dDeltaFlag==0.&&g_uiDeltaTick==-1) LatticeFreeGreen(uii);
                if(uii!=g_uiNumSlices-1) dOption.at(uii+1).clear();
            } // loop on slices
        } // if g_uiProductModelCode == 5
        static DKMaille<double> ret(3);
        if(g_uiProductModelCode<4) ret.at(0)=dLatticeDF;
        else if(g_uiProductModelCode==4) ret.at(0)=123;
        else ret.at(0)=dOption.at(0).at(0);
        ret.at(1)=0.0;
        ret.at(2)=0.0;
        return ret;
    } // if on g_dType == 1,2
    else if(g_dType>=3.)
    {
        // callable fx booster pricing
        // Get the product data first
        DKMaille<double> dFXCouponResetDates(dBoosterData.rows()); DKMaille<double> dFXCouponPaymentDates(dBoosterData.rows());
        DKMaille<double> dIndexStartDates(dBoosterData.rows());	DKMaille<double> dIndexEndDates(dBoosterData.rows());	DKMaille<double> dIndexPaymentDates(dBoosterData.rows());
        DKMaille<double> dFundingSpread(dBoosterData.rows());	DKMaille<double> dFXNotionalMultiplier(dBoosterData.rows());	DKMaille<double> dFXStrike(dBoosterData.rows());
        DKMaille<double> dFXCap(dBoosterData.rows()); DKMaille<double> dFXCouponAccrualBasis(dBoosterData.rows()); DKMaille<double> dFundingAccrualBasis(dBoosterData.rows());
        DKMaille<double> dCouponDomestic(dBoosterData.rows()); DKMaille<double> dCouponForeign(dBoosterData.rows()); DKMaille<double> dIndexResetDates(dBoosterData.rows());
        // Notice dates managed at entry level
        for(unsigned uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dFXCouponResetDates.at(uiBD)=(dBoosterData.at(uiBD,1)-g_dSpotDate)/365.;
        for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dFXCouponPaymentDates.at(uiBD)=(dBoosterData.at(uiBD,2)-g_dSpotDate)/365.;
        for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dIndexStartDates.at(uiBD)=(dBoosterData.at(uiBD,3)-g_dSpotDate)/365.;
        for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dIndexEndDates.at(uiBD)=(dBoosterData.at(uiBD,4)-g_dSpotDate)/365.;
        for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dIndexPaymentDates.at(uiBD)=(dBoosterData.at(uiBD,5)-g_dSpotDate)/365.;
        for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dFundingSpread.at(uiBD)=(dBoosterData.at(uiBD,6));
        for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dFXNotionalMultiplier.at(uiBD)=(dBoosterData.at(uiBD,7));
        for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dFXStrike.at(uiBD)=(dBoosterData.at(uiBD,8));
        for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dFXCap.at(uiBD)=(dBoosterData.at(uiBD,9));
        for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dFXCouponAccrualBasis.at(uiBD)=(dBoosterData.at(uiBD,10));
        for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dFundingAccrualBasis.at(uiBD)=(dBoosterData.at(uiBD,11));
        for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dCouponDomestic.at(uiBD)=(dBoosterData.at(uiBD,12));
        for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dCouponForeign.at(uiBD)=(dBoosterData.at(uiBD,13));

        // Management of reset dates
        for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dIndexResetDates.at(uiBD)=(dCouponForeign.at(uiBD)-g_dSpotDate)/365.;

        // Xccy callable swap description
        if(g_dType==5.||g_dType==6.)
        {
            for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
                dIndexResetDates.at(uiBD)=dIndexStartDates.at(uiBD);
            for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
                dCouponDomestic.at(uiBD)=(dBoosterData.at(uiBD,12));
            for(uiBD=0;uiBD<dBoosterData.rows();uiBD++)
                dCouponForeign.at(uiBD)=(dBoosterData.at(uiBD,13));
        }



        DKMaille<bool> bFXCouponIsCalculated(dFXCouponPaymentDates.entries());
        DKMaille<bool> bLIBORCouponIsCalculated(dIndexPaymentDates.entries());
        // initialise
        for(unsigned int uiB=0;uiB<dFXCouponPaymentDates.entries();uiB++)
            bFXCouponIsCalculated.at(uiB)=false;
        for(uiB=0;uiB<dIndexPaymentDates.entries();uiB++)
            bLIBORCouponIsCalculated.at(uiB)=false;
        // Tests on the resulting Green's function
        double dLatticeDF =0.;
        double dLatticeDFForeign =0.;
        double dMarketDFForeign=0.;
        double dGreenOption=0.;

        DKMaille< DKMaille<double> > dOption;
        dOption.resize(g_uiNumSlices);
        DKMaille< DKMaille<unsigned int> > uiExerciseDecision;
        DKMaille< DKMaille<double> > dSurvivalProbability;
        DKMaille< DKMaille<double> > dADProbability;
        if(g_dSurvivalCalculation==1.0)
        {
            uiExerciseDecision.resize(g_uiNumSlices);
            dSurvivalProbability.resize(g_uiNumSlices);
            dADProbability.resize(g_uiNumSlices);
        }
        DKMaille< DKMaille<double> > dUnderlying;
        dUnderlying.resize(g_uiNumSlices);
        DKMaille< DKMaille<double> > dEuropeanOption;
        dEuropeanOption.resize(g_uiNumSlices);
        DKMaille< DKMaille<double> > dExisting;
        dExisting.resize(g_uiNumSlices);
        DKMaille< DKMaille<double> > dFunctionToSmooth;
        dFunctionToSmooth.resize(g_uiNumSlices);
        if(g_LATTICE_DATE.at(g_uiNumSlices-1)!=g_dNoticeDates.at(g_dNoticeDates.entries()-1))
            throw("Last notice is not last date in the lattice. Not acceptable lattice configuration");
        int uii=0;

        double dResult1=0.;
        double dResult2=0.;
        double dResult3=0.;
        double dResult4=0.;

        if(g_dType==4.4) // path-dependent product
        {
            DKMaille< DKMaille<double> > dTotalCoupon;
            dTotalCoupon.resize(g_uiNumSlices);
            DKMaille< DKMaille<double> > dActualCoupon;
            dActualCoupon.resize(g_uiNumSlices);
            DKMaille< DKMaille<double> > dTotalCouponBeforeReset;
            dTotalCouponBeforeReset.resize(g_uiNumSlices);
            // initialise total coupon
            dTotalCoupon.at(0).resize(1);
            dTotalCoupon.at(0).at(0)=0.;
            for(uii=0;uii<g_uiNumSlices;uii++)
            {
                unsigned int uiNextSlice=uii+1;
                // Preliminaries for each slice
                unsigned int uiiNumberOfStates=0;
                if(g_uiNumFactors==3) uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uii)*g_LATTICE_NUMBER_OF_STATES2.at(uii)*g_LATTICE_NUMBER_OF_STATES3.at(uii);
                if(g_uiNumFactors==2) uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uii)*g_LATTICE_NUMBER_OF_STATES2.at(uii);
                if(g_uiNumFactors==1) uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uii);
                dActualCoupon.at(uii).resize(uiiNumberOfStates);
                // zero the coupon holders
                for(unsigned int uf=0;uf<uiiNumberOfStates;uf++) dActualCoupon.at(uii).at(uf)=0.;
                int iNoticeDate=iIsNoticeDate(uii,g_dNoticeDates);
                // If we are on a notice date calculate Coupon
                if(iNoticeDate>-1)
                {
                    dTotalCouponBeforeReset.at(uii).resize(uiiNumberOfStates);
                    // Calculate the index of relevant cashflows that consitute the underlying at this date
                    unsigned int uiNoticeDateIndex=(unsigned int)iNoticeDate;
                    // unsigned int uiNFX=uiFirstFXCouponAfterNoticeDate(g_dNoticeDates.at(iNoticeDate),dFXCouponResetDates,dFXCouponPaymentDates);
                    // We have the right to set the first FX coupon by looking the funding leg since the funding leg has by definition the same frequency
                    unsigned int uiNFX=uiFirstIndexCouponAfterNoticeDate(g_dNoticeDates.at(iNoticeDate),dIndexResetDates);
                    unsigned int uiNSwap=uiFirstIndexCouponAfterNoticeDate(g_dNoticeDates.at(iNoticeDate),dIndexResetDates);
                    if(dIndexStartDates.at(uiNFX)-g_dNoticeDates.at(iNoticeDate)>2.1/365.)
                        throw("Reset date vs Start date configuration not supported");
                    if(g_dADPLimit!=1.0)
                    {
                        // 1. Coupon Calculation
                        unsigned int uiNumberOfBranches;
                        if(g_uiNumFactors==1) uiNumberOfBranches=3;
                        if(g_uiNumFactors==2) uiNumberOfBranches=9;
                        if(g_uiNumFactors==3) uiNumberOfBranches=27;
                        // Side 1 calcs
                        if(g_dType==4.4) // callable inverse FRN TARN structure
                        {
                            // no loop over coupon cashflows in forward mode
                            // for(unsigned int uiFXCashFlow=uiNFX;uiFXCashFlow<dFXCouponPaymentDates.entries();uiFXCashFlow++)
                            unsigned int uiFXCashFlow=uiNFX;
                            // no test for coupon calcs since they all have to be calculated
                            // if(bFXCouponIsCalculated.at(uiFXCashFlow)==false)
                            // Preparation for node calculation
                            static DKMaille<double> dNorm(3);
                            static DKMaille2D<double> dDecay(0,0);
                            if(g_uiProductModelCode>3) dDecay.resize(3,1);
                            if(g_uiProductModelCode<=3) dDecay.resize(3,g_uiNumFactors);
                            LatticeSliceGreenLIBORCalcs(uii,
                                                        dIndexStartDates.at(uiFXCashFlow),
                                                        dIndexEndDates.at(uiFXCashFlow),
                                                        dIndexPaymentDates.at(uiFXCashFlow),
                                                        dNorm,
                                                        dDecay);
                            // Initialise local volatility calculation
                            g_uiVolIsCalculated.at(uii)=0;
                            // Loop over nodes
                            for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                            {
                                if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                {
                                    dTotalCouponBeforeReset.at(uii).at(uk)=dTotalCoupon.at(uii).at(uk);
                                    // The function returns the value of the underlying.
                                    // In that case it is a simple FX Option cashflow
                                    dActualCoupon.at(uii).at(uk)-=LatticeInverseFRNAnalytics(uii,
                                                                  uk,
                                                                  dDecay,
                                                                  dIndexStartDates.at(uiFXCashFlow),
                                                                  dIndexEndDates.at(uiFXCashFlow),
                                                                  dIndexPaymentDates.at(uiFXCashFlow),
                                                                  dFXCouponAccrualBasis.at(uiFXCashFlow),
                                                                  dFundingSpread.at(uiFXCashFlow),
                                                                  dFXNotionalMultiplier.at(uiFXCashFlow),
                                                                  dFXStrike.at(uiFXCashFlow),
                                                                  dCouponDomestic.at(uiFXCashFlow),
                                                                  dNorm);
                                    dTotalCoupon.at(uii).at(uk)+=dActualCoupon.at(uii).at(uk);
                                } // if node exists
                            }	// loop on nodes
                            // bFXCouponIsCalculated.at(uiFXCashFlow)=true;
                        }
                        // End of Side 1 calcs
                        // Side 2 calcs
                        if(g_dType==4.4) // funding leg
                        {
                            unsigned int uiLIBORCashFlow=uiNSwap;
                            // Preparation for node calculation
                            static DKMaille<double> dNorm(3);
                            static DKMaille2D<double> dDecay(0,0);
                            if(g_uiProductModelCode>3) dDecay.resize(3,1);
                            if(g_uiProductModelCode<=3) dDecay.resize(3,g_uiNumFactors);
                            LatticeSliceGreenLIBORCalcs(uii,
                                                        dIndexStartDates.at(uiLIBORCashFlow),
                                                        dIndexEndDates.at(uiLIBORCashFlow),
                                                        dIndexPaymentDates.at(uiLIBORCashFlow),
                                                        dNorm,
                                                        dDecay);
                            // Loop over nodes
                            for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                            {
                                if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                {
                                    // The function returns the value of the underlying.
                                    // In that case it is a simple FX Option cashflow
                                    dActualCoupon.at(uii).at(uk)+=LatticeFloatSwapIndexAnalytics(uii,
                                                                  uk,
                                                                  dDecay,
                                                                  dIndexStartDates.at(uiLIBORCashFlow),
                                                                  dIndexEndDates.at(uiLIBORCashFlow),
                                                                  dIndexPaymentDates.at(uiLIBORCashFlow),
                                                                  dFundingAccrualBasis.at(uiLIBORCashFlow),
                                                                  dFundingSpread.at(uiLIBORCashFlow),
                                                                  dNorm);
                                } // if node exists
                            }	// loop on nodes
                            bLIBORCouponIsCalculated.at(uiLIBORCashFlow)=true;
                        }
                        unsigned int uiFXCashFlow=uiNFX;
                        double dPayOrNot=1.;
                        // 2. PV calculation
                        for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                        {
                            if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                            {
                                // Heading into the reset date the total coupon is greater than the redemption amount
                                if(dTotalCouponBeforeReset.at(uii).at(uk)<=-dFXCap.at(uiFXCashFlow))
                                    dPayOrNot=0.; // deal has already ended
                                else
                                    dPayOrNot=1.; // pay up
                                dResult1+=dActualCoupon.at(uii).at(uk)*(double)g_LATTICE_ARROW_DEBREU.at(uii).at(uk)*dPayOrNot;
                            } // if node exists
                        }
                        if(uii!=g_uiNumSlices-1)
                            PropagateTotalCoupon(uii,dTotalCoupon);
                    } // if Green's function exists
                    else
                        throw("You are missing the Green's function for such calculation");
                } // if on notice date
                else // if not a notice date just discount the option and the underlying
                {
                    if(uii!=g_uiNumSlices-1)
                        PropagateTotalCoupon(uii,dTotalCoupon);
                } // if on NOT notice date
                LatticeFreeGreen(uii);
                if(uii!=g_uiNumSlices-1)
                {
                    dActualCoupon.at(uii).clear();
                    dTotalCoupon.at(uii).clear();
                }
            } // loop on slices
        } // if g_dType==4.4 path dependent product
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if(g_dType!=4.4) // callable product
        {
            size_t noticeIdx=0;
            AnalyticFxOptionGlobalExtraction.Resize(g_dNoticeDates.entries(),3);
            for(size_t iii=0;iii<g_dNoticeDates.entries();++iii)
            {
                AnalyticFxOptionGlobalExtraction.Elt(iii,0)=0.0;
                AnalyticFxOptionGlobalExtraction.Elt(iii,1)=0.0;
                AnalyticFxOptionGlobalExtraction.Elt(iii,2)=0.0;
            }

            for(uii=g_uiNumSlices-1;uii>=0;uii--)
            {
                // Preliminaries for each slice
                unsigned int uiiNumberOfStates=0;
                if(g_uiNumFactors==3) uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uii)*g_LATTICE_NUMBER_OF_STATES2.at(uii)*g_LATTICE_NUMBER_OF_STATES3.at(uii);
                if(g_uiNumFactors==2) uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uii)*g_LATTICE_NUMBER_OF_STATES2.at(uii);
                if(g_uiNumFactors==1) uiiNumberOfStates=g_LATTICE_NUMBER_OF_STATES1.at(uii);
                dUnderlying.at(uii).resize(uiiNumberOfStates);
                dOption.at(uii).resize(uiiNumberOfStates);
                dEuropeanOption.at(uii).resize(uiiNumberOfStates);
                dExisting.at(uii).resize(uiiNumberOfStates);
                dFunctionToSmooth.at(uii).resize(uiiNumberOfStates);
                if(g_dSurvivalCalculation==1.0)
                {
                    uiExerciseDecision.at(uii).resize(uiiNumberOfStates);
                    dSurvivalProbability.at(uii).resize(uiiNumberOfStates);
                    dADProbability.at(uii).resize(uiiNumberOfStates);
                }
                // zero the product holders
                for(unsigned uf=0;uf<uiiNumberOfStates;uf++) dOption.at(uii).at(uf)=0.;
                for(uf=0;uf<uiiNumberOfStates;uf++) dUnderlying.at(uii).at(uf)=0.;
                for(uf=0;uf<uiiNumberOfStates;uf++) dEuropeanOption.at(uii).at(uf)=0.;
                for(uf=0;uf<uiiNumberOfStates;uf++) dExisting.at(uii).at(uf)=0.;
                for(uf=0;uf<uiiNumberOfStates;uf++) dFunctionToSmooth.at(uii).at(uf)=0.;
                if(g_dSurvivalCalculation==1.0)
                {
                    for(uf=0;uf<uiiNumberOfStates;uf++) uiExerciseDecision.at(uii).at(uf)=0;
                    for(uf=0;uf<uiiNumberOfStates;uf++) dSurvivalProbability.at(uii).at(uf)=0.;
                    for(uf=0;uf<uiiNumberOfStates;uf++) dADProbability.at(uii).at(uf)=0.;
                }
                int iNoticeDate=iIsNoticeDate(uii,g_dNoticeDates);
                // If we are on a notice date calculate both Option and Underlying
                if(iNoticeDate>-1)
                {
                    ++noticeIdx;
                    // Calculate the index of relevant cashflows that consitute the underlying at this date
                    unsigned int uiNoticeDateIndex=(unsigned int)iNoticeDate;
                    // unsigned int uiNFX=uiFirstFXCouponAfterNoticeDate(g_dNoticeDates.at(iNoticeDate),dFXCouponResetDates,dFXCouponPaymentDates);
                    // We have the right to set the first FX coupon by looking the funding leg since the funding leg has by definition the same frequency
                    unsigned int uiNFX=uiFirstIndexCouponAfterNoticeDate(g_dNoticeDates.at(iNoticeDate),dIndexResetDates);
                    unsigned int uiNSwap=uiFirstIndexCouponAfterNoticeDate(g_dNoticeDates.at(iNoticeDate),dIndexResetDates);
                    // if(g_dType==4.)
                    // {
                    // uiNFX=uiFirstCouponAfterNoticeDate(g_dNoticeDates.at(iNoticeDate),dFXCouponPaymentDates);
                    //  uiNSwap=uiFirstCouponAfterNoticeDate(g_dNoticeDates.at(iNoticeDate),dIndexPaymentDates);
                    // }

                    if(g_dADPLimit!=1.0)
                    {
                        // 1. Underlying Calculation
                        // First get discounted value if there is any
                        unsigned int uiNumberOfBranches;
                        if(g_uiNumFactors==1) uiNumberOfBranches=3;
                        if(g_uiNumFactors==2) uiNumberOfBranches=9;
                        if(g_uiNumFactors==3) uiNumberOfBranches=27;
                        if((g_dNoticeDates.entries()!=1)&&(g_LATTICE_DATE.at(uii)!=g_dNoticeDates.at(g_dNoticeDates.entries()-1))) // if it is not european and if not last notice
                        {
                            LatticeGetDiscountedSecurity(uii,dUnderlying);
                        } // if on european option

                        // Side 1 calcs
                        if(g_dType==3.) // callable PRCS
                        {
                            if(g_uiNumFactors==1||g_uiNumFactors==2) throw("Hybrid product cannot be priced with 1 or 2 factors");

                            // FX swap leg
                            // loop over FX Cashflows
                            for(unsigned int uiFXCashFlow=uiNFX;uiFXCashFlow<dFXCouponPaymentDates.entries();uiFXCashFlow++)
                            {
                                if(bFXCouponIsCalculated.at(uiFXCashFlow)==false)
                                {
                                    // Base case : coupon is an fx option
                                    // recognition pattern index is in domestic coupon
                                    if(dCouponDomestic.at(uiFXCashFlow)==-1.)
                                    {
                                        // Preparation for node calculation
                                        // This is general with regards to fx reset against fx payment date.
                                        double dNormBase1; 	double dNormForeign1; double dDecayDomestic1; double dDecayForeign1;
                                        double dNormBase2; 	double dNormForeign2; double dDecayDomestic2; double dDecayForeign2;
                                        LatticeSliceGreenFXCalcs(uii,
                                                                 dFXCouponResetDates.at(uiFXCashFlow),
                                                                 &dNormBase1,
                                                                 &dNormForeign1,
                                                                 &dDecayDomestic1,
                                                                 &dDecayForeign1);
                                        if(dFXCouponPaymentDates.at(uiFXCashFlow)!=dFXCouponResetDates.at(uiFXCashFlow))
                                        {
                                            LatticeSliceGreenFXCalcs(uii,
                                                                     dFXCouponPaymentDates.at(uiFXCashFlow),
                                                                     &dNormBase2,
                                                                     &dNormForeign2,
                                                                     &dDecayDomestic2,
                                                                     &dDecayForeign2);
                                        }
                                        else
                                        {
                                            dNormBase2=dNormBase1;
                                            dNormForeign2=dNormForeign1;
                                            dDecayDomestic2=dDecayDomestic1;
                                            dDecayForeign2=dDecayForeign1;
                                        }
                                        DKMaille<double> dNormBase(2); DKMaille<double> dNormForeign(2);
                                        DKMaille<double> dDecayDomestic(2); DKMaille<double> dDecayForeign(2);
                                        dNormBase.at(0)=dNormBase1; dNormBase.at(1)=dNormBase2;
                                        dNormForeign.at(0)=dNormForeign1; dNormForeign.at(1)=dNormForeign2;
                                        dDecayForeign.at(0)=dDecayForeign1; dDecayForeign.at(1)=dDecayForeign2;
                                        dDecayDomestic.at(0)=dDecayDomestic1; dDecayDomestic.at(1)=dDecayDomestic2;

/*********************  JMP : Extraction of all analytical prices of Fx options *************************/
                                        DKMaille<double> dNormBaseA(2); DKMaille<double> dNormForeignA(2);
                                        DKMaille<double> dDecayDomesticA(2); DKMaille<double> dDecayForeignA(2);
                                        LatticeSliceGreenFXCalcs(0,
                                                                 dFXCouponResetDates.at(uiFXCashFlow),
                                                                 &dNormBase1,
                                                                 &dNormForeign1,
                                                                 &dDecayDomestic1,
                                                                 &dDecayForeign1);
                                        if(dFXCouponPaymentDates.at(uiFXCashFlow)!=dFXCouponResetDates.at(uiFXCashFlow))
                                        {
                                            LatticeSliceGreenFXCalcs(0,
                                                                     dFXCouponPaymentDates.at(uiFXCashFlow),
                                                                     &dNormBase2,
                                                                     &dNormForeign2,
                                                                     &dDecayDomestic2,
                                                                     &dDecayForeign2);
                                        }
                                        else
                                        {
                                            dNormBase2=dNormBase1;
                                            dNormForeign2=dNormForeign1;
                                            dDecayDomestic2=dDecayDomestic1;
                                            dDecayForeign2=dDecayForeign1;
                                        }
                                        dNormBaseA.at(0)=dNormBase1; dNormBaseA.at(1)=dNormBase2;
                                        dNormForeignA.at(0)=dNormForeign1; dNormForeignA.at(1)=dNormForeign2;
                                        dDecayDomesticA.at(0)=dDecayDomestic1; dDecayDomesticA.at(1)=dDecayDomestic2;
                                        dDecayForeignA.at(0)=dDecayForeign1; dDecayForeignA.at(1)=dDecayForeign2;
                                        AnalyticFxOptionGlobalExtraction.Elt(g_dNoticeDates.entries()-noticeIdx,0)=dFXCouponResetDates.at(uiFXCashFlow)*365.0;
                                        AnalyticFxOptionGlobalExtraction.Elt(g_dNoticeDates.entries()-noticeIdx,1)+=10000
                                                                                *dFXNotionalMultiplier.at(uiFXCashFlow)
                                                                                *dFXCouponAccrualBasis.at(uiFXCashFlow)
                                                                                *LatticeFXOptionAnalytics(  0,0,
                                                                                                            dDecayDomesticA,
                                                                                                            dDecayForeignA,
                                                                                                            dFXCouponResetDates.at(uiFXCashFlow),
                                                                                                            dFXCouponPaymentDates.at(uiFXCashFlow),
                                                                                                            dFXStrike.at(uiFXCashFlow),
                                                                                                            1000000.0,
                                                                                                            dNormBaseA,
                                                                                                            dNormForeignA,
                                                                                                            0.0);
                                        AnalyticFxOptionGlobalExtraction.Elt(g_dNoticeDates.entries()-noticeIdx,2)-=10000
                                                                                *dFXNotionalMultiplier.at(uiFXCashFlow)
                                                                                *dFXCouponAccrualBasis.at(uiFXCashFlow)
                                                                                *LatticeFXOptionAnalytics(  0,0,
                                                                                                            dDecayDomesticA,
                                                                                                            dDecayForeignA,
                                                                                                            dFXCouponResetDates.at(uiFXCashFlow),
                                                                                                            dFXCouponPaymentDates.at(uiFXCashFlow),
                                                                                                            1000000.0,
                                                                                                            dFXCap.at(uiFXCashFlow),
                                                                                                            dNormBaseA,
                                                                                                            dNormForeignA,
                                                                                                            0.0);
/*********************  JMP : Extraction of all analytical prices of Fx options *************************/



                                        // Initialise local FX volatility calculation
                                        g_uiVolIsCalculated.at(uii)=0;


                                        // Loop over nodes
                                        for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                                        {
                                            if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                            {
                                                // The function returns the value of the underlying.
                                                // In that case it is a simple FX Option cashflow
/*****/
                                                dUnderlying.at(uii).at(uk)-=dFXNotionalMultiplier.at(uiFXCashFlow)
                                                                            *dFXCouponAccrualBasis.at(uiFXCashFlow)
                                                                            *LatticeFXOptionAnalytics(uii,
                                                                                                      uk,
                                                                                                      dDecayDomestic,
                                                                                                      dDecayForeign,
                                                                                                      dFXCouponResetDates.at(uiFXCashFlow),
                                                                                                      dFXCouponPaymentDates.at(uiFXCashFlow),
                                                                                                      dFXStrike.at(uiFXCashFlow),
                                                                                                      dFXCap.at(uiFXCashFlow),
                                                                                                      dNormBase,
                                                                                                      dNormForeign,
                                                                                                      0.0);
/*****/

/*********************  JMP TESTS : isolate only one fx option *************************
                                                if(FxOptIdxGlobal==g_dNoticeDates.entries() || g_dNoticeDates.entries()-noticeIdx == FxOptIdxGlobal)
                                                {   /// All fx options....                           .... or the selected one !
                                                    dUnderlying.at(uii).at(uk)-=dFXNotionalMultiplier.at(uiFXCashFlow)
                                                                            *dFXCouponAccrualBasis.at(uiFXCashFlow)
                                                                            *LatticeFXOptionAnalytics(uii,
                                                                                                      uk,
                                                                                                      dDecayDomestic,
                                                                                                      dDecayForeign,
                                                                                                      dFXCouponResetDates.at(uiFXCashFlow),
                                                                                                      dFXCouponPaymentDates.at(uiFXCashFlow),
                                                                                                      dFXStrike.at(uiFXCashFlow),
                                                                                                      dFXCap.at(uiFXCashFlow),
                                                                                                      dNormBase,
                                                                                                      dNormForeign,
                                                                                                      0.0);
                                                }
*********************  JMP TESTS : isolate only one fx option *************************/

                                            } // if node exists
                                        }	// loop on nodes
                                        bFXCouponIsCalculated.at(uiFXCashFlow)=true;
                                    } // if on coupon being an fx option

                                    // Rare case : coupon is a fixed coupon inside the call's underlying
                                    // recognition pattern index is in domestic coupon
                                    // fixed coupon is in strike
                                    // obligation to equal reset date to payment date
                                    if(dCouponDomestic.at(uiFXCashFlow)==+1.)
                                    {
                                        // Preparation for node calculation
                                        // if(dFXCouponResetDates.at(uiFXCashFlow)!=dFXCouponPaymentDates.at(uiFXCashFlow))
                                        // throw("Fixed Coupon Reset Date (irrelevant) is different from Payment Date");
                                        double dNormBase1; 	double dNormForeign1; double dDecayDomestic1; double dDecayForeign1;
                                        double dNormBase2; 	double dNormForeign2; double dDecayDomestic2; double dDecayForeign2;
                                        LatticeSliceGreenFXCalcs(uii,
                                                                 dFXCouponResetDates.at(uiFXCashFlow),
                                                                 &dNormBase1,
                                                                 &dNormForeign1,
                                                                 &dDecayDomestic1,
                                                                 &dDecayForeign1);
                                        if(dFXCouponPaymentDates.at(uiFXCashFlow)!=dFXCouponResetDates.at(uiFXCashFlow))
                                        {
                                            LatticeSliceGreenFXCalcs(uii,
                                                                     dFXCouponPaymentDates.at(uiFXCashFlow),
                                                                     &dNormBase2,
                                                                     &dNormForeign2,
                                                                     &dDecayDomestic2,
                                                                     &dDecayForeign2);
                                        }
                                        else
                                        {
                                            dNormBase2=dNormBase1;
                                            dNormForeign2=dNormForeign1;
                                            dDecayDomestic2=dDecayDomestic1;
                                            dDecayForeign2=dDecayForeign1;
                                        }
                                        DKMaille<double> dNormBase(2); DKMaille<double> dNormForeign(2);
                                        DKMaille<double> dDecayDomestic(2); DKMaille<double> dDecayForeign(2);
                                        dNormBase.at(0)=dNormBase1; dNormBase.at(1)=dNormBase2;
                                        dNormForeign.at(0)=dNormForeign1; dNormForeign.at(1)=dNormForeign2;
                                        dDecayForeign.at(0)=dDecayForeign1; dDecayForeign.at(1)=dDecayForeign2;
                                        dDecayDomestic.at(0)=dDecayDomestic1; dDecayDomestic.at(1)=dDecayDomestic2;
                                        // Loop over nodes
                                        for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                                        {
                                            if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                            {
                                                // The function returns the value of the underlying.
                                                // In that case it is a simple  cashflow
                                                dUnderlying.at(uii).at(uk)-=dFXNotionalMultiplier.at(uiFXCashFlow)
                                                                            *dFXCouponAccrualBasis.at(uiFXCashFlow)
                                                                            *LatticeCashFlowAnalytics(uii,
                                                                                                      uk,
                                                                                                      dDecayDomestic,
                                                                                                      dDecayForeign,
                                                                                                      dFXCouponResetDates.at(uiFXCashFlow),
                                                                                                      dNormBase,
                                                                                                      dNormForeign);
                                            } // if node exists
                                        }	// loop on nodes
                                        bFXCouponIsCalculated.at(uiFXCashFlow)=true;
                                    } // if on coupon being a fixed coupon
                                } // if on calculated cashflow
                            }	// loop on cashflows

                            // FX swap leg
                            // Calculation of simple redemption value only for last notice date
                            if(g_dRedemption==1.&&g_dRedemptionOption==0.&&g_LATTICE_DATE.at(uii)==g_dNoticeDates.at(g_dNoticeDates.entries()-1))
                            {
                                // Preparation for node calculation
                                // This is general with regards to fx reset against fx payment date.
                                // fx reset date
                                double dRedemptionDate=(g_dRedemptionDate-g_dSpotDate)/365.;
                                // double dRedemptionNoticeDate=(g_dRedemptionNoticeDate-g_dSpotDate)/365.;
                                double dRedemptionNoticeDate=(g_dRedemptionDate-g_dSpotDate)/365.;
                                double dResetDate=dRedemptionNoticeDate;
                                double dPaymentDate=dRedemptionDate;
                                double dNormBase1; 	double dNormForeign1; double dDecayDomestic1; double dDecayForeign1;
                                double dNormBase2; 	double dNormForeign2; double dDecayDomestic2; double dDecayForeign2;
                                LatticeSliceGreenFXCalcs(uii,
                                                         dResetDate,
                                                         &dNormBase1,
                                                         &dNormForeign1,
                                                         &dDecayDomestic1,
                                                         &dDecayForeign1);
                                if(dResetDate!=dPaymentDate)
                                {
                                    LatticeSliceGreenFXCalcs(uii,
                                                             dPaymentDate,
                                                             &dNormBase2,
                                                             &dNormForeign2,
                                                             &dDecayDomestic2,
                                                             &dDecayForeign2);
                                }
                                else
                                {
                                    dNormBase2=dNormBase1;
                                    dNormForeign2=dNormForeign1;
                                    dDecayDomestic2=dDecayDomestic1;
                                    dDecayForeign2=dDecayForeign1;
                                }
                                DKMaille<double> dNormBase(2); DKMaille<double> dNormForeign(2);
                                DKMaille<double> dDecayDomestic(2); DKMaille<double> dDecayForeign(2);
                                dNormBase.at(0)=dNormBase1; dNormBase.at(1)=dNormBase2;
                                dNormForeign.at(0)=dNormForeign1; dNormForeign.at(1)=dNormForeign2;
                                dDecayForeign.at(0)=dDecayForeign1; dDecayForeign.at(1)=dDecayForeign2;
                                dDecayDomestic.at(0)=dDecayDomestic1; dDecayDomestic.at(1)=dDecayDomestic2;

                                // Initialise local FX volatility calculation
                                g_uiVolIsCalculated.at(uii)=0;
                                // Loop over nodes
                                for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                                {
                                    if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                    {
                                        // The function returns the value of the underlying.
                                        // In that case it is a simple FX Option cashflow
                                        dUnderlying.at(uii).at(uk)-=(1./g_dRedemptionStrike)
                                                                    *LatticeFXForwardAnalytics(uii,
                                                                                               uk,
                                                                                               dDecayDomestic,
                                                                                               dDecayForeign,
                                                                                               dResetDate,
                                                                                               dPaymentDate,
                                                                                               g_dRedemptionStrike,
                                                                                               dNormBase,
                                                                                               dNormForeign,
                                                                                               0.0);
                                    } // if node exists
                                }	// loop on nodes
                            } // if on Redemption

                            // FX swap leg
                            // Calculation of option redemption value only for last notice date
                            if(g_dRedemption==1.&&g_dRedemptionOption==1.&&g_LATTICE_DATE.at(uii)==g_dNoticeDates.at(g_dNoticeDates.entries()-1))
                            {
                                // Preparation for node calculation
                                // This is general with regards to fx reset against fx payment date.
                                // fx reset date
                                double dRedemptionDate=(g_dRedemptionDate-g_dSpotDate)/365.;
                                double dRedemptionNoticeDate=(g_dRedemptionNoticeDate-g_dSpotDate)/365.;
                                double dResetDate=dRedemptionNoticeDate;
                                double dPaymentDate=dRedemptionDate;
                                double dNormBase1; 	double dNormForeign1; double dDecayDomestic1; double dDecayForeign1;
                                double dNormBase2; 	double dNormForeign2; double dDecayDomestic2; double dDecayForeign2;
                                LatticeSliceGreenFXCalcs(uii,
                                                         dResetDate,
                                                         &dNormBase1,
                                                         &dNormForeign1,
                                                         &dDecayDomestic1,
                                                         &dDecayForeign1);
                                if(dResetDate!=dPaymentDate)
                                {
                                    LatticeSliceGreenFXCalcs(uii,
                                                             dPaymentDate,
                                                             &dNormBase2,
                                                             &dNormForeign2,
                                                             &dDecayDomestic2,
                                                             &dDecayForeign2);
                                }
                                else
                                {
                                    dNormBase2=dNormBase1;
                                    dNormForeign2=dNormForeign1;
                                    dDecayDomestic2=dDecayDomestic1;
                                    dDecayForeign2=dDecayForeign1;
                                }
                                DKMaille<double> dNormBase(2); DKMaille<double> dNormForeign(2);
                                DKMaille<double> dDecayDomestic(2); DKMaille<double> dDecayForeign(2);
                                dNormBase.at(0)=dNormBase1; dNormBase.at(1)=dNormBase2;
                                dNormForeign.at(0)=dNormForeign1; dNormForeign.at(1)=dNormForeign2;
                                dDecayForeign.at(0)=dDecayForeign1; dDecayForeign.at(1)=dDecayForeign2;
                                dDecayDomestic.at(0)=dDecayDomestic1; dDecayDomestic.at(1)=dDecayDomestic2;
                                // Initialise local FX volatility calculation
                                g_uiVolIsCalculated.at(uii)=0;
                                // Loop over nodes
                                for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                                {
                                    if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                    {
                                        // The function returns the value of the underlying.
                                        // In that case it is a simple FX Option cashflow
                                        dUnderlying.at(uii).at(uk)+=(1./g_dRedemptionStrike)
                                                                    *LatticeFXOptionAnalytics(uii,
                                                                                              uk,
                                                                                              dDecayDomestic,
                                                                                              dDecayForeign,
                                                                                              dResetDate,
                                                                                              dPaymentDate,
                                                                                              g_dRedemptionStrike,
                                                                                              100000.,
                                                                                              dNormBase,
                                                                                              dNormForeign,
                                                                                              1.0);
                                    } // if node exists
                                }	// loop on nodes
                            } // if on RedemptionOption
                        }
                        else if(g_dType==4.) // callable swap
                        {
                            if(g_uiProductModelCode!=g_uiNumFactors) throw("Hybrid classes derive from single currency fixed income product and not the other way round");
                            // Fixed leg of the swap
                            // loop over swap cashflows
                            for(unsigned int uiFXCashFlow=uiNFX;uiFXCashFlow<dFXCouponPaymentDates.entries();uiFXCashFlow++)
                            {
                                if(bFXCouponIsCalculated.at(uiFXCashFlow)==false)
                                {
                                    // Preparation for node calculation
                                    double dNorm;
                                    DKMaille<double> dDecay(g_uiNumFactors);
                                    LatticeSliceGreenCashCalcs(uii,
                                                               dFXCouponPaymentDates.at(uiFXCashFlow),
                                                               &dNorm,
                                                               dDecay);

                                    // Loop over nodes
                                    for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                                    {
                                        if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                        {
                                            // The function returns the value of the underlying.
                                            // In that case it is a simple FX Option cashflow
                                            dUnderlying.at(uii).at(uk)-=LatticeFixedSwapIndexAnalytics(uii,
                                                                        uk,
                                                                        dNorm,
                                                                        dDecay,
                                                                        dFXCouponPaymentDates.at(uiFXCashFlow),
                                                                        dFXCouponAccrualBasis.at(uiFXCashFlow),
                                                                        dCouponDomestic.at(uiFXCashFlow));
                                        } // if node exists
                                    }	// loop on nodes
                                    bFXCouponIsCalculated.at(uiFXCashFlow)=true;
                                } // if on calculated coupon
                            }	// loop on cashflows
                        }
                        else if(g_dType==5.||g_dType==6.) // callable cross currency swap
                        {
                            if(g_uiNumFactors==1||g_uiNumFactors==2) throw("Hybrid product cannot be priced with 1 or 2 factors");
                            // domestic swap leg
                            // loop over domestic cashflows
                            for(unsigned int uiFXCashFlow=uiNFX;uiFXCashFlow<dFXCouponPaymentDates.entries();uiFXCashFlow++)
                            {
                                if(bFXCouponIsCalculated.at(uiFXCashFlow)==false)
                                {
                                    // Preparation for node calculation
                                    // This is with regards to payment date.
                                    double dNormBase1; 	double dNormForeign1; double dDecayDomestic1; double dDecayForeign1;
                                    LatticeSliceGreenFXCalcs(uii,
                                                             dFXCouponPaymentDates.at(uiFXCashFlow),
                                                             &dNormBase1,
                                                             &dNormForeign1,
                                                             &dDecayDomestic1,
                                                             &dDecayForeign1);
                                    // Loop over nodes
                                    double dSide=1.;
                                    if(g_dType==6.0) dSide=-1.;
                                    for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                                    {
                                        if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                        {
                                            // The function returns the value of the underlying.
                                            // In that case it is a simple FX Option cashflow
                                            dUnderlying.at(uii).at(uk)-=dSide
                                                                        *LatticeFixedSwapIndexAnalyticsDomestic(uii,
                                                                                                                uk,
                                                                                                                dNormBase1,
                                                                                                                dDecayDomestic1,
                                                                                                                dFXCouponPaymentDates.at(uiFXCashFlow),
                                                                                                                dFXCouponAccrualBasis.at(uiFXCashFlow),
                                                                                                                dCouponDomestic.at(uiFXCashFlow));
                                        } // if node exists
                                    }	// loop on nodes
                                    bFXCouponIsCalculated.at(uiFXCashFlow)=true;
                                } // if on calculated cashflow
                            }	// loop on cashflows
                        }
                        else if(g_dType==4.33) // callable inverse FRN // initially on LIBOR to extend to 1Y and 2Y
                        {
                            // loop over coupon cashflows
                            for(unsigned int uiFXCashFlow=uiNFX;uiFXCashFlow<dFXCouponPaymentDates.entries();uiFXCashFlow++)
                            {
                                if(bFXCouponIsCalculated.at(uiFXCashFlow)==false)
                                {
                                    // Preparation for node calculation
                                    static DKMaille<double> dNorm(3);
                                    static DKMaille2D<double> dDecay(0,0);
                                    if(g_uiProductModelCode>3) dDecay.resize(3,1);
                                    if(g_uiProductModelCode<=3) dDecay.resize(3,g_uiNumFactors);
                                    LatticeSliceGreenLIBORCalcs(uii,
                                                                dIndexStartDates.at(uiFXCashFlow),
                                                                dIndexEndDates.at(uiFXCashFlow),
                                                                dIndexPaymentDates.at(uiFXCashFlow),
                                                                dNorm,
                                                                dDecay);

                                    // Initialise local volatility calculation
                                    g_uiVolIsCalculated.at(uii)=0;
                                    // Loop over nodes
                                    for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                                    {
                                        if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                        {
                                            // The function returns the value of the underlying.
                                            // In that case it is a simple FX Option cashflow
                                            dUnderlying.at(uii).at(uk)-=LatticeInverseFRNAnalytics(uii,
                                                                        uk,
                                                                        dDecay,
                                                                        dIndexStartDates.at(uiFXCashFlow),
                                                                        dIndexEndDates.at(uiFXCashFlow),
                                                                        dIndexPaymentDates.at(uiFXCashFlow),
                                                                        dFXCouponAccrualBasis.at(uiFXCashFlow),
                                                                        dFundingSpread.at(uiFXCashFlow),
                                                                        dFXNotionalMultiplier.at(uiFXCashFlow),
                                                                        dFXStrike.at(uiFXCashFlow),
                                                                        dCouponDomestic.at(uiFXCashFlow),
                                                                        dNorm);
                                        } // if node exists
                                    }	// loop on nodes
                                    bFXCouponIsCalculated.at(uiFXCashFlow)=true;
                                } // if on calculated cashflows
                            } // loop on cashflows
                        }
                        // End of Side 1 calcs

                        // Side 2 calcs
                        if(g_dType==3.||g_dType==4.||g_dType==4.33||g_dType==4.66) // either callable PRCS or callable single ccy swap therefore funded by floating leg
                        {
                            // loop over float cashflows
                            for(unsigned int uiLIBORCashFlow=uiNSwap;uiLIBORCashFlow<dIndexPaymentDates.entries();uiLIBORCashFlow++)
                            {
                                if(bLIBORCouponIsCalculated.at(uiLIBORCashFlow)==false)
                                {
                                    // Preparation for node calculation
                                    static DKMaille<double> dNorm(3);
                                    static DKMaille2D<double> dDecay(0,0);
                                    if(g_uiProductModelCode>3) dDecay.resize(3,1);
                                    if(g_uiProductModelCode<=3) dDecay.resize(3,g_uiNumFactors);
                                    LatticeSliceGreenLIBORCalcs(uii,
                                                                dIndexStartDates.at(uiLIBORCashFlow),
                                                                dIndexEndDates.at(uiLIBORCashFlow),
                                                                dIndexPaymentDates.at(uiLIBORCashFlow),
                                                                dNorm,
                                                                dDecay);

                                    // Loop over nodes
                                    for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                                    {
                                        if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                        {
                                            // The function returns the value of the underlying.
                                            // In that case it is a simple FX Option cashflow
                                            dUnderlying.at(uii).at(uk)+=LatticeFloatSwapIndexAnalytics(uii,
                                                                        uk,
                                                                        dDecay,
                                                                        dIndexStartDates.at(uiLIBORCashFlow),
                                                                        dIndexEndDates.at(uiLIBORCashFlow),
                                                                        dIndexPaymentDates.at(uiLIBORCashFlow),
                                                                        dFundingAccrualBasis.at(uiLIBORCashFlow),
                                                                        dFundingSpread.at(uiLIBORCashFlow),
                                                                        dNorm);
                                        } // if node exists
                                    }	// loop on nodes
                                    bLIBORCouponIsCalculated.at(uiLIBORCashFlow)=true;
                                } // if on calculated cashflows
                            } // loop on cashflows
                        }

                        if(g_dType==5.||g_dType==6.) // callable cross currency swap
                        {
                            // loop over foreign cashflows
                            for(unsigned int uiLIBORCashFlow=uiNSwap;uiLIBORCashFlow<dIndexPaymentDates.entries();uiLIBORCashFlow++)
                            {
                                if(bLIBORCouponIsCalculated.at(uiLIBORCashFlow)==false)
                                {
                                    // Preparation for node calculation
                                    // This is with regards to payment date.
                                    // Foreign leg as opposed to LIBOR leg that is why the names of the variables refer to LIBOR
                                    // tidy this up later on
                                    double dNormBase1;
                                    double dNormForeign1;
                                    double dDecayDomestic1;
                                    double dDecayForeign1;
                                    LatticeSliceGreenFXCalcs(uii,
                                                             dIndexPaymentDates.at(uiLIBORCashFlow),
                                                             &dNormBase1,
                                                             &dNormForeign1,
                                                             &dDecayDomestic1,
                                                             &dDecayForeign1);
                                    double dSide=1.;
                                    if(g_dType==6.0) dSide=-1.;
                                    // Loop over nodes
                                    for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                                    {
                                        if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                        {
                                            // The function returns the value of the underlying.
                                            // In that case it is a simple FX Option cashflow
                                            dUnderlying.at(uii).at(uk)+=dSide
                                                                        *dFXNotionalMultiplier.at(uiLIBORCashFlow)
                                                                        *LatticeFixedSwapIndexAnalyticsForeign(uii,
                                                                                                               uk,
                                                                                                               dNormBase1,
                                                                                                               dDecayDomestic1,
                                                                                                               dNormForeign1,
                                                                                                               dDecayForeign1,
                                                                                                               dIndexPaymentDates.at(uiLIBORCashFlow),
                                                                                                               dFundingAccrualBasis.at(uiLIBORCashFlow),
                                                                                                               dCouponForeign.at(uiLIBORCashFlow));
                                        } // if node exists
                                    }	// loop on nodes
                                    bLIBORCouponIsCalculated.at(uiLIBORCashFlow)=true;
                                } // if on calculated cashflows
                            } // loop on cashflows
                        }


                        // 2. Option calculation
                        // Case of last notice
                        if(g_LATTICE_DATE.at(uii)==g_dNoticeDates.at(g_dNoticeDates.entries()-1))
                        {
                            // only one notice case then fill the European Holder
                            if(g_dNoticeDates.entries()==1)
                            {
                                for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                                {
                                    if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                    {
                                        if(dUnderlying.at(uii).at(uk)<=0) dEuropeanOption.at(uii).at(uk)=-dUnderlying.at(uii).at(uk);
                                        else dEuropeanOption.at(uii).at(uk)=0.;
                                        // Record of exercise decision
                                        if(g_dSurvivalCalculation==1.0)
                                        {
                                            if(dUnderlying.at(uii).at(uk)<=0) uiExerciseDecision.at(uii).at(uk)=1;
                                            else uiExerciseDecision.at(uii).at(uk)=0;
                                        }
                                    } // if node exists
                                }
                            }
                            for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                            {
                                if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                {
                                    if(dUnderlying.at(uii).at(uk)<=0) dOption.at(uii).at(uk)=-dUnderlying.at(uii).at(uk);
                                    else dOption.at(uii).at(uk)=0.;
                                    dFunctionToSmooth.at(uii).at(uk)=-dUnderlying.at(uii).at(uk);
                                    // Record of exercise decision
                                    if(g_dSurvivalCalculation==1.0)
                                    {
                                        if(dUnderlying.at(uii).at(uk)<=0) uiExerciseDecision.at(uii).at(uk)=1;
                                        else uiExerciseDecision.at(uii).at(uk)=0;
                                    }
                                } // if node exists
                            }
                            if(g_uiSmoothing==1) LatticeApplySmoothing(uii,1,dOption,dFunctionToSmooth);
                        } // if on last notice
                        else // previous notices
                        {
                            // fill the European Holder
                            if(g_LATTICE_DATE.at(uii)==g_dNoticeDates.at(0))
                            {
                                for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                                {
                                    if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                    {
                                        if(dUnderlying.at(uii).at(uk)<=0) dEuropeanOption.at(uii).at(uk)=-dUnderlying.at(uii).at(uk);
                                        else dEuropeanOption.at(uii).at(uk)=0.;
                                    } // if node exists
                                }
                            }
                            unsigned int uiNumberOfBranches;
                            if(g_uiNumFactors==1) uiNumberOfBranches=3;
                            if(g_uiNumFactors==2) uiNumberOfBranches=9;
                            if(g_uiNumFactors==3) uiNumberOfBranches=27;
                            if(g_dNoticeDates.entries()!=1) // if it is not european
                            {
                                // Get discounted option
                                LatticeGetDiscountedSecurity(uii,dOption);
                                for(unsigned uk=0;uk<uiiNumberOfStates;uk++)
                                {
                                    if(g_LATTICE_NODE_EXISTENCE.at(uii).at(uk)==true)
                                    {
                                        dExisting.at(uii).at(uk)=dOption.at(uii).at(uk);
                                        dFunctionToSmooth.at(uii).at(uk)=-dUnderlying.at(uii).at(uk)-dExisting.at(uii).at(uk);
                                        double dHold=0.;
                                        if(dUnderlying.at(uii).at(uk)<=0) dHold=-dUnderlying.at(uii).at(uk);
                                        // Record of exercise decision
                                        if(g_dSurvivalCalculation==1.0)
                                        {
                                            if(dUnderlying.at(uii).at(uk)<=0)
                                            {
                                                if(dOption.at(uii).at(uk)<=dHold) uiExerciseDecision.at(uii).at(uk)=1;
                                                else uiExerciseDecision.at(uii).at(uk)=0;
                                            }
                                            else uiExerciseDecision.at(uii).at(uk)=0;
                                        }
                                        // if optimal replace option by early exercise price
                                        if(dOption.at(uii).at(uk)<=dHold) dOption.at(uii).at(uk)=dHold;
                                    } // if node exists
                                } // loop on nodes
                                if(g_uiSmoothing==1) LatticeApplySmoothing(uii,1,dOption,dFunctionToSmooth);
                            } // if on european option
                            // Do nothing
                        }
                    } // if Green's function exists
                    else
                        throw("You are missing the Green's function for such calculation");
                } // if on notice date
                else // if not a notice date just discount the option and the underlying
                {
                    unsigned int uiNumberOfBranches;
                    if(g_uiNumFactors==1) uiNumberOfBranches=3;
                    if(g_uiNumFactors==2) uiNumberOfBranches=9;
                    if(g_uiNumFactors==3) uiNumberOfBranches=27;

                    LatticeGetDiscountedSecurity(uii,dOption);
                    LatticeGetDiscountedSecurity(uii,dUnderlying);
                    // if before first notice get the european option and discount it
                    if(g_LATTICE_DATE.at(uii)<g_dNoticeDates.at(0))
                    {
                        LatticeGetDiscountedSecurity(uii,dEuropeanOption);
                    } // if on before first notice
                } // if on NOT notice date

                if(g_dSurvivalCalculation==0.0)
                {
                    if(g_dDeltaFlag==1.&&g_uiDeltaTick==1) LatticeFreeGreen(uii);
                    if(g_dDeltaFlag==0.&&g_uiDeltaTick==-1) LatticeFreeGreen(uii);
                }
                if(uii!=g_uiNumSlices-1)
                {
                    dOption.at(uii+1).clear();
                    dUnderlying.at(uii+1).clear();
                    dExisting.at(uii+1).clear();
                    dFunctionToSmooth.at(uii+1).clear();
                    dEuropeanOption.at(uii+1).clear();
                }
            } // loop on slices
        } // if g_dType!=4.4 callable and non-path dependent product

        static DKMaille<double> ret;
        if(g_dSurvivalCalculation==0.0)
        {
            ret.resize(3);
            if(g_dType==4.4)
            {
                ret.at(0)=10000.*dResult1;
                ret.at(1)=10000.*0.0;
                ret.at(2)=10000.*0.0;
            }
            else
            {
                ret.at(0)=10000.*dUnderlying.at(0).at(0);
                ret.at(1)=10000.*dOption.at(0).at(0);
                ret.at(2)=10000.*dEuropeanOption.at(0).at(0);
            }
        }

        DKMaille<double> dOutput(g_uiNumSlices);
        if(g_dSurvivalCalculation==1.0)
        {
            ret.resize(3+2.*g_uiNumSlices);
            ret.at(0)=10000.*dUnderlying.at(0).at(0);
            ret.at(1)=10000.*dOption.at(0).at(0);
            ret.at(2)=10000.*dEuropeanOption.at(0).at(0);

            dOutput.at(0)=1.;
            dSurvivalProbability.at(0).at(0)=1.;
            dADProbability.at(0).at(0)=1.;

            for(unsigned int uiSlice=0;uiSlice<g_uiNumSlices-1;uiSlice++)
            {
                unsigned int uiNextSlice=uiSlice+1;
                double dt = g_LATTICE_PERIOD_TIME_SPAN.at(uiSlice);
                DKMaille<long> currentNodeStates(g_uiNumFactors);
                DKMaille<long> targetNodeStates(g_uiNumFactors);
                DKMaille<long> iBranch(g_uiNumFactors);
                double dProbability=0.;
                double dContribution=0.;
                double dExerciseFactor=0.;

                dOutput.at(uiNextSlice)=0.;
                switch(g_uiNumFactors)
                {
                case 1:
                    {
                        for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
                        {
                            currentNodeStates.at(0)=i1;
                            unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                            // Get the states from the next slice that we will need to use in connecting
                            if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                            {
                                double dDiscountConnection=LatticeDiscount(uiSlice,uiIndex);
                                for(unsigned int um=0;um<3;um++)
                                {
                                    dProbability=(1./(double)g_dConvert)*(double)g_LATTICE_GREEN.at(uiSlice,um).at(uiIndex)/dDiscountConnection;
                                    if(uiExerciseDecision.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))==1) dExerciseFactor=0.;
                                    else dExerciseFactor=1.;
                                    dContribution=dSurvivalProbability.at(uiSlice,um).at(uiIndex)*dProbability*dExerciseFactor;
                                    dSurvivalProbability.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))+=dContribution;
                                    dADProbability.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))+=dADProbability.at(uiSlice,um).at(uiIndex)
                                            *dProbability;
                                    dOutput.at(uiNextSlice)+=dContribution;
                                }
                            }
                        } // while loop over nodes
                    }
                    break;
                case 2:
                    {
                        for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
                        {
                            for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                            {
                                currentNodeStates.at(0)=i1;
                                currentNodeStates.at(1)=i2;
                                unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                                // Get the states from the next slice that we will need to use in connecting
                                if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                                {
                                    double dDiscountConnection=LatticeDiscount(uiSlice,uiIndex);
                                    for(unsigned int um=0;um<9;um++)
                                    {
                                        dProbability=(1./(double)g_dConvert)*(double)g_LATTICE_GREEN.at(uiSlice,um).at(uiIndex)/dDiscountConnection;
                                        if(uiExerciseDecision.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))==1) dExerciseFactor=0.;
                                        else dExerciseFactor=1.;
                                        dContribution=dSurvivalProbability.at(uiSlice,um).at(uiIndex)*dProbability*dExerciseFactor;
                                        dSurvivalProbability.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))+=dContribution;
                                        dADProbability.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))+=dADProbability.at(uiSlice,um).at(uiIndex)
                                                *dProbability;
                                        dOutput.at(uiNextSlice)+=dContribution;
                                    }
                                }
                            }
                        } // while loop over nodes
                    }
                    break;
                case 3:
                    {
                        for(int i1=g_LATTICE_MIN_STATE1.at(uiSlice);i1<=g_LATTICE_MAX_STATE1.at(uiSlice);i1++)
                        {
                            for(int i2=g_LATTICE_MIN_STATE2.at(uiSlice);i2<=g_LATTICE_MAX_STATE2.at(uiSlice);i2++)
                            {
                                for(int i3=g_LATTICE_MIN_STATE3.at(uiSlice);i3<=g_LATTICE_MAX_STATE3.at(uiSlice);i3++)
                                {
                                    currentNodeStates.at(0)=i1;
                                    currentNodeStates.at(1)=i2;
                                    currentNodeStates.at(2)=i3;
                                    unsigned int uiIndex=LatticeIndexConvert(uiSlice,currentNodeStates);
                                    // Get the states from the next slice that we will need to use in connecting
                                    if(g_LATTICE_NODE_EXISTENCE.at(uiSlice).at(uiIndex))
                                    {
                                        double dDiscountConnection=LatticeDiscount(uiSlice,uiIndex);
                                        for(unsigned int um=0;um<27;um++)
                                        {
                                            dProbability=(1./(double)g_dConvert)*(double)g_LATTICE_GREEN.at(uiSlice,um).at(uiIndex)/dDiscountConnection;
                                            if(uiExerciseDecision.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))==1) dExerciseFactor=0.;
                                            else dExerciseFactor=1.;
                                            dContribution=dSurvivalProbability.at(uiSlice,um).at(uiIndex)*dProbability*dExerciseFactor;
                                            dSurvivalProbability.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))+=dContribution;
                                            dADProbability.at(uiNextSlice).at(g_LATTICE_CONNECTION.at(uiSlice,um).at(uiIndex))+=dADProbability.at(uiSlice,um).at(uiIndex)
                                                    *dProbability;
                                            dOutput.at(uiNextSlice)+=dContribution;
                                        }
                                    }
                                }
                            }
                        } // while loop over nodes
                    }
                    break;
                }

                if(g_dDeltaFlag==1.&&g_uiDeltaTick==1) LatticeFreeGreen(uiSlice);
                if(g_dDeltaFlag==0.&&g_uiDeltaTick==-1) LatticeFreeGreen(uiSlice);

                dSurvivalProbability.at(uiSlice).clear();
                dADProbability.at(uiSlice).clear();
                uiExerciseDecision.at(uiSlice).clear();

            } // loop over slices

            for(unsigned int uiM=0;uiM<g_uiNumSlices;uiM++)
            {
                ret.at(3+uiM)=g_LATTICE_DATE.at(uiM);
            }
            for(uiM=0;uiM<g_uiNumSlices;uiM++)
            {
                ret.at(3+uiM+g_uiNumSlices)=dOutput.at(uiM);
            }
        } // if on calculation of survival probability
        return ret;
    } // if on g_dType >= 3
    else
    {
        static DKMaille<double> ret(3);
        ret.at(0)=1.;
        ret.at(1)=1.;
        ret.at(2)=1.;
        return ret;
    }
}



void SelectBootstrappingDates(DKMaille<double> &dNewNoticeDates,DKMaille<double> &dNoticeDates,DKMaille<double> &dSurfaceX)
{
    DKMaille<double> dRet;
    unsigned int uiRows=dSurfaceX.entries();
    DKMaille<double>  OptionTenors;
    if(dSurfaceX.at(0)!=0.)
    {
        OptionTenors.resize(uiRows+1);
        OptionTenors.at(0)  = 0.0;
        for (unsigned int uiY = 0;  uiY <  uiRows; uiY++ )
        {
            OptionTenors.at(uiY+1)  = dSurfaceX[uiY];
        }
    }
    else
    {
        OptionTenors.resize(uiRows);
        for (unsigned int uiY = 0;  uiY <  uiRows; uiY++ )
        {
            OptionTenors.at(uiY)  = dSurfaceX[uiY];
        }
    }

    double  dModelDate    = 0.;
    double  dLowerModelDateBound  = 0.;
    double  dUpperModelDateBound  = 0.;

    if(dNoticeDates.entries()==0)
        throw("No notice dates. Unhandled case");

    unsigned int uk=0;
    dRet.insert(dNoticeDates.at(0));
    uk++;
    unsigned int uiPrevious=0;
    // if not European
    if(dNoticeDates.entries()>2)
    {
        //Loop over the option expiry inputs (including the zero model time)
        unsigned  uiLowWindowIndex  = 0;
        unsigned  uiLowNoticeIndex  = uk;

        for (unsigned int ui  = uiLowWindowIndex;  ui  <= uiRows; ui++  )
        {

            dLowerModelDateBound  = OptionTenors.at(ui);
            if(ui<uiRows) dUpperModelDateBound  = OptionTenors.at(ui+1);
            else dUpperModelDateBound=100.;


            for ( uk  = uiLowNoticeIndex; uk  < dNoticeDates.entries()-1;  uk++  )
            {

                dModelDate    = dNoticeDates.at(uk);

                if  ( ( dLowerModelDateBound  <= dModelDate  ) &&  ( dModelDate  < dUpperModelDateBound  ) )
                {
                    if  ( (uk  > 0) && (dRet.at(uiPrevious) < dLowerModelDateBound) )
                    {
                        uiPrevious++;
                        dRet.insert(dNoticeDates.at(uk));//Insert the first date in any given window
                        uiLowNoticeIndex++;
                    }
                    break;
                }
            }
        }
    }
    dNewNoticeDates=dRet;
}

double Correction(DKMaille<double> &dVolStrip,DKMaille<double> &dVolStripDates)
{

    // Distance from middle between strip dates
    double dX1=0.5*(dVolStripDates.at(1)-dVolStripDates.at(0))/365.;
    if(dVolStripDates.entries()<=3) throw("You should not be here");
    // Calculate first derivative
    double dX10=dX1;
    double dX20=0.5*(dVolStripDates.at(2)-dVolStripDates.at(1))/365.;
    double dXX1=dX10;
    double dXX2=(dVolStripDates.at(1)-dVolStripDates.at(0))/365.+dX20;
    double dYY1=dVolStrip.at(0);
    double dYY2=dVolStrip.at(1);
    double dDerivative=(dYY1-dYY2)/(dXX1-dXX2);
    return dVolStrip.at(0)-dDerivative*dX1;

}




void LatticeFillBootstrappedVolsSingleCurrencyFixedIncome(DKMaille<double> dBaseDates,
        DKMaille<double> dBaseRates,
        DKMaille<double> dBaseRatesNoBasis,
        DKMaille<double> dStdDevBaseX,
        DKMaille<double> dStdDevBaseY,
        DKMaille2D<double> dStdDevBaseZ,
        double dSpotDate,
        DKMaille2D<double> dBoosterData,
        double dMeanReversionBase, // MeanReversion1
        double dMeanReversionForeign, // MeanReversion2
        double dSpotFXVol, // // MeanReversion3
        double dStdDevForeign, // RelativeFactor12
        double dSpotFX, // RelativeFactor13
        double dBaseForeignCorrelation, // Correlation12
        double dBaseSpotFXCorrelation, // Correlation13
        double dForeignSpotFXCorrelation) // Correlation23

{
    // Sort Notice Dates for domestic currency bootstrapping
    DKMaille<double> dVolStripDates;
    SelectBootstrappingDates(dVolStripDates,g_dNoticeDates,dStdDevBaseX);
    DKMaille<double> dNewVolStripDates(dVolStripDates.entries()+2);
    DKMaille<double> dSwapStartDates(dVolStripDates.entries());
    DKMaille<double> dSwapEndDates(dVolStripDates.entries());
    dNewVolStripDates.at(0)=dSpotDate;
    double dFinalMaturity=dBoosterData.at(dBoosterData.rows()-1,2);
    for(unsigned int uiD=0;uiD<dVolStripDates.entries();uiD++)
    {
        dNewVolStripDates.at(uiD+1)=dSpotDate+dVolStripDates.at(uiD)*365.;
        // dNoticePeriod is set to 0 days for cash settlement in the bootstrapping procedure. It will be questionable to change it
        dSwapStartDates.at(uiD)=dNewVolStripDates.at(uiD+1)+0.;
        // Only final maturity diagonal is allowed for calibration purposes for the moment
        // Change the below line to make it accept any date
        dSwapEndDates.at(uiD)=dFinalMaturity;
    }
    dNewVolStripDates.at(dNewVolStripDates.entries()-1)=dFinalMaturity;
    DKMaille<double> dModelParameters(9);
    dModelParameters.at(0)=(double)g_uiNumFactors;
    dModelParameters.at(1)=dMeanReversionBase; dModelParameters.at(2)=dMeanReversionForeign; dModelParameters.at(3)=dSpotFXVol;
    dModelParameters.at(4)=dStdDevForeign; dModelParameters.at(5)=dSpotFX; // If fixed income these parameters represent the relative factor
    dModelParameters.at(6)=dBaseForeignCorrelation; dModelParameters.at(7)=dBaseSpotFXCorrelation; dModelParameters.at(8)=dForeignSpotFXCorrelation;

    DKMaille<double> dVolStrip;
    // Check if constant
    double dIsConstant=isConstant(dStdDevBaseZ);
    if(dIsConstant==0.)
    {

        if(g_uiIsDK==0)
        {
            // VFDK for HW
            dVolStrip=Bootstrapping_VFDK_HW1To3F(dStdDevBaseZ,
                                                 dStdDevBaseX,
                                                 dStdDevBaseY,
                                                 dBaseDates,
                                                 dBaseRatesNoBasis,
                                                 dBaseRates,
                                                 dNewVolStripDates,
                                                 dSwapStartDates,
                                                 dSwapEndDates,
                                                 dModelParameters,
                                                 dNewVolStripDates.at(0));
        }
        else
        {
            if(g_uiIsSP==0)
            {
                throw("Unauthorised input !");
            }
            else
            {
                dVolStrip=Bootstrapping_DK3F_Numerical_SP( dStdDevBaseZ,
                          dStdDevBaseX,
                          dStdDevBaseY,
                          dBaseDates,
                          dBaseRatesNoBasis,
                          dBaseRates,
                          dNewVolStripDates,
                          dSwapStartDates,
                          dSwapEndDates,
                          dModelParameters,
                          dNewVolStripDates.at(0),
                          g_dNumSlicesBeforeFirstNotice,
                          g_ppy1,
                          g_ppy,
                          g_dOptimal,
                          0.01); // instead of g_dQParameter_Base
            }
        }


    }
    else
    {
        for(unsigned int uiVol=0;uiVol<dNewVolStripDates.entries()-2;uiVol++)
        {
            dVolStrip.insert(dIsConstant);
        }
    }
    if(dVolStrip.entries()!=dNewVolStripDates.entries()-2)
        throw("Error in input bootstrapping dates");
    DKMaille<double> dOutputStrip(dVolStrip.entries()+2);
    for(uiD=0;uiD<dVolStrip.entries();uiD++)
        dOutputStrip.at(uiD)=dVolStrip.at(uiD);
    // fill the last two holders
    dOutputStrip.at(dVolStrip.entries())=dVolStrip.at(dVolStrip.entries()-1);
    dOutputStrip.at(dVolStrip.entries()+1)=dVolStrip.at(dVolStrip.entries()-1);

    // Write on top of the pre-processed global member variables
    g_S1=dVolStrip.at(0);
    g_MR1=dMeanReversionBase;
    g_S2=dStdDevForeign*(dVolStrip.at(0));
    g_MR2=dMeanReversionForeign;
    g_RHO=dBaseForeignCorrelation;
    g_S3=dSpotFX*(dVolStrip.at(0));
    g_MR3=dSpotFXVol;
    g_RHO1=dBaseSpotFXCorrelation;
    g_RHO2=dForeignSpotFXCorrelation;

    // Write on top of the vol holders
    for(uiD=0;uiD<dNewVolStripDates.entries();uiD++)
        dNewVolStripDates.at(uiD)=(dNewVolStripDates.at(uiD)-dSpotDate)/365.;
    g_dStdDevBaseX=dNewVolStripDates;
    g_dStdDevBaseZ=dOutputStrip;
}


void LatticeFillBootstrappedVolsSingleCurrencyFixedIncome1F(DKMaille<double> dDates,
        DKMaille<double> dRates,
        DKMaille<double> dRatesNoBasis,
        DKMaille<double> dStdDevX,
        DKMaille<double> dStdDevY,
        DKMaille2D<double> dStdDevZ,
        double dSpotDate,
        DKMaille2D<double> dBoosterData,
        double dMeanReversion,
        double dDomesticOrForeign)

{
    // Sort Notice Dates for domestic currency bootstrapping
    DKMaille<double> dVolStripDates;
    SelectBootstrappingDates(dVolStripDates,g_dNoticeDates,dStdDevX);
    DKMaille<double> dNewVolStripDates(dVolStripDates.entries()+2);
    DKMaille<double> dSwapStartDates(dVolStripDates.entries());
    DKMaille<double> dSwapEndDates(dVolStripDates.entries());
    dNewVolStripDates.at(0)=dSpotDate;
    double dFinalMaturity=dBoosterData.at(dBoosterData.rows()-1,2);
    for(unsigned int uiD=0;uiD<dVolStripDates.entries();uiD++)
    {
        dNewVolStripDates.at(uiD+1)=dSpotDate+dVolStripDates.at(uiD)*365.;
        // dNoticePeriod is set to 0 days for cash settlement in the bootstrapping procedure. It will be questionable to change it
        dSwapStartDates.at(uiD)=dNewVolStripDates.at(uiD+1)+0.;
        // Only final maturity diagonal is allowed for calibration purposes for the moment
        // Change the below line to make it accept any date
        dSwapEndDates.at(uiD)=dFinalMaturity;
    }
    dNewVolStripDates.at(dNewVolStripDates.entries()-1)=dFinalMaturity;
    DKMaille<double> dModelParameters(2);
    dModelParameters.at(0)=1.;
    dModelParameters.at(1)=dMeanReversion;
    DKMaille<double> dVolStrip;
    double dIsConstant=isConstant(dStdDevZ);
    if(dIsConstant==0.)
    {
        if(g_uiIsDK==0)
        {
            dVolStrip=Bootstrapping_VFDK_HW1To3F(dStdDevZ,
                                                 dStdDevX,
                                                 dStdDevY,
                                                 dDates,
                                                 dRatesNoBasis,
                                                 dRates,
                                                 dNewVolStripDates,
                                                 dSwapStartDates,
                                                 dSwapEndDates,
                                                 dModelParameters,
                                                 dNewVolStripDates.at(0));
        }
        else
        {
            if(g_uiIsSP==0)
            {
                throw("Unauthorised input !");
            }
            else
            {
                double dSkew=0.;
                if(dDomesticOrForeign==0) dSkew=g_dQParameter_Base;
                if(dDomesticOrForeign==1) dSkew=g_dQParameter_Foreign;
                dVolStrip=Bootstrapping_DK3F_Numerical_SP(dStdDevZ,
                          dStdDevX,
                          dStdDevY,
                          dDates,
                          dRatesNoBasis,
                          dRates,
                          dNewVolStripDates,
                          dSwapStartDates,
                          dSwapEndDates,
                          dModelParameters,
                          dNewVolStripDates.at(0),
                          g_dNumSlicesBeforeFirstNotice,
                          g_ppy1,
                          g_ppy,
                          g_dOptimal,
                          0.01); // instead of dSkew
            }
        }
    }
    else
    {
        for(unsigned int uiVol=0;uiVol<dNewVolStripDates.entries()-2;uiVol++)
        {
            dVolStrip.insert(dIsConstant);
        }
    }
    if(dVolStrip.entries()!=dNewVolStripDates.entries()-2)
        throw("Error in input bootstrapping dates");
    DKMaille<double> dOutputStrip(dVolStrip.entries()+2);
    for(uiD=0;uiD<dVolStrip.entries();uiD++)
        dOutputStrip.at(uiD)=dVolStrip.at(uiD);
    // fill the last two holders
    dOutputStrip.at(dVolStrip.entries())=dVolStrip.at(dVolStrip.entries()-1);
    dOutputStrip.at(dVolStrip.entries()+1)=dVolStrip.at(dVolStrip.entries()-1);

    // Write on top of the vol holders
    for(uiD=0;uiD<dNewVolStripDates.entries();uiD++)
        dNewVolStripDates.at(uiD)=(dNewVolStripDates.at(uiD)-dSpotDate)/365.;

    // Write on top of the pre-processed global member variables
    if(dDomesticOrForeign==0.)
    {
        g_S1=dVolStrip.at(0);
        g_MR1=dMeanReversion;
        g_dStdDevBaseX=dNewVolStripDates;
        g_dStdDevBaseZ=dOutputStrip;
    }
    if(dDomesticOrForeign==1.)
    {
        g_S2=dVolStrip.at(0);
        g_MR2=dMeanReversion;
        g_dStdDevForeignX=dNewVolStripDates;
        g_dStdDevForeignZ=dOutputStrip;
    }


    /// JMP : Extraction of bootstrapped vol to global matrixes
    /// -------------------------------------------------------
    if(dDomesticOrForeign==0.)
    {
        DomVolGlobalExtraction.Resize(dNewVolStripDates.entries(),2);
        for(int i=0;i<dNewVolStripDates.entries();i++)
        {
            DomVolGlobalExtraction.Elt(i,0) = 365.0*g_dStdDevBaseX.at(i);
            DomVolGlobalExtraction.Elt(i,1) = g_dStdDevBaseZ.at(i);
        }
    }
    if(dDomesticOrForeign==1.)
    {
        ForVolGlobalExtraction.Resize(dNewVolStripDates.entries(),2);
        for(int i=0;i<dNewVolStripDates.entries();i++)
        {
            ForVolGlobalExtraction.Elt(i,0) = 365.0*g_dStdDevForeignX.at(i);
            ForVolGlobalExtraction.Elt(i,1) = g_dStdDevForeignZ.at(i);
        }
    }
}

void InterpolationFXVol( DKMaille<double> &dNewVolStripDates,
                         DKMaille<double> &dFXVolDatesPrime,
                         DKMaille<double> &dFXVolPrime,
                         DKMaille<double> &dFXVolDates,
                         DKMaille<double> &dFXVol)
{
    unsigned int uf;

    dFXVolDatesPrime.resize(dNewVolStripDates.entries()-3);
    dFXVolPrime.resize(dNewVolStripDates.entries()-3);


    for(unsigned int ui=1;ui<dNewVolStripDates.entries()-2;ui++)
    {
        dFXVolDatesPrime.at(ui-1)=dNewVolStripDates.at(ui);
        if(dNewVolStripDates.at(ui)<dFXVolDates.at(0))
            dFXVolPrime.at(ui-1)=dFXVol.at(0);
        else if(dNewVolStripDates.at(ui)>dFXVolDates.at(dFXVolDates.entries()-1))
            dFXVolPrime.at(ui-1)=dFXVol.at(dFXVolDates.entries()-1);
        else
        {
            uf=0;
            while(uf<dFXVolDates.entries() && dFXVolDates.at(uf)<dNewVolStripDates.at(ui))
                uf++;
            dFXVolPrime.at(ui-1)=dFXVol.at(uf-1)
                                 +(dNewVolStripDates.at(ui)-dFXVolDates.at(uf-1))/(dFXVolDates.at(uf)-dFXVolDates.at(uf-1))*(dFXVol.at(uf)-dFXVol.at(uf-1));
        }

    }
}


// New function developed by S Pannetier
void LatticeSetSpotFXVol(DKMaille<double> dFXVolDates,
                         DKMaille<double> dFXVol,
                         double dSpotDate,
                         DKMaille2D<double> dBoosterData)

{
    if(dFXVolDates.at(0)==0.)
        throw("Unauthorised first date in spot FX vol");

    if(dFXVol.entries()!=dFXVolDates.entries()+2)
        throw("Unauthorised input for FX volatility");

    //double dFinalMaturity=(dBoosterData.at(dBoosterData.rows()-1,2)-dSpotDate)/365.;

    double dCutOff=dFXVol.at(dFXVol.entries()-2);


    DKMaille<double> dNewVolStripDates;
    DKMaille<double> dNewVols;


    if(dCutOff==0.)
    {
        // no bootstrapping
        // Add anchor to the date strip
        DKMaille<double> dRawSpotFXVol;
        DKMaille<double> dStrip;
        dStrip.insert(0.);
        dRawSpotFXVol.insert(dFXVol.at(0));
        // Long dated historical Spot FX Vol
        double dLongDatedSpotFXVol=dFXVol.at(dFXVol.entries()-1);

        for(unsigned int ui=0;ui<dFXVolDates.entries();ui++)
        {
            dStrip.insert(dFXVolDates.at(ui));
            if(ui<dFXVolDates.entries()-1) dRawSpotFXVol.insert(dFXVol.at(ui+1));
            else dRawSpotFXVol.insert(dLongDatedSpotFXVol);
        }

        dStrip.insert(40.);
        dRawSpotFXVol.insert(dLongDatedSpotFXVol);
        dStrip.insert(50.);
        dRawSpotFXVol.insert(dLongDatedSpotFXVol);

        dNewVolStripDates=dStrip;
        dNewVols=dRawSpotFXVol;

    }
    else
    {
        double dCutOffBis=0.;
        //verification du positionnement du CutOff
        for(unsigned int ui=0;ui<dFXVolDates.entries();ui++)
        {
            if(fabs(dCutOff-dFXVolDates.at(ui))<=0.1)
            {
                dCutOff=dFXVolDates.at(ui);
                if(ui<dFXVolDates.entries()-1) dCutOffBis=dFXVolDates.at(ui+1);
                else dCutOffBis=40.;
                ui=dFXVolDates.entries();
            }
        }

        // si le CutOff n'est pas dans la liste des dates on sort
        if(dCutOffBis==0.)
            throw("Mauvais CutOff");

        // on va calibrer la vol spot sur la premiere date de notice et sur tous les plots de marches qui suivent
        DKMaille<double> dFXVolDatesPrime;
        DKMaille<double> dFXVolPrime;

        //Determination de la premiere Notice Date
        double FirstNotice;
        ui=0;
        while(ui<dBoosterData.rows() && dBoosterData.at(ui,0)<=dSpotDate) ui++;
        // s'il n'y a pas de notice date apres la spotdate, on envoie un message
        if(ui==dBoosterData.rows())
            throw("LatticeSetSpotFXVol : no Notice Date after SpotDate, you should not be here");
        FirstNotice=(dBoosterData.at(ui,0)-dSpotDate)/365.;

        // si la first Notice est avant le plot overnight on ne l'inclut pas
        if(FirstNotice>=dFXVolDates.at(0))
        {
            dFXVolDatesPrime.insert(FirstNotice);
            dFXVolPrime.insert(LinearInterpolation(FirstNotice,dFXVolDates,dFXVol));
        }

        for(ui=0;ui<dFXVolDates.entries();ui++)
        {
            if(dFXVolDates.at(ui)<=dCutOff && dFXVolDates.at(ui)>FirstNotice)
            {
                dFXVolDatesPrime.insert(dFXVolDates.at(ui));
                dFXVolPrime.insert(dFXVol.at(ui));
            }
        }

        DKMaille<double> dStartDate(dFXVolDatesPrime.entries());
        for(ui=0;ui<dFXVolDatesPrime.entries();ui++)
        {
            dStartDate.at(ui)=0.;
        }

        DKMaille<double> dLocalLatticeDates;
        DKMaille<double> dLocalStdDevBaseZ;
        DKMaille<double> dLocalStdDevForeignZ;

        dLocalLatticeDates=g_dDates;
        dLocalLatticeDates.insert(g_dStdDevBaseX.at(g_dStdDevBaseX.entries()-1));

        dLocalStdDevBaseZ=g_dStdDevBaseZ;
        dLocalStdDevForeignZ=g_dStdDevForeignZ;

        TransformVolatilitiesSingle(g_dStdDevBaseX,dLocalLatticeDates, dLocalStdDevBaseZ);
        TransformVolatilitiesSingle(g_dStdDevForeignX,dLocalLatticeDates, dLocalStdDevForeignZ);

        dNewVols=BootstrappingSpotFXVolatility3FTD(0.,
                 dStartDate,
                 dFXVolDatesPrime,
                 dFXVolDatesPrime,
                 dFXVolPrime,
                 dLocalLatticeDates,
                 dLocalLatticeDates,
                 dLocalStdDevBaseZ,
                 dLocalStdDevForeignZ,
                 g_MR1,
                 g_MR2,
                 g_RHO1,
                 g_RHO2,
                 g_RHO);

        // Modification par D Kalafatis & S Pannetier.
        // si la vol spot FX sort negative a cause des donnees de marche
        // renvoyer a l'utilisateur un message d'erreur
        for(unsigned int uiJ=0;uiJ<dNewVols.entries();uiJ++)
            if(dNewVols.at(uiJ)<=0.) throw("Implied Spot FX Vol Is Negative. Check Input FX Vols and Input Swaption Vols");

        dNewVolStripDates.insert(0.);
        for(ui=0;ui<dFXVolDatesPrime.entries();ui++)
            dNewVolStripDates.insert(dFXVolDatesPrime.at(ui));
        dNewVolStripDates.insert(dCutOffBis);
        dNewVolStripDates.insert(50);

        if(fabs(dCutOff-30.)<=0.1||fabs(dFXVol.at(dFXVol.entries()-1))<0.00001)
        {
            double dLastVol=dNewVols.at(dNewVols.entries()-1);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
        }
        else if(fabs(dFXVol.at(dFXVol.entries()-1))<0.011)
        {
            double dLastVol=SpotFXVolIntegral(dNewVolStripDates,dNewVols,0.,dCutOff);
            dLastVol=sqrt(dLastVol/dCutOff);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
        }
        else
        {
            dNewVols.insert(dFXVol.at(dFXVol.entries()-1));
            dNewVols.insert(dFXVol.at(dFXVol.entries()-1));
            dNewVols.insert(dFXVol.at(dFXVol.entries()-1));
        }
    }


    //Write on top of the pre-processed global menber variables
    g_S3=dNewVols.at(0);
    g_dSpotFXVolDatesTD=dNewVolStripDates;
    g_dSpotFXVolTD=dNewVols;


    /// JMP : Extraction of bootstrapped vol to global matrix
    /// -----------------------------------------------------
    FxVolGlobalExtraction.Resize(g_dSpotFXVolDatesTD.entries(),2);
    for(int i=0;i<g_dSpotFXVolDatesTD.entries();i++)
    {
        FxVolGlobalExtraction.Elt(i,0) = 365.0*g_dSpotFXVolDatesTD.at(i);
        FxVolGlobalExtraction.Elt(i,1) = g_dSpotFXVolTD.at(i);
    }
}

// New function developed by S Pannetier
void LatticeSetSpotFXVol_DK(DKMaille<double> dFXVolDates,
                            DKMaille<double> dFXVol,
                            double dSpotDate,
                            DKMaille2D<double> dBoosterData)

{
    if(dFXVolDates.at(0)==0.)
        throw("Unauthorised first date in spot FX vol");

    if(dFXVol.entries()!=dFXVolDates.entries()+2)
        throw("Unauthorised input for FX volatility");

    //double dFinalMaturity=(dBoosterData.at(dBoosterData.rows()-1,2)-dSpotDate)/365.;

    double dCutOff=dFXVol.at(dFXVol.entries()-2);


    DKMaille<double> dNewVolStripDates;
    DKMaille<double> dNewVols;


    if(dCutOff==0.)
    {
        // no bootstrapping
        // Add anchor to the date strip
        DKMaille<double> dRawSpotFXVol;
        DKMaille<double> dStrip;
        dStrip.insert(0.);
        dRawSpotFXVol.insert(dFXVol.at(0));
        // Long dated historical Spot FX Vol
        double dLongDatedSpotFXVol=dFXVol.at(dFXVol.entries()-1);

        for(unsigned int ui=0;ui<dFXVolDates.entries();ui++)
        {
            dStrip.insert(dFXVolDates.at(ui));
            if(ui<dFXVolDates.entries()-1) dRawSpotFXVol.insert(dFXVol.at(ui+1));
            else dRawSpotFXVol.insert(dLongDatedSpotFXVol);
        }

        dStrip.insert(40.);
        dRawSpotFXVol.insert(dLongDatedSpotFXVol);
        dStrip.insert(50.);
        dRawSpotFXVol.insert(dLongDatedSpotFXVol);

        dNewVolStripDates=dStrip;
        dNewVols=dRawSpotFXVol;

    }
    else
    {
        double dCutOffBis=0.;
        //verification du positionnement du CutOff
        for(unsigned int ui=0;ui<dFXVolDates.entries();ui++)
        {
            if(fabs(dCutOff-dFXVolDates.at(ui))<=0.1)
            {
                dCutOff=dFXVolDates.at(ui);
                if(ui<dFXVolDates.entries()-1) dCutOffBis=dFXVolDates.at(ui+1);
                else dCutOffBis=40.;
                ui=dFXVolDates.entries();
            }
        }

        // si le CutOff n'est pas dans la liste des dates on sort
        if(dCutOffBis==0.)
            throw("Mauvais CutOff");

        // on va calibrer la vol spot sur la premiere date de notice et sur tous les plots de marches qui suivent
        DKMaille<double> dFXVolDatesPrime;
        DKMaille<double> dFXVolPrime;

        //Determination de la premiere Notice Date
        double FirstNotice;
        ui=0;
        while(ui<dBoosterData.rows() && dBoosterData.at(ui,0)<=dSpotDate) ui++;
        // s'il n'y a pas de notice date apres la spotdate, on envoie un message
        if(ui==dBoosterData.rows())
            throw("LatticeSetSpotFXVol : no Notice Date after SpotDate, you should not be here");
        FirstNotice=(dBoosterData.at(ui,0)-dSpotDate)/365.;

        // si la first Notice est avant le plot overnight on ne l'inclut pas
        if(FirstNotice>=dFXVolDates.at(0))
        {
            dFXVolDatesPrime.insert(FirstNotice);
            dFXVolPrime.insert(LinearInterpolation(FirstNotice,dFXVolDates,dFXVol));
        }

        for(ui=0;ui<dFXVolDates.entries();ui++)
        {
            if(dFXVolDates.at(ui)<=dCutOff && dFXVolDates.at(ui)>FirstNotice)
            {
                dFXVolDatesPrime.insert(dFXVolDates.at(ui));
                dFXVolPrime.insert(dFXVol.at(ui));
            }
        }

        DKMaille<double> dStartDate(dFXVolDatesPrime.entries());
        for(ui=0;ui<dFXVolDatesPrime.entries();ui++)
        {
            dStartDate.at(ui)=0.;
        }

        DKMaille<double> dLocalLatticeDates;
        DKMaille<double> dLocalStdDevBaseZ;
        DKMaille<double> dLocalStdDevForeignZ;

        dLocalLatticeDates=g_dDates;
        dLocalLatticeDates.insert(g_dStdDevBaseX.at(g_dStdDevBaseX.entries()-1));

        dLocalStdDevBaseZ=g_dStdDevBaseZ;
        dLocalStdDevForeignZ=g_dStdDevForeignZ;

        TransformVolatilitiesSingle(g_dStdDevBaseX,dLocalLatticeDates, dLocalStdDevBaseZ);
        TransformVolatilitiesSingle(g_dStdDevForeignX,dLocalLatticeDates, dLocalStdDevForeignZ);

        dNewVols=BootstrappingSpotFXVolatility3FTD_DK(0.,
                 dStartDate,
                 dFXVolDatesPrime,
                 dFXVolDatesPrime,
                 dFXVolPrime,
                 dLocalLatticeDates,
                 dLocalLatticeDates,
                 dLocalStdDevBaseZ,
                 dLocalStdDevForeignZ,
                 g_MR1,
                 g_MR2,
                 g_RHO1,
                 g_RHO2,
                 g_RHO,
                 g_dCurveDates,
                 g_dCurve,
                 g_dForeignCurveDates,
                 g_dForeignCurve);

        // Modification par D Kalafatis & S Pannetier.
        // si la vol spot FX sort negative a cause des donnees de marche
        // renvoyer a l'utilisateur un message d'erreur
        for(unsigned int uiJ=0;uiJ<dNewVols.entries();uiJ++)
            if(dNewVols.at(uiJ)<=0.) throw("Implied Spot FX Vol Is Negative. Check Input FX Vols and Input Swaption Vols");

        dNewVolStripDates.insert(0.);
        for(ui=0;ui<dFXVolDatesPrime.entries();ui++)
            dNewVolStripDates.insert(dFXVolDatesPrime.at(ui));
        dNewVolStripDates.insert(dCutOffBis);
        dNewVolStripDates.insert(50);

        if(fabs(dCutOff-30.)<=0.1||fabs(dFXVol.at(dFXVol.entries()-1))<0.00001)
        {
            double dLastVol=dNewVols.at(dNewVols.entries()-1);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
        }
        else if(fabs(dFXVol.at(dFXVol.entries()-1))<0.011)
        {
            double dLastVol=SpotFXVolIntegral(dNewVolStripDates,dNewVols,0.,dCutOff);
            dLastVol=sqrt(dLastVol/dCutOff);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
        }
        else
        {
            dNewVols.insert(dFXVol.at(dFXVol.entries()-1));
            dNewVols.insert(dFXVol.at(dFXVol.entries()-1));
            dNewVols.insert(dFXVol.at(dFXVol.entries()-1));
        }
    }


    //Write on top of the pre-processed global menber variables
    g_S3=dNewVols.at(0);
    g_dSpotFXVolDatesTD=dNewVolStripDates;
    g_dSpotFXVolTD=dNewVols;

}


DKMaille<double>  Lattice_HWVFDK_3F_(DKMaille<double> dLatticeGeometryData,
                                     double dNumTimeLinesBeforeFirstNotice,
                                     double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                                     double dNumTimeLinesPerYear,
                                     double dSpotDate,
                                     DKMaille<double> dBaseDates,
                                     DKMaille<double> dBaseRates,
                                     DKMaille<double> dForeignDates,
                                     DKMaille<double> dForeignRates,
                                     DKMaille<double> dRedemptionData,
                                     double dStrike,
                                     double dType,
                                     double dOptionExpiry,
                                     DKMaille<double> dStdDevBaseX,
                                     DKMaille<double> dStdDevBaseY,
                                     DKMaille2D<double> dStdDevBaseZ,
                                     double dMeanReversionBase,
                                     DKMaille<double> dStdDevForeignX,
                                     DKMaille<double> dStdDevForeignY,
                                     DKMaille2D<double> dStdDevForeignZ,
                                     double dMeanReversionForeign,
                                     double dSpotFX,
                                     DKMaille<double> dSpotFXVolDatesTD,
                                     DKMaille<double> dSpotFXVolTD,
                                     double dBaseForeignCorrelation,
                                     double dBaseSpotFXCorrelation,
                                     double dForeignSpotFXCorrelation,
                                     double dProductModelCode,
                                     double dX1Limit,
                                     double dX2Limit,
                                     double dX3Limit,
                                     double dI1Limit,
                                     double dI2Limit,
                                     double dI3Limit,
                                     double dADPLimit,
                                     double dOptimal,
                                     double dTimeBoost,
                                     double dDeltaFlag,
                                     double dStringModel,
                                     double dSmoothing,
                                     double dSurvivalCalculation,
                                     DKMaille2D<double> dBoosterData,
                                     DKMaille<double> dBaseRatesNoBasis,
                                     DKMaille<double> dForeignRatesNoBasis,
                                     double dSmileParameterBase,
                                     double dSmileParameterForeign,
                                     double dIsSwaptionCalibrationWithBasis) 
{


    g_dCurveDates=dBaseDates;
    g_dCurve=dBaseRates;
    g_dForeignCurveDates=dForeignDates;
    g_dForeignCurve=dForeignRates;


    g_uiIsDK=0;
    g_uiIsSP=1;

    if(g_uiIsDK==0)
    {
        dSmileParameterBase=0.;
        dSmileParameterForeign=0.;
    }

    g_dQParameter_Base=dSmileParameterBase;
    g_dQParameter_Foreign=dSmileParameterForeign;

    g_dIsFundingLegStochastic=1.;
    // Tests on data consistency
    // To be written more extensively for final release
    // Currently check that fictive leg is generated properly
    // ie payment date of the fictive funding leg is within 15 days of payment of the FX coupon leg
    for(unsigned int uiCheck=0;uiCheck<dBoosterData.rows();uiCheck++)
    {
        if(fabs(dBoosterData.at(uiCheck,2)-dBoosterData.at(uiCheck,5))>15.)
            throw("Unauthorised lag between fictive funding leg payment date and FX coupon date");
    }

    if(dSurvivalCalculation>0.) g_dSurvivalCalculation=dSurvivalCalculation;
    else g_dSurvivalCalculation=0.;

    if(dRedemptionData.at(0)>0.) g_dRedemption=dRedemptionData.at(0);
    else g_dRedemption=0.;
    if(dRedemptionData.at(1)>0.) g_dRedemptionOption=dRedemptionData.at(1);
    else g_dRedemptionOption=0.;
    if(dRedemptionData.at(2)>0.) g_dRedemptionStrike=dRedemptionData.at(2);
    else g_dRedemptionStrike=0.;
    if(dRedemptionData.at(3)>0.) g_dRedemptionDate=dRedemptionData.at(3);
    else g_dRedemptionDate=0.;
    if(dRedemptionData.at(4)>0.) g_dRedemptionNoticeDate=dRedemptionData.at(4);
    else g_dRedemptionNoticeDate=0.;



    g_dIsSwaptionCalibrationWithBasis=dIsSwaptionCalibrationWithBasis;


    if(dRedemptionData.at(0)>0.&&g_dRedemptionStrike==0.0)
        throw("Unathorised strike input for dual redemption deal");

    // double dStdDevBase=dStdDevBaseZ.at(0,0);
    double dStdDevBase=0.;
    double dStdDevForeign=dStdDevForeignZ.at(0,0);
    // double dStdDevForeign=0.;
    double dSpotFXVol=dSpotFXVolTD.at(0);
    // g_dStdDevBaseZ.resize(dStdDevBaseZ.rows());
    // g_dStdDevForeignX=dStdDevForeignX;
    // g_dStdDevForeignZ.resize(dStdDevForeignZ.rows());
    // for(unsigned int ui=0;ui<dStdDevBaseZ.rows();ui++) g_dStdDevBaseZ.at(ui)=dStdDevBaseZ.at(ui,0);

    // for(unsigned int ui=0;ui<dStdDevForeignZ.rows();ui++) g_dStdDevForeignZ.at(ui)=dStdDevForeignZ.at(ui,0);

    unsigned int ui=0;
    g_dSpotFXVolDatesTD=dSpotFXVolDatesTD;
    g_dSpotFXVolTD=dSpotFXVolTD;


    g_dDeltaFlag=dDeltaFlag;
    g_uiSmoothing=(unsigned int)dSmoothing;
    g_dType=dType;
    g_dOptionExpiry=dOptionExpiry;
    g_dStrike=dStrike;
    g_bStringModel=false;
    if(dStringModel==1.) g_bStringModel=true;
    g_dOptimal=dOptimal;
    g_dTimeBoost=dTimeBoost;
    g_dADPLimit=dADPLimit;
    g_dNumSlicesBeforeFirstNotice=dNumTimeLinesBeforeFirstNotice;
    g_dSpotDate=dSpotDate;

    g_S1=dStdDevBase;
    g_MR1=dMeanReversionBase;
    g_S2=dStdDevForeign;
    g_MR2=dMeanReversionForeign;
    g_RHO=dBaseForeignCorrelation;
    g_S3=dSpotFX;
    g_MR3=dSpotFXVol;
    g_RHO1=dBaseSpotFXCorrelation;
    g_RHO2=dForeignSpotFXCorrelation;

    g_dLimitXX1=dX1Limit;
    g_dLimitXX2=dX2Limit;
    g_dLimitXX3=dX3Limit;

    g_uiI1Limit=(unsigned int)dI1Limit;
    g_uiI2Limit=(unsigned int)dI2Limit;
    g_uiI3Limit=(unsigned int)dI3Limit;

    // Single currency fixed income exotics
    g_uiNumFactors=(unsigned int)dProductModelCode;

    // Dual currency fixed income only
    if(g_uiNumFactors==4) g_uiNumFactors=2;

    // Dual currency fixed income and foreign exchange
    if(g_uiNumFactors==5) g_uiNumFactors=3;

    g_uiProductModelCode=(unsigned int)dProductModelCode;

    if(g_uiProductModelCode==5)
    {
        g_dSpotFX=dSpotFX;
        g_S3=dSpotFXVol;
    }

    g_dDiscountCurveBase.resize(dBaseDates.entries());
    g_dDiscountCurveForeign.resize(dForeignDates.entries());
    g_dDiscountCurveBaseDates.resize(dBaseDates.entries());
    g_dDiscountCurveForeignDates.resize(dForeignDates.entries());
    g_dZCCurveBase.resize(dBaseDates.entries());
    g_dZCCurveForeign.resize(dForeignDates.entries());

    g_dDiscountCurveBaseNoBasis.resize(dBaseDates.entries());

    if(dBaseRates.entries()!=dBaseRatesNoBasis.entries()||dForeignRates.entries()!=dForeignRatesNoBasis.entries())
        throw("Unauthorised attempt for interest rate inputs");

    for(ui=0;ui<dBaseDates.entries();ui++)
    {
        g_dDiscountCurveBase.at(ui)=dBaseRates.at(ui);
        g_dDiscountCurveBaseNoBasis.at(ui)=dBaseRatesNoBasis.at(ui);
        g_dDiscountCurveBaseDates.at(ui)=(dBaseDates.at(ui)-dSpotDate)/365.;
        if(g_dDiscountCurveBaseDates.at(ui)>0.) g_dZCCurveBase.at(ui)=-log(g_dDiscountCurveBase.at(ui))/g_dDiscountCurveBaseDates.at(ui);
        else g_dZCCurveBase.at(ui)=-log(g_dDiscountCurveBase.at(ui+1))/g_dDiscountCurveBaseDates.at(ui+1);
    }

    g_dDiscountCurveForeignNoBasis.resize(dForeignDates.entries());

    for(ui=0;ui<dForeignDates.entries();ui++)
    {
        g_dDiscountCurveForeign.at(ui)=dForeignRates.at(ui);
        g_dDiscountCurveForeignNoBasis.at(ui)=dForeignRatesNoBasis.at(ui);
        g_dDiscountCurveForeignDates.at(ui)=(dForeignDates.at(ui)-dSpotDate)/365.;
        if(g_dDiscountCurveForeignDates.at(ui)>0.) g_dZCCurveForeign.at(ui)=-log(g_dDiscountCurveForeign.at(ui))/g_dDiscountCurveForeignDates.at(ui);
        else g_dZCCurveForeign.at(ui)=-log(g_dDiscountCurveForeign.at(ui+1))/g_dDiscountCurveForeignDates.at(ui+1);
    }





    // Calculate relevant dates for lattice/notice date management
    // if(g_dType<3.0)
    {
        /*
        g_dNoticeDates.resize(1);
        g_dNoticeDates.at(0)=(dNoticeDates.at(dNoticeDates.entries()-1)-dSpotDate)/365.;
        // local date version
        for(ui=0;ui<dNoticeDates.entries();ui++)
        	dNoticeDates.at(ui)=(dNoticeDates.at(ui)-dSpotDate)/365.;
        g_dSpecialDates=dNoticeDates;
        */
    }
    // else
    {
        // It's a non-standard product, ie IR or FX exotic other than European or American FX option
        // First sort the notice dates from the product input
        DKMaille<double> dNoticeDates_(dBoosterData.rows());
        for(unsigned int uiBD=0;uiBD<dBoosterData.rows();uiBD++)
            dNoticeDates_.at(uiBD)=(dBoosterData.at(uiBD,0)-g_dSpotDate)/365.;
        DKMaille<double> dLocalNoticeDates;
        unsigned int uiNoticeSize=0;
        for(uiBD=0;uiBD<dNoticeDates_.entries();uiBD++)
        {
            if(dNoticeDates_.at(uiBD)==0.)
                throw("Notice time is the same as evaluation time. Not authorised");
            // only value the future portion
            if(dNoticeDates_.at(uiBD)>0.)
            {
                if(dNoticeDates_.at(uiBD)!=(0.-g_dSpotDate)/365.)
                {
                    dLocalNoticeDates.insert(dNoticeDates_.at(uiBD));
                    uiNoticeSize++;
                }
            }
        }
        g_dSpecialDates.resize(uiNoticeSize+1);
        g_dSpecialDates.at(0)=0.;
        if(uiNoticeSize==0) throw("There is no notice (call) date left after today. You should use dk_analytics instead of 3FLatticeModel");
        g_dNoticeDates.resize(uiNoticeSize);
        for(uiBD=0;uiBD<uiNoticeSize;uiBD++)
        {
            g_dNoticeDates.at(uiBD)=dLocalNoticeDates.at(uiBD);
            g_dSpecialDates.at(uiBD+1)=g_dNoticeDates.at(uiBD);
        }
    }




    g_uiMinimum=10;
    g_ppy1=(unsigned int)dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal;
    g_ppy=(unsigned int)dNumTimeLinesPerYear;

    g_dTreeCutOffDate=g_dSpecialDates.at(g_dSpecialDates.entries()-1);
    g_dTreeFirstExpiry=g_dSpecialDates.at(1);

    g_uiBranching=3;

    DKMaille<double> dates;

    if(dLatticeGeometryData.at(1)!=dSpotDate) throw("Error Input Lattice Dates");
    if(dLatticeGeometryData.at(0)==0.) CalcLatticeDates(dates,g_uiMinimum);
    else CalcSimpleLatticeDates(dLatticeGeometryData,dates);
    g_uiNumSlices=dates.entries();
    g_dDates=dates;



    DKMaille<double> dShortRateOnSlicesBase;
    DKMaille<double> dShortRateKOnSlicesBase;
    DKMaille<double> dShortRateOnSlicesForeign;
    DKMaille<double> dShortRateKOnSlicesForeign;
    DKMaille<double> dShortRateOnSlices;
    DKMaille<double> dShortRateKOnSlices;
    if(g_uiProductModelCode<4)
    {
        DKMaille<double> dShortRates;
        DKMaille<double> dShortRateDates;
        // Single Currency
        CalcShortRates(dSpotDate,dBaseDates,dBaseRates,dBaseRates,dShortRateDates,dShortRates);
        GetShortRateOnSlice(dShortRateDates,dShortRates,dates,dShortRateOnSlices,dShortRateKOnSlices);
        if(g_uiIsDK==0)
        {
            // HW
            for(unsigned int ui=0;ui<dShortRateOnSlices.entries();ui++)
            {
                dShortRateKOnSlices.at(ui)=1.;
                dShortRateOnSlices.at(ui)=1.;
            }
        }
        g_dEquilibriumForwardBase=dShortRateKOnSlices;
        g_dEquilibriumForwardShortBase=dShortRateOnSlices;
    }
    else
    {
        // Cross Currency
        DKMaille<double> dShortRates;
        DKMaille<double> dShortRateDates;
        DKMaille<double> dShortRatesForeign;
        CalcXCCYShortRates(dSpotDate,
                           dBaseDates,
                           dBaseRates,
                           dBaseRates,
                           dForeignDates,
                           dForeignRates,
                           dForeignRates,
                           dShortRateDates,
                           dShortRates,
                           dShortRatesForeign);
        GetShortRateOnSlice(dShortRateDates,dShortRates,dates,dShortRateOnSlicesBase,dShortRateKOnSlicesBase);
        GetShortRateOnSlice(dShortRateDates,dShortRatesForeign,dates,dShortRateOnSlicesForeign,dShortRateKOnSlicesForeign);
        if(g_uiIsDK==0)
        {
            // HW
            for(unsigned int ui=0;ui<dShortRateKOnSlicesBase.entries();ui++)
            {
                dShortRateKOnSlicesBase.at(ui)=1.;
                dShortRateOnSlicesBase.at(ui)=1.;
                dShortRateKOnSlicesForeign.at(ui)=1.;
                dShortRateOnSlicesForeign.at(ui)=1.;
            }
        }
        g_dEquilibriumForwardBase=dShortRateKOnSlicesBase;
        g_dEquilibriumForwardShortBase=dShortRateOnSlicesBase;
        g_dEquilibriumForwardForeign=dShortRateKOnSlicesForeign;
        g_dEquilibriumForwardShortForeign=dShortRateOnSlicesForeign;

    }


    // Deal with American options, silly case
    if(g_dType>=2.0&&g_dType<3.0)
    {
        g_dNoticeDates.resize(g_uiNumSlices-1);
        for(unsigned int uiNotice=0;uiNotice<g_uiNumSlices-1;uiNotice++)
        {
            // avoid today as a notice date
            g_dNoticeDates.at(uiNotice)=g_dDates.at(uiNotice+1);
        }
    }

    g_dOptionExpiry=(dOptionExpiry-dSpotDate)/365.;


    // Checks on the lattice structure
    if(g_dType==2.0||g_dType==2.5)
    {
        if(g_dDates.at(g_dDates.entries()-1)!=g_dOptionExpiry)
            throw("Lattice must extend to option maturity for American options");
    }

    if(g_dType==3.0||g_dType==4.0||g_dType==5.0||g_dType==6.0||g_dType==4.33)
    {
        if(g_dDates.at(g_dDates.entries()-1)!=g_dNoticeDates.at(g_dNoticeDates.entries()-1))
            throw("Lattice must extend and not exceed to last notice date for PRCS or swaption products");
    }

    if(g_dType<3.) g_dLastDate=(dOptionExpiry-dSpotDate)/365.; // Option Expiry
    else g_dLastDate=(dBoosterData.at(dBoosterData.rows()-1,1)-dSpotDate)/365.; // Last Payment
    // Upgrade this in the future to differentiate between fixed and funding legs

    // If FX option calculate the FX forward at terminal maturity
    if(g_dType<3.)
    {
        double dDomesticDiscount=exp(-g_dLastDate*rateinterpolation_dk_maille(2,g_dLastDate,g_dDiscountCurveBaseDates,
                                     g_dZCCurveBase,g_dDiscountCurveBaseDates.entries()-1));
        double dForeignDiscount=exp(-g_dLastDate*rateinterpolation_dk_maille(2,g_dLastDate,g_dDiscountCurveForeignDates,
                                    g_dZCCurveForeign,g_dDiscountCurveForeignDates.entries()-1));
        g_dFXForwardAtMaturity=g_dSpotFX*dForeignDiscount/dDomesticDiscount;
    }

    g_GridScaling1=2./3.;
    g_GridScaling2=2./3.;
    g_GridScaling3=2./3.;

    if(g_uiProductModelCode<4)
    {
        if(g_dIsSwaptionCalibrationWithBasis==1.0)
        {
            LatticeFillBootstrappedVolsSingleCurrencyFixedIncome(dBaseDates,
                    dBaseRates,
                    dBaseRatesNoBasis,
                    dStdDevBaseX,
                    dStdDevBaseY,
                    dStdDevBaseZ,
                    dSpotDate,
                    dBoosterData,
                    dMeanReversionBase, // MeanReversion1
                    dMeanReversionForeign, // MeanReversion2
                    dSpotFXVol, // MeanReversion3
                    dStdDevForeign, // RelativeFactor12
                    dSpotFX, // RelativeFactor13
                    dBaseForeignCorrelation, // Correlation12
                    dBaseSpotFXCorrelation, // Correlation13
                    dForeignSpotFXCorrelation); // Correlation23
        }
        else if(g_dIsSwaptionCalibrationWithBasis==0.0)
        {
            LatticeFillBootstrappedVolsSingleCurrencyFixedIncome(dBaseDates,
                    dBaseRatesNoBasis,
                    dBaseRatesNoBasis,
                    dStdDevBaseX,
                    dStdDevBaseY,
                    dStdDevBaseZ,
                    dSpotDate,
                    dBoosterData,
                    dMeanReversionBase, // MeanReversion1
                    dMeanReversionForeign, // MeanReversion2
                    dSpotFXVol, // MeanReversion3
                    dStdDevForeign, // RelativeFactor12
                    dSpotFX, // RelativeFactor13
                    dBaseForeignCorrelation, // Correlation12
                    dBaseSpotFXCorrelation, // Correlation13
                    dForeignSpotFXCorrelation); // Correlation23

        }
        else
        {
            throw("Unauthorised Input in Swaption Calibration Method");
        }
    }

    if(g_uiProductModelCode>=4)
    {

        double dDummy=0.;

        if(dIsSwaptionCalibrationWithBasis==1.0)
        {
            LatticeFillBootstrappedVolsSingleCurrencyFixedIncome1F(dBaseDates,
                    dBaseRates,
                    dBaseRatesNoBasis,
                    dStdDevBaseX,
                    dStdDevBaseY,
                    dStdDevBaseZ,
                    dSpotDate,
                    dBoosterData,
                    dMeanReversionBase, // DomesticMeanReversion
                    0.);
            LatticeFillBootstrappedVolsSingleCurrencyFixedIncome1F(dForeignDates,
                    dForeignRates,
                    dForeignRatesNoBasis,
                    dStdDevForeignX,
                    dStdDevForeignY,
                    dStdDevForeignZ,
                    dSpotDate,
                    dBoosterData,
                    dMeanReversionForeign, // ForeignMeanReversion
                    1.);
        }
        else if(dIsSwaptionCalibrationWithBasis==0.0)
        {
            LatticeFillBootstrappedVolsSingleCurrencyFixedIncome1F(dBaseDates,
                    dBaseRatesNoBasis,
                    dBaseRatesNoBasis,
                    dStdDevBaseX,
                    dStdDevBaseY,
                    dStdDevBaseZ,
                    dSpotDate,
                    dBoosterData,
                    dMeanReversionBase, // DomesticMeanReversion
                    0.);
            LatticeFillBootstrappedVolsSingleCurrencyFixedIncome1F(dForeignDates,
                    dForeignRatesNoBasis,
                    dForeignRatesNoBasis,
                    dStdDevForeignX,
                    dStdDevForeignY,
                    dStdDevForeignZ,
                    dSpotDate,
                    dBoosterData,
                    dMeanReversionForeign, // ForeignMeanReversion
                    1.);
        }
        else
        {
            throw("Unauthorised Input in Swaption Calibration Method");
        }


        g_RHO=dBaseForeignCorrelation;

        if(g_uiProductModelCode==5)
        {
            if(g_uiIsDK==0)
            {
                // Desk Version perform bootstrapping up to some date then converge to long term spot FX volatility
                LatticeSetSpotFXVol(dSpotFXVolDatesTD,
                                    dSpotFXVolTD,
                                    dSpotDate,
                                    dBoosterData);
            }
            if(g_uiIsDK==1)
            {
                LatticeSetSpotFXVol_DK(dSpotFXVolDatesTD,
                                       dSpotFXVolTD,
                                       dSpotDate,
                                       dBoosterData);
            }

            g_MR3=dSpotFX; // dummy variable
            g_RHO1=dBaseSpotFXCorrelation;
            g_RHO2=dForeignSpotFXCorrelation;
        }
    }

    // create local version, important for hybrid classes
    DKMaille<double> local_dStdDevBaseZ=g_dStdDevBaseZ;
    DKMaille<double> local_dStdDevBaseX=g_dStdDevBaseX;
    g_dShiftShortGRate.resize(g_uiNumSlices);
    LabelLatticeWithDates();
    // Perform U(3)xU(3)xU(1) projection to calc the drifts in the dual currency and/or the three-factor fx model
    if(g_uiProductModelCode>3)
    {
        // Project into 1D Base dimension to find the base diffusion drift
        U3U3U1_Projection_Base();
        // Project into 1D Foreign dimension to find the foreign diffusion drift
        U3U3U1_Projection_Foreign(dForeignDates.entries());
        // Now move back to three-factor mode
        U3U3U3_Reinstate(dSpotDate,
                         dStdDevBase,
                         dMeanReversionBase,
                         dX1Limit,
                         dBaseDates,
                         dBaseRates,
                         local_dStdDevBaseX,
                         local_dStdDevBaseZ,
                         dProductModelCode,
                         dShortRateKOnSlicesBase,
                         dShortRateOnSlicesBase,
                         dSmileParameterBase);
    }


    //////////////////////////////////////////////////////
    // This is the compact version of the old lattice class //
    //////////////////////////////////////////////////////
    LatticeInitialise();
    LatticeInitialiseOrigin();
    LabelLatticeWithDerivatives();
    InitLatticeModel();
    LabelLatticeLangevinDiffusion();
    SetLatticeMarketDf(0);
    if(g_uiProductModelCode>3) SetLatticeMarketDfForeign(0);
    SetLatticeRateStep(0, 0.0);
    DKMaille<long> indexStates(g_uiNumFactors,0L);
    PrepareTheLatticeNode(0, indexStates);
    LatticeForwardConstructionOfSlices();
    /////////////////////////////////////////

    g_uiVolIsCalculated.resize(g_uiNumSlices);


    // Analytics for FI/FX hybrids
    if(g_uiProductModelCode==5)
    {
        DKMaille<double> dLocalStrip;
        g_dStripDomesticStdDev.resize(g_dStdDevBaseX.entries());
        g_dStripForeignStdDev.resize(g_dStdDevForeignZ.entries());
        g_dStripSpotFXVol.resize(g_dSpotFXVolDatesTD.entries());
        g_dStripDomesticStdDev=g_dStdDevBaseZ;
        g_dStripForeignStdDev=g_dStdDevForeignZ;
        g_dStripSpotFXVol=g_dSpotFXVolTD;
        g_dStrip=CreateNewStripDate();
        TransformVolatilities(g_dStdDevBaseX,
                              g_dStdDevForeignX,
                              g_dSpotFXVolDatesTD,
                              g_dStrip,
                              g_dStripDomesticStdDev,
                              g_dStripForeignStdDev,
                              g_dStripSpotFXVol);

        size_t i;
        VolAnalyticsGlobalExtraction.Resize(g_dStrip.length(),4);
        for(i=0;i<g_dStrip.length();++i)
        {
            VolAnalyticsGlobalExtraction.Elt(i,0)=365.0*g_dStrip.at(i);
            VolAnalyticsGlobalExtraction.Elt(i,1)=g_dStripDomesticStdDev.at(i);
            VolAnalyticsGlobalExtraction.Elt(i,2)=g_dStripForeignStdDev.at(i);
            VolAnalyticsGlobalExtraction.Elt(i,3)=g_dStripSpotFXVol.at(i);
        }

    }
    g_dModelParameters_FixedIncomeExotics.resize(9);
    if(g_uiProductModelCode<4)
    {
        g_dDiscountCurveBase_FixedIncomeExotics.resize(dBaseDates.entries());
        g_dDiscountCurveBaseNoBasis_FixedIncomeExotics.resize(dBaseDates.entries());
        g_dDiscountCurveBaseDates_FixedIncomeExotics.resize(dBaseDates.entries());

        // Reshuffle data for analytics prices
        for(ui=0;ui<dBaseDates.entries();ui++)
        {
            g_dDiscountCurveBase_FixedIncomeExotics.at(ui)=dBaseRates.at(ui);
            g_dDiscountCurveBaseNoBasis_FixedIncomeExotics.at(ui)=dBaseRatesNoBasis.at(ui);
            g_dDiscountCurveBaseDates_FixedIncomeExotics.at(ui)=dBaseDates.at(ui);
        }
        g_dModelParameters_FixedIncomeExotics.at(0)=(double) g_uiNumFactors;
        g_dModelParameters_FixedIncomeExotics.at(1)=dMeanReversionBase;
        g_dModelParameters_FixedIncomeExotics.at(2)=dMeanReversionForeign;
        g_dModelParameters_FixedIncomeExotics.at(3)=dSpotFXVol;
        g_dModelParameters_FixedIncomeExotics.at(4)=dStdDevForeign;
        g_dModelParameters_FixedIncomeExotics.at(5)=dSpotFX;
        g_dModelParameters_FixedIncomeExotics.at(6)=dBaseForeignCorrelation;
        g_dModelParameters_FixedIncomeExotics.at(7)=dBaseSpotFXCorrelation;
        g_dModelParameters_FixedIncomeExotics.at(8)=dForeignSpotFXCorrelation;
        g_dStrip=CreateNewStripDate();
        g_dStripDomesticStdDev.resize(g_dStdDevBaseX.entries());
        g_dStripDomesticStdDev=g_dStdDevBaseZ;
        TransformVolatilitiesSingle(g_dStdDevBaseX,g_dStrip,g_dStripDomesticStdDev);
        g_dStripJulian.resize(g_dStrip.entries());
        for(unsigned int ui=0;ui<g_dStripJulian.entries();ui++)
            g_dStripJulian.at(ui)=g_dSpotDate+g_dStrip.at(ui)*365.;
    }
    g_uiDeltaTick=-1;
    if(g_uiProductModelCode==5)
    {
        LatticeFillFX(g_dSpotFX);
    }
    g_dSpotFXTilde=g_dSpotFX;
    unsigned int uiSizeResult=0;
    if(g_dSurvivalCalculation==0.) uiSizeResult=3;
    if(g_dSurvivalCalculation==1.) uiSizeResult=3+2*g_uiNumSlices;

    static DKMaille<double> dResult;
    dResult.resize(uiSizeResult);

    dResult=Lattice_HWVFDK_3F_Calc(dBoosterData);

    static DKMaille<double> dReturn;
    dReturn.resize(uiSizeResult+1);
    dReturn.at(0)=dResult.at(0);
    dReturn.at(1)=dResult.at(1);
    dReturn.at(2)=dResult.at(2);
    dReturn.at(3)=0.;

    for(unsigned int uiK=4;uiK<uiSizeResult+1;uiK++)
    {
        if(g_dSurvivalCalculation==1.)
        {
            dReturn.at(uiK)=dResult.at(uiK-1);
        }
        else dReturn.at(uiK)=0.;
    }

    if(g_dDeltaFlag==1.)
    {
        g_dSurvivalCalculation=0.;
        double dSpotFXUp=g_dSpotFX+0.01*g_dSpotFX;
        g_dSpotFXTilde=dSpotFXUp;
        LatticeFillFX(dSpotFXUp);
        // Calc with a different SpotFX
        static DKMaille<double> dResultShifted(3);
        dResultShifted=Lattice_HWVFDK_3F_Calc(dBoosterData);
        if(dType>2.5) dReturn.at(3)=(dResultShifted.at(1)+dResultShifted.at(0)-dResult.at(1)-dResult.at(0))/(0.01*g_dSpotFX);
        else dReturn.at(3)=(dResultShifted.at(0)-dResult.at(0))/(0.01*g_dSpotFX);
    }

    deleteAll();
    return dReturn;
}

