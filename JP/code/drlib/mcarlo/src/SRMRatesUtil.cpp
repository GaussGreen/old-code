//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMRatesUtil.cpp
//
//   Description : Helper for SRM - used for holding intermediate data
//
//   Author      : Mark A Robson
//
//   Date        : 14 June 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SRMSwaption.hpp"
#include "edginc/IRCalib.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include <cassert>

DRLIB_BEGIN_NAMESPACE

///// constructs and populates SRMRatesUtil object
SRMRatesUtil::SRMRatesUtil(
    const DateTime&      baseDate,
    int                  numFactors,
    const string&        modelParamsKey,
    const string&        smileParamsKey,
    CVolProcessedSP      _processedVol,
    IYieldCurveConstSP   discYC,
    IYieldCurveConstSP   diffYC,
    bool                 _skipFlag,
    double               _flatVolIr,
    const string&        cutoffChoice,
    double               constCutoffValue,
    const string&        _corrSwapStart, // eg 1Y  (offset to yc spot date)
    const string&        _corrSwapMat,   // eg 10Y (offset to start)
    const string&        _corrSwapDCC,   // eg Act/365F
    const string&        _corrSwapFreq): // eg 6M
        QMCRatesUtil(baseDate, discYC, diffYC),
        corrSwapDCC(DayCountConventionFactory::make(_corrSwapDCC)),
        corrSwapFreq(_corrSwapFreq), 
        cutoffChoice(cutoffChoice), 
        constCutoffValue(constCutoffValue),
        processedVol(_processedVol),
        flatVolIr(_flatVolIr),
        momentMatching(false),
        skipFlag(_skipFlag)
{
        // then process data for swap to correlate against other IR's
        MaturityTimePeriod start(_corrSwapStart, baseDate.getTime());
        MaturityTimePeriod mat(_corrSwapMat, baseDate.getTime());
        corrSwapStart = start.toDate(diffYC->getSpotDate());
        corrSwapMat = mat.toDate(corrSwapStart);
}


/** Constructs and populates SRMRatesUtil object suitable using the 'VF' style
    calibration. qLeft and qRight are hard coded to 1, as is alpha. Note
    that only single factor IR is supported */
SRMRatesUtil::SRMRatesUtil(
    const DateTime&      baseDate,
    FXAssetConstSP       fx, /* optional - if specified will add FX smile
                                dates to timeline */
    const string&        modelParamsKey,
    CVolProcessedSP      processedVol,
    IYieldCurveConstSP   discYC,
    IYieldCurveConstSP   diffYC,
    bool                 skipFlag,
    const string&        _corrSwapStart, // eg 1Y  (offset to today)
    const string&        _corrSwapMat,   // eg 10Y (offset to start)
    const string&        _corrSwapDCC,   // eg Act/365F
    const string&        _corrSwapFreq): // eg 6M
        QMCRatesUtil(baseDate, discYC, diffYC),
        corrSwapDCC(DayCountConventionFactory::make(_corrSwapDCC)),
        corrSwapFreq(_corrSwapFreq)
{ 
        initialized = true; // dates are known

        // then process data for swap to correlate against other IR's
        MaturityTimePeriod start(_corrSwapStart, baseDate.getTime());
        MaturityTimePeriod mat(_corrSwapMat, baseDate.getTime());
        corrSwapStart = start.toDate(diffYC->getSpotDate());
        corrSwapMat = mat.toDate(corrSwapStart);
}


/*****  From irdiffuse::IR_SpotVol   ****************************************/
void SRMRatesUtil::spotVol(double flatVolIr) 
{
    // this case (as invoked by srm3 code) collapses to this
    SpotVol = vector<double>(extendedTimeLine.size()-1, flatVolIr);
    swaptionSpotVolAtSimDates = vector<double>((*dates).size()-1, flatVolIr);
}

//// for ease until we move everything over
void SRMRatesUtil::extendSpotVol(const DoubleArray& SrcSpotVol) 
{
    vector<double> tmp(SrcSpotVol.begin(), SrcSpotVol.end());
    this->extendSpotVol(tmp);
}

void SRMRatesUtil::extendSpotVol(const vector<double>& SrcSpotVol) 
{
    SpotVol = SRMUtil::extendVol(swaptionExpiries,
                        SrcSpotVol, 
                        extendedTimeLine);
}

/** wrapper around extendSpotVol to extend spot vols from this object on
    supplied (*dates) */
vector<double> SRMRatesUtil::extendSpotVol(const DateTimeArray& newTimePoints) const
{
    return SRMUtil::extendVol(swaptionExpiries, SpotVol, newTimePoints);
}

DRLIB_END_NAMESPACE
