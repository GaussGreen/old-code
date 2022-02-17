//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : SVQmcImplemented.cpp
//
//   Description : State variables for accessing information in verious
//                 IQMCDiffusibleAsset objects
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/MCPathConfigSRMGen.hpp"
#include <map>

DRLIB_BEGIN_NAMESPACE


SVQmcDiscFactor::SVQmcDiscFactor(
    IQMCDiffusibleInterestRate*  _pAsset,
    const DateTimeArray&        _measurementDates) :
  IQMCStateVariableSpot(_measurementDates), pAsset(_pAsset)
{

    if (pAsset != NULL)
    {
        pAsset->addSpotDates(measurementDates);
    }
}

double SVQmcPathWeightTimeIndependent::element(int /*idx*/ ) const
    { return pPathConfig->getPathWeight(); }


void SVQmcDiscFactor::prepare(bool mm)
{
    measurementDatesIdx.resize(measurementDates.size());

    DateTime today = pAsset->getBaseDate();

    int beginIdx = 0;
    for(int i=0; i<measurementDates.size(); ++i)
    {
        measurementDatesIdx[i] = pAsset->getSpotIndex(measurementDates[i]);
        if (measurementDates[i] <= today)
            ++beginIdx;
    }

    // keep measurements dates, so we can call prepare() again if pDF changed
    // measurementDates = DateTimeArray(); // dispose memory
	// Q: handle the situation when pDF is different across paths? A: it's impossible
    double* pDF = pAsset->getInternalDFArray();
    if (pDF && !mm) {// important performance boost is available in the model (use old style path object)
        vector<int> v(measurementDatesIdx.size());
        std::copy(measurementDatesIdx.begin(),measurementDatesIdx.end(),v.begin());
        initialize(pDF,v,
                   beginIdx,
                   measurementDatesIdx.size());
    }
    else
        initialize(beginIdx, measurementDatesIdx.size());
}

void SVQmcDiscFactorMM::prepare(bool mm)
{
    SVQmcDiscFactor::prepare(true);
    correction = vector<double>(measurementDatesIdx.size(), 1.0);
    accumulation = vector<double>(measurementDatesIdx.size(), 0.0);
    curvedf.resize(measurementDatesIdx.size(), 0.0);
}


void SVQmcDiscFactorMM::resetMMCorrection()
{
    nAcc = 0;
    correction = vector<double>(measurementDatesIdx.size(), 1.0);
    for (size_t i = 0; i < measurementDatesIdx.size(); ++i)
    {
        curvedf[i] = exp(pAsset->getOriginalLnDiscFactor(measurementDatesIdx[i]));
    }
}

void SVQmcDiscFactorMM::accumulateMMCorrection()
{
    ++nAcc;
    for(size_t i = 0; i < measurementDatesIdx.size(); ++i)
    {
        double df = pAsset->getDiscFactor(measurementDatesIdx[i]);
        accumulation[i] +=  df;
    }
}

void SVQmcDiscFactorMM::setMMCorrection()
{
    for(size_t i = 0; i < measurementDatesIdx.size(); ++i)
    {
        correction[i] = nAcc * curvedf.at(i)/accumulation[i];
   }
}


SVQmcExpectedDiscFactor::SVQmcExpectedDiscFactor(
    IQMCDiffusibleInterestRate*    _pAsset,
    const DateTime&                _measurementDate,
    const DateTimeArray&           _futureDates,
    bool                           _doLog,
    YieldCurveConstSP			   _yc) :
        IQMCStateVariableExpected(_measurementDate, _futureDates, _doLog),
        pAsset(_pAsset), ycIdx((size_t)-1)
{
    if (pAsset != NULL)
    {
        pAsset->addForwardDates(measurementDate, futureDates);
        // select between diffusing and discounting YCs
        ycIdx = pAsset->registerYCFlavor(_yc);
    }
}

void SVQmcExpectedDiscFactor::prepare(bool mm)
{
    measurementDateIdx = pAsset->getForwardForwardIndex(measurementDate);
    futureDatesIdx.resize(futureDates.size());

    for(int i=0; i<futureDates.size(); ++i) {
        futureDatesIdx[i] = pAsset->getForwardForwardIndex(futureDates[i]);
    }
    initialize(0, futureDatesIdx.size());

    //// possibly -- futureDates.clear();
}

SVQmcExpectedDiscFactorMM::SVQmcExpectedDiscFactorMM(
    IQMCDiffusibleInterestRate*    _pAsset,
    const DateTime&                _measurementDate,
    const DateTimeArray&           _futureDates,
    bool                           _doLog,
    YieldCurveConstSP			   _yc) :
    SVQmcExpectedDiscFactor(_pAsset,_measurementDate,_futureDates,_doLog,_yc),
        sumDF(0.0), futureDateFXIdx(-1), futureDateDomDFIdx(-1)
{
    if (pAsset != NULL)
    {
	    IQMCDiffusibleFX* pFX = pAsset->getFXAsset();
	    if (pFX) {
            pFX->addSpotDates(DateTimeArray(1, _measurementDate));
            pFX->getDomesticIRAsset()->addSpotDates(DateTimeArray(1, _measurementDate));
	    }
	    else {
            pAsset->addSpotDates(DateTimeArray(1, _measurementDate));
	    }

    }

}

void SVQmcExpectedDiscFactorMM::prepare(bool mm)
{
    IQMCDiffusibleFX* pFX = pAsset->getFXAsset();
    if (pFX)
    {
        futureDateFXIdx = pFX->getSpotIndex(measurementDate);
        futureDateDomDFIdx = pFX->getDomesticIRAsset()->getSpotIndex(measurementDate);
    }
    else
    {
        futureDateDomDFIdx = pAsset->getSpotIndex(measurementDate);
    }
    if (doLog)
    {
        correction = vector<double>(futureDates.size(), 0.0);
    }
    else
    {
        correction = vector<double>(futureDates.size(), 1.0);
    }

    SVQmcExpectedDiscFactor::prepare(true);
}


void SVQmcExpectedDiscFactorMM::resetMMCorrection()
{
    sumDF = 0.0;
    accumulation = vector<double>(futureDatesIdx.size(), 0.0);
    curveLnEDF.resize(futureDatesIdx.size());

    for (size_t i = 0; i < futureDatesIdx.size(); ++i)
    {
        curveLnEDF[i] = pAsset->getOriginalLnExpectedDiscFactor(ycIdx, measurementDateIdx, futureDatesIdx.at(i));
    }

    if (doLog)
    {
        correction = vector<double>(futureDatesIdx.size(), 0.0);
    }
    else
    {
        correction = vector<double>(futureDatesIdx.size(), 1.0);
    }

}

void SVQmcExpectedDiscFactorMM::accumulateMMCorrection()
{
    IQMCDiffusibleFX* pFX = pAsset->getFXAsset();
    if (pFX)
    {
        IQMCDiffusibleInterestRateSP domPtr = pFX->getDomesticIRAsset();
        double domdf = domPtr->getDiscFactor(futureDateDomDFIdx);
        double fx = pFX->getSpotFX(futureDateFXIdx);
        sumDF += domdf * fx;
        for (size_t i = 0; i < futureDatesIdx.size(); ++i)
        {
            double forEDF = pAsset->getExpectedDiscFactor(ycIdx, measurementDateIdx, futureDatesIdx.at(i));
            accumulation[i] += domdf * fx * forEDF;
        }
    }
    else
    {
        double df = pAsset->getDiscFactor(futureDateDomDFIdx);
        sumDF += df;
        for (size_t i = 0; i < futureDatesIdx.size(); ++i)
        {
            double edf = pAsset->getExpectedDiscFactor(ycIdx, measurementDateIdx, futureDatesIdx.at(i));
            accumulation[i] += df * edf;
        }
    }
}
void SVQmcExpectedDiscFactorMM::setMMCorrection()
{
    if (doLog)
    {
        for (size_t i = 0; i < futureDatesIdx.size(); ++i)
        {
            correction[i] = log(sumDF / accumulation.at(i)) + curveLnEDF.at(i);
        }
    }
    else
    {
        for (size_t i = 0; i < futureDatesIdx.size(); ++i)
        {
            correction[i] = sumDF * exp(curveLnEDF.at(i)) / accumulation.at(i);
        }
    }
}


//////////////////////////////////////////////////////////////////////////

/** Credit-spread specific function names realized in these state variables */


SVQmcSurvivalDiscFactor::SVQmcSurvivalDiscFactor(
    IQMCDiffusibleCreditSpreadBase*  _pAsset,
    const DateTimeArray&        _measurementDates) :
        IQMCStateVariableSpot(_measurementDates),
        pAsset(_pAsset)
{
    if (pAsset)
    {
        pAsset->addSpotDates(measurementDates);
    }
}

void SVQmcSurvivalDiscFactor::prepare(bool mm)
{
    measurementDatesIdx.resize(measurementDates.size());

    DateTime today = pAsset->getBaseDate();

    int beginIdx = 0;
    for(int i=0; i<measurementDates.size(); ++i)
    {
        measurementDatesIdx[i] = pAsset->getSpotIndex(measurementDates[i]);
        if (measurementDates[i] <= today)
         ++beginIdx;
    }

    double* pSDF = pAsset->getInternalSDFArray();
    if (pSDF && !mm) {// important performance boost is available in the model
        vector<int> v(measurementDatesIdx.size());
        std::copy(measurementDatesIdx.begin(),measurementDatesIdx.end(),v.begin());
        initialize(pSDF,v,
                   beginIdx,
                   measurementDatesIdx.size());
    }
    else
        initialize(beginIdx, measurementDatesIdx.size());
}

SVQmcSurvivalDiscFactorMM::SVQmcSurvivalDiscFactorMM(
    IQMCDiffusibleCreditSpreadBase*  _pAsset,
    const DateTimeArray&        _measurementDates) :
        SVQmcSurvivalDiscFactor(_pAsset, _measurementDates)
{
    if (pAsset && pAsset->getUnderlyingIRAsset().get())
        pAsset->getUnderlyingIRAsset()->addSpotDates(measurementDates);
}

void SVQmcSurvivalDiscFactorMM::prepare(bool mm)
{
    correction = vector<double>(measurementDates.size(), 1.0);
    IQMCDiffusibleInterestRateSP irPtr = pAsset->getUnderlyingIRAsset();
    if (irPtr.get())
    {
        measurementDFDatesIdx.resize(measurementDates.size());
        for(int i=0; i<measurementDates.size(); ++i)
        {
            measurementDFDatesIdx[i] = irPtr->getSpotIndex(measurementDates[i]);
        }
    }
    SVQmcSurvivalDiscFactor::prepare(true);
}

void SVQmcSurvivalDiscFactorMM::resetMMCorrection()
{
    correction = vector<double>(measurementDatesIdx.size(), 1.0);
    accumulation.resize(measurementDatesIdx.size(), 0.0);
    accumulationDF.resize(measurementDatesIdx.size(), 0.0);
    curveSDF.resize(measurementDatesIdx.size());
    for (size_t i = 0; i < measurementDatesIdx.size(); ++i)
    {
        curveSDF[i] = exp(pAsset->getOriginalLnSurvivalDiscFactor(measurementDatesIdx[i]));
    }
}
void SVQmcSurvivalDiscFactorMM::accumulateMMCorrection()
{
    IQMCDiffusibleInterestRateSP irPtr = pAsset->getUnderlyingIRAsset();
    for (size_t i = 0; i < measurementDatesIdx.size(); ++i)
    {
        double df =  irPtr.get() ?  irPtr->getDiscFactor(measurementDFDatesIdx.at(i))  : 1.0;
        accumulationDF[i] += df;
        accumulation[i] += df * pAsset->getSurvivalDiscFactor(measurementDatesIdx.at(i));
    }
}

void SVQmcSurvivalDiscFactorMM::setMMCorrection()
{
    for (size_t i = 0; i < measurementDatesIdx.size(); ++i)
    {
        correction[i] = accumulationDF.at(i) * curveSDF.at(i) / accumulation.at(i);
    }
}





SVQmcDateOfDefault::SVQmcDateOfDefault(IQMCDiffusibleDefaultableCreditSpread*  _pAsset)
:
    IQMCStateVariableSpot(DateTimeArray()),
    pAsset(_pAsset)
{
    if (pAsset)
    {
        pAsset->setDateOfDefaultEnquiry();
    }
}

void SVQmcDateOfDefault::prepare(bool mm) {}

/////////////////////////////////////////////////////////////////////////

SVQmcExpSurvDiscFactor::SVQmcExpSurvDiscFactor(
    IQMCDiffusibleCreditSpreadBase*     _pAsset,
    const DateTime&                _measurementDate,
    const DateTimeArray&           _futureDates,
    bool                           _doLog) :
        IQMCStateVariableExpected(_measurementDate, _futureDates, _doLog),
        pAsset(_pAsset)
        {
	        if (pAsset)
            {
	            pAsset->addForwardDates(measurementDate, futureDates);
	        }
        }

void SVQmcExpSurvDiscFactor::prepare(bool mm)
{
    measurementDateIdx = pAsset->getForwardForwardIndex(measurementDate);
    futureDatesIdx.resize(futureDates.size());

    for(int i=0; i<futureDates.size(); ++i)
        futureDatesIdx[i] = pAsset->getForwardForwardIndex(futureDates[i]);

    initialize(0, futureDatesIdx.size());

    //// possibly -- futureDates.clear();
}

double SVQmcExpSurvDiscFactor::getExpRecoveryRate(int idx) const 
{ 
    return pAsset->getExpRecoveryRate(
                            measurementDateIdx,
                            futureDatesIdx[idx]); 
}


SVQmcExpSurvDiscFactorMM::SVQmcExpSurvDiscFactorMM(
    IQMCDiffusibleCreditSpreadBase*     _pAsset,
    const DateTime&                _measurementDate,
    const DateTimeArray&           _futureDates,
    bool                           _doLog) :
        SVQmcExpSurvDiscFactor(_pAsset, _measurementDate, _futureDates, _doLog),
        accumulationDF(0.0), accumulationDFSP(0.0),
        measureDateSDFIdx(-1), measureDateEDFIdx(-1), measureDateDFIdx(-1)
        {
	        if (pAsset)
            {
                pAsset->addSpotDates(DateTimeArray(1, measurementDate)); // for MM
                if (pAsset->getUnderlyingIRAsset().get())
                {
                    pAsset->getUnderlyingIRAsset()->addForwardDates(measurementDate, futureDates);
                    pAsset->getUnderlyingIRAsset()->addSpotDates(DateTimeArray(1, measurementDate));
                }
	        }
        }

void SVQmcExpSurvDiscFactorMM::prepare(bool mm)
{
    IQMCDiffusibleInterestRateSP irPtr = pAsset->getUnderlyingIRAsset();
    measureDateSDFIdx = pAsset->getSpotIndex(measurementDate);
    if(irPtr.get())
    {
        measureDateDFIdx = irPtr->getSpotIndex(measurementDate);
        measureDateEDFIdx = irPtr->getForwardForwardIndex(measurementDate);
        futureDatesDFIdx.resize(futureDates.size());
        for (int i = 0; i < futureDates.size(); ++i)
        {
            futureDatesDFIdx[i] = irPtr->getForwardForwardIndex(futureDates[i]);
        }
    }
    if (!doLog)
    {
        correction = vector<double>(futureDates.size(), 1.0);
    }
    else
    {
        correction = vector<double>(futureDates.size(), 0.0);
    }
    SVQmcExpSurvDiscFactor::prepare(true);

}

void SVQmcExpSurvDiscFactorMM::resetMMCorrection()
{
    accumulationDF = 0.0;
    accumulationDFSP = 0.0;
    accumulation = vector<double>(futureDatesIdx.size(), 0.0);
    accumulationDFEDF = vector<double>(futureDatesIdx.size(), 0.0);
    curveESDF.resize(futureDatesIdx.size());
    for (size_t i = 0; i < futureDatesIdx.size(); ++i)
    {
        curveESDF[i] = pAsset->getOriginalLnExpectedSurvivalDiscFactor(measurementDateIdx, futureDatesIdx.at(i));
    }

    if (!doLog)
    {
        correction = vector<double>(futureDatesIdx.size(), 1.0);
    }
    else
    {
        correction = vector<double>(futureDatesIdx.size(), 0.0);
    }
}

void SVQmcExpSurvDiscFactorMM::accumulateMMCorrection()
{
    IQMCDiffusibleInterestRateSP irPtr = pAsset->getUnderlyingIRAsset();

    double sp = pAsset->getSurvivalDiscFactor(measureDateSDFIdx);
    double df = irPtr.get() ? irPtr->getDiscFactor(measureDateDFIdx) : 1.0 ;

    accumulationDF += df;
    accumulationDFSP += df*sp;
    int ycIdx = 0;

    for (size_t i = 0; i < futureDatesIdx.size(); ++i)
    {
        double edf = irPtr.get() ? irPtr->getExpectedDiscFactor(ycIdx, measureDateEDFIdx, futureDatesDFIdx.at(i)) : 1.0;
        double esp = pAsset->getExpectedSurvivalDiscFactor(measurementDateIdx, futureDatesIdx.at(i));
        accumulationDFEDF[i] += df * edf;
        accumulation[i] += df * sp * edf * esp;
    }
}

void SVQmcExpSurvDiscFactorMM::setMMCorrection()
{
    if (!doLog)
    {
        for (size_t i = 0; i < futureDatesIdx.size(); ++i)
        {
            correction[i] = exp(curveESDF.at(i)) * accumulationDFSP * accumulationDFEDF.at(i)/ (accumulationDF * accumulation.at(i));
        }
    }
    else
    {
        for (size_t i = 0; i < futureDatesIdx.size(); ++i)
        {
            correction[i] = log(accumulationDFSP * accumulationDFEDF.at(i)/ (accumulationDF * accumulation.at(i))) + curveESDF.at(i);
        }
    }
}


SVQmcAggregatedSurvDiscFactor::SVQmcAggregatedSurvDiscFactor(
        IQMCDiffusibleCreditSpreadBase*     _pAsset,
        DateTimeArraySP   _sdfDates,            // set of dates for discount factors
        SpotIdxArraySP    _sdfIdxSP,
        DateTimeArraySP   _esdfRequestedDates,  //  union of {t_i} for all expected factors
        FwdIdxArraySP     _esdfReqIdxSP,
        DateTimeArraySP   _esdfForwardDates,    // union of all {T_j} for all expected factors
        FwdIdxArraySP     _esdfForIdxSP,
        const DateTime&   maxDiffDate,
        const DateTime&   maxCurveMat,
        bool _doLog) :
    pAsset(_pAsset),
    sdfDates(_sdfDates),
    sdfIdxSP(_sdfIdxSP),
    esdfRequestedDates(_esdfRequestedDates),
    esdfReqIdxSP(_esdfReqIdxSP),
    esdfForwardDates(_esdfForwardDates),
    esdfForIdxSP(_esdfForIdxSP)
{
    ASSERT(pAsset);
    pAsset->addAggregatedDates(sdfDates, esdfRequestedDates, esdfForwardDates, maxDiffDate, maxCurveMat);
}

void SVQmcAggregatedSurvDiscFactor::prepare(bool mm) {
    ASSERT(pAsset != NULL);
    DateTime maxDiffDate = pAsset->getDiffusionBound()->getMaxDiffDate();
    DateTime maxCurveMat = pAsset->getDiffusionBound()->getMaxCurveMat();
    //Shared index is the same, but as different assets have different end dates they might not be able to resolve indices
    for(int i= sdfIdxSP->size(); i < sdfDates->size() && (*sdfDates)[i] <= maxDiffDate; ++i)
        sdfIdxSP->push_back(pAsset->getSpotIndex((*sdfDates)[i]));
    // FIXME: precompute loop's upper bound
    for(int i= esdfReqIdxSP->size(); i < esdfRequestedDates->size() && (*esdfRequestedDates)[i] <= maxDiffDate; ++i)
        esdfReqIdxSP->push_back(pAsset->getForwardForwardIndex((*esdfRequestedDates)[i]));
    // FIXME: precompute loop's upper bound
    for(int i= esdfForIdxSP->size(); i < esdfForwardDates->size() && (*esdfForwardDates)[i] <= maxCurveMat; ++i)
        esdfForIdxSP->push_back(pAsset->getForwardForwardIndex((*esdfForwardDates)[i]));
} //



//////////////////////////////////////////////////////////////////////////

/** FX specific function names realized in these state variables */

SVSpotFX::SVSpotFX(
    IQMCDiffusibleFX*            _pAsset,
    const DateTimeArray&        _measurementDates) :
        IQMCStateVariableSpot(_measurementDates),
        pAsset(_pAsset)
{
    if (pAsset)
    {
        pAsset->addSpotDates(measurementDates);
    }
}

void SVSpotFX::prepare(bool mm)
{
    measurementDatesIdx.resize(measurementDates.size());

    DateTime today = pAsset->getBaseDate();

    int beginIdx = 0;
    for(int i=0; i<measurementDates.size(); ++i)
    {
        measurementDatesIdx[i] = pAsset->getSpotIndex(measurementDates[i]);
        if (measurementDates[i] <= today)
         ++beginIdx;
    }

    double* pFX = pAsset->getInternalSpotFXArray();
    if (pFX && !mm) {// important performance boost is available in the model
        vector<int> v(measurementDatesIdx.size());
        std::copy(measurementDatesIdx.begin(),measurementDatesIdx.end(),v.begin());
        initialize(pFX,v,
                   beginIdx,
                   measurementDatesIdx.size());
    }
    else
        initialize(beginIdx, measurementDatesIdx.size());
}

SVSpotFXMM::SVSpotFXMM(
    IQMCDiffusibleFX*            _pAsset,
    const DateTimeArray&        _measurementDates) :
        SVSpotFX(_pAsset, _measurementDates)
{
    if (pAsset)
    {
        pAsset->getDomesticIRAsset()->addSpotDates(measurementDates);
    }
}

void SVSpotFXMM::prepare(bool mm)
{
    measurementDFDatesIdx.resize(measurementDates.size());

    for(int i=0; i<measurementDates.size(); ++i)
    {
        measurementDFDatesIdx[i] = pAsset->getDomesticIRAsset()->getSpotIndex(measurementDates[i]);
    }

    SVSpotFX::prepare(true);
}

void SVSpotFXMM::resetMMCorrection()
{
    correction = vector<double>(measurementDatesIdx.size(), 1.0);
    accumulation = vector<double>(measurementDatesIdx.size(), 0.0);
    accumulationDomDF = vector<double>(measurementDatesIdx.size(), 0.0);

    expFX.resize(measurementDatesIdx.size());

    for (size_t i = 0; i < measurementDatesIdx.size(); ++i)
    {
        expFX[i] = exp(pAsset->getOriginalLnFwdSpot(measurementDatesIdx.at(i)));
    }
}

void SVSpotFXMM::accumulateMMCorrection()
{
    IQMCDiffusibleInterestRateSP domPtr = pAsset->getDomesticIRAsset();
    assert(domPtr.get());
    for(size_t i = 0; i < measurementDatesIdx.size(); ++i)
    {
        double df = domPtr->getDiscFactor(measurementDFDatesIdx[i]);
        accumulationDomDF[i] += df;
        double fx = pAsset->getSpotFX(measurementDatesIdx[i]);
        accumulation[i] +=  fx * df;
    }
}

void SVSpotFXMM::setMMCorrection()
{
    for(size_t i = 0; i < measurementDFDatesIdx.size(); ++i)
    {
        correction[i] = expFX.at(i) *  accumulationDomDF.at(i) / accumulation.at(i);
    }
}

/////////////////////////////////////////////////////////////////////////

SVQmcExpectedFX::SVQmcExpectedFX(
    IQMCDiffusibleFX*               _pAsset,
    const DateTime&                _measurementDate,
    const DateTimeArray&           _futureDates,
    bool                           _doLog) :
        IQMCStateVariableExpected(_measurementDate, _futureDates, _doLog),
        pAsset(_pAsset)
{
  if (pAsset) {
    pAsset->addForwardDates(measurementDate, futureDates);
  }
}


void SVQmcExpectedFX::prepare(bool mm)
{
    measurementDateIdx = pAsset->getForwardForwardIndex(measurementDate);

    futureDatesIdx.resize(futureDates.size());
    for(int i=0; i<futureDates.size(); ++i)
        futureDatesIdx[i] = pAsset->getForwardForwardIndex(futureDates[i]);
    //// possibly -- futureDates.clear();
    initialize(0, futureDatesIdx.size());
}

SVQmcExpectedFXMM::SVQmcExpectedFXMM(
    IQMCDiffusibleFX*              _pAsset,
    const DateTime&                _measurementDate,
    const DateTimeArray&           _futureDates,
    bool                           _doLog) :
        SVQmcExpectedFX(_pAsset, _measurementDate, _futureDates, _doLog), dfDateIdx(-1)
{
  if (pAsset) {
    pAsset->getDomesticIRAsset()->addForwardDates(measurementDate, futureDates);
    pAsset->getForeignIRAsset()->addForwardDates(measurementDate, futureDates);
    pAsset->getDomesticIRAsset()->addSpotDates(DateTimeArray(1,measurementDate));
  }
}


void SVQmcExpectedFXMM::prepare(bool mm)
{
    IQMCDiffusibleInterestRateSP domPtr = pAsset->getDomesticIRAsset();

    originalExpectedFX.resize(futureDates.size());
    edfDomDateIdx.resize(futureDates.size());
    dfDateIdx = domPtr->getSpotIndex(measurementDate);
    edfDomMeasureDateIdx = pAsset->getDomesticIRAsset()->getForwardForwardIndex(measurementDate);
    for (int i = 0; i < futureDates.size(); ++i)
    {
        edfDomDateIdx[i] = pAsset->getDomesticIRAsset()->getForwardForwardIndex(futureDates[i]);
        FwdIdx idx =  pAsset->getForwardForwardIndex(futureDates[i]);
        originalExpectedFX[i] = exp(pAsset->getOriginalLnExpectedFwdSpot(idx));
    }
    if (doLog)
    {
        correction = vector<double>(futureDates.size(), 0.0);
    }
    else
    {
        correction = vector<double>(futureDates.size(), 1.0);
    }

    SVQmcExpectedFX::prepare(true);
}

void SVQmcExpectedFXMM::resetMMCorrection()
{
    accumulation = vector<double>(futureDatesIdx.size(), 0.0);
    accumulationDF = vector<double>(futureDatesIdx.size(), 0.0);

    if (doLog)
    {
        correction = vector<double>(futureDatesIdx.size(), 0.0);
    }
    else
    {
        correction = vector<double>(futureDatesIdx.size(), 1.0);
    }
}

void SVQmcExpectedFXMM::accumulateMMCorrection()
{
    IQMCDiffusibleInterestRateSP domPtr = pAsset->getDomesticIRAsset();
    double df = domPtr->getDiscFactor(dfDateIdx);
    for (size_t i = 0; i < futureDatesIdx.size(); ++i)
    {
        double edf = domPtr->getExpectedDiscFactor(domPtr->getDiscYCIdx(),edfDomMeasureDateIdx, edfDomDateIdx[i]);
        double fwdFX = exp(pAsset->getLnExpectedFX(measurementDateIdx, futureDatesIdx.at(i)));
        accumulation[i] += df * edf * fwdFX;
        accumulationDF[i] += df * edf;
    }
}

void SVQmcExpectedFXMM::setMMCorrection()
{
    for (size_t i = 0; i < futureDatesIdx.size(); ++i)
    {
        correction[i] = originalExpectedFX.at(i) * accumulationDF.at(i) / accumulation.at(i);
    }
    if (doLog)
    {
        for (size_t i = 0; i < futureDatesIdx.size(); ++i)
        {
            correction[i] = log(correction.at(i));
        }
    }
}

//////////////////////////////////////////////////////////////////////////


SVSpotEQ::SVSpotEQ(
    IQMCDiffusibleEQ*            _pAsset,
    const DateTimeArray&        _measurementDates) :
        IQMCStateVariableSpot(_measurementDates) ,
        pAsset(_pAsset)
{
    if (pAsset)
        pAsset->addSpotDates(measurementDates);
}


void SVSpotEQ::prepare(bool mm)
{
    measurementDatesIdx.resize(measurementDates.size());

    DateTime today = pAsset->getBaseDate();

    int beginIdx = 0;
    for(int i=0; i<measurementDates.size(); ++i)
    {
        measurementDatesIdx[i] = pAsset->getSpotIndex(measurementDates[i]);
        if (measurementDates[i] <= today)
         ++beginIdx;
    }

    double* pEQ = pAsset->getInternalSpotPriceArray();
    if (pEQ && !mm) {// important performance boost is available in the model
        vector<int> v(measurementDatesIdx.size());
        std::copy(measurementDatesIdx.begin(),measurementDatesIdx.end(),v.begin());
        initialize(pEQ,v,
                   beginIdx,
                   measurementDatesIdx.size());
    }
    else
        initialize(beginIdx, measurementDatesIdx.size());

}

SVSpotEQMM::SVSpotEQMM(
    IQMCDiffusibleEQ*            _pAsset,
    const DateTimeArray&        _measurementDates) :
        SVSpotEQ(_pAsset, _measurementDates), fxPtr(NULL)
{
   domPtr = pAsset->getDomesticIRAsset();
   if (domPtr.get())
   {
       fxPtr = domPtr->getFXAsset();
       if (fxPtr)
        {
            fxPtr->addSpotDates(_measurementDates);
            domPtr = fxPtr->getDomesticIRAsset();
        }
        domPtr->addSpotDates(_measurementDates);
   }
}


void SVSpotEQMM::prepare(bool mm)
{
    correction = vector<double>(measurementDates.size(), 1.0);
    if (fxPtr)
    {
        measurementFXIdx.resize(measurementDates.size());
        for (int i = 0; i < measurementDates.size(); ++i)
            measurementFXIdx[i] = fxPtr->getSpotIndex(measurementDates[i]);
    }
    if (domPtr.get())
    {
        measurementDFIdx.resize(measurementDates.size());
        for (int i = 0; i < measurementDates.size(); ++i)
            measurementDFIdx[i] = domPtr->getSpotIndex(measurementDates[i]);
    }

    SVSpotEQ::prepare(true);
}

void SVSpotEQMM::resetMMCorrection()
{
    accumulation = vector<double>(measurementDatesIdx.size(), 0.0);
    accumulationDF = vector<double>(measurementDatesIdx.size(), 0.0);
    origFwdPrice.resize(measurementDatesIdx.size());
    for (size_t i = 0; i < measurementDatesIdx.size(); ++i)
        origFwdPrice[i] = exp(pAsset->getOriginalLnFwdSpot(measurementDatesIdx[i]));
    correction = vector<double>(measurementDates.size(), 1.0);
}

void SVSpotEQMM::accumulateMMCorrection()
{
    for (size_t i = 0; i < measurementDatesIdx.size(); ++i)
    {
        double df = domPtr.get() ? domPtr->getDiscFactor(measurementDFIdx.at(i)) : 1.0;
        double fx = fxPtr ? fxPtr->getSpotFX(measurementFXIdx.at(i)) : 1.0;
        double spot = pAsset->getSpotPrice(measurementDatesIdx.at(i));
        accumulation[i] += spot * df * fx;
        accumulationDF[i] += df * fx;
    }
}


void SVSpotEQMM::setMMCorrection()
{
    for (size_t i = 0; i < measurementDatesIdx.size(); ++i)
    {
        correction[i] = origFwdPrice.at(i) * accumulationDF.at(i) / accumulation.at(i);
    }
}


//////////////////////////////////////////////////////////////////////////

SVQmcExpectedEQ::SVQmcExpectedEQ(
    IQMCDiffusibleEQ*               _pAsset,
    const DateTime&                _measurementDate,
    const DateTimeArray&           _futureDates,
    bool                           _doLog) :
        IQMCStateVariableExpected(_measurementDate, _futureDates, _doLog),
        pAsset(_pAsset)
{
    if (pAsset)
    {
        pAsset->addForwardDates(measurementDate, futureDates);
    }
}

void SVQmcExpectedEQ::prepare(bool mm)
{

    measurementDateIdx = pAsset->getForwardForwardIndex(measurementDate);

    futureDatesIdx.resize(futureDates.size());
    for(int i=0; i<futureDates.size(); ++i)
        futureDatesIdx[i] = pAsset->getForwardForwardIndex(futureDates[i]);
    initialize(0, futureDatesIdx.size());
    //// possibly -- futureDates.clear();

}

SVQmcExpectedEQMM::SVQmcExpectedEQMM(
    IQMCDiffusibleEQ*               _pAsset,
    const DateTime&                _measurementDate,
    const DateTimeArray&           _futureDates,
    bool                           _doLog) :
        SVQmcExpectedEQ(_pAsset, _measurementDate, _futureDates, _doLog),
            fxPtr(NULL)
{
   undPtr = pAsset->getDomesticIRAsset(); // should be underlying ccy for given stock
   if (undPtr.get())
   {
       fxPtr = undPtr->getFXAsset();
       if (fxPtr)
        {
            fxPtr->addSpotDates(DateTimeArray(1,_measurementDate));
            domPtr = fxPtr->getDomesticIRAsset(); // should be domestic
        }
       else
           domPtr = undPtr;
       domPtr->addSpotDates(DateTimeArray(1,_measurementDate));
       undPtr->addForwardDates(_measurementDate, _futureDates);
   }
}

void SVQmcExpectedEQMM::prepare(bool mm)
{
    if (doLog) {
        correction = vector<double>(futureDates.size(), 0.0);
    }
    else {
        correction = vector<double>(futureDates.size(), 1.0);
    }

    if (domPtr.get()) dfDomSpotIdx = domPtr->getSpotIndex(measurementDate);
    if (fxPtr) fxSpotIdx = fxPtr->getSpotIndex(measurementDate);
    if (undPtr.get())
    {
        edfForSpotIdx = undPtr->getForwardForwardIndex(measurementDate);
        edfForFwdIdx.resize(futureDates.size());
        for(int i=0; i<futureDates.size(); ++i)
            edfForFwdIdx[i] = undPtr->getForwardForwardIndex(futureDates[i]);
    }

    SVQmcExpectedEQ::prepare(true);
}

void SVQmcExpectedEQMM::resetMMCorrection()
{
    accumulation = vector<double>(futureDatesIdx.size(), 0.0);
    accumulationDF = vector<double>(futureDatesIdx.size(), 0.0);

    origExpectedFwdPrice.resize(futureDatesIdx.size());

    for (size_t i = 0; i < futureDatesIdx.size(); ++i)
    {
        origExpectedFwdPrice.at(i) = exp(pAsset->getOriginalLnExpectedFwdSpot(futureDatesIdx.at(i)));
    }

    if (doLog)
    {
        correction = vector<double>(futureDatesIdx.size(), 0.0);
    }
    else
    {
        correction = vector<double>(futureDatesIdx.size(), 1.0);
    }
}

void SVQmcExpectedEQMM::accumulateMMCorrection()
{
    double df = domPtr.get() ? domPtr->getDiscFactor(dfDomSpotIdx) : 1.0;
    double fx = fxPtr ? fxPtr->getSpotFX(fxSpotIdx) : 1.0;
    for (size_t i = 0; i < futureDatesIdx.size(); ++i)
    {
        double edf = undPtr.get() ?
            undPtr->getExpectedDiscFactor(undPtr->getDiscYCIdx(), edfForSpotIdx, edfForFwdIdx.at(i)) : 1.0;
        double espot = pAsset->getExpectedPrice(measurementDateIdx, futureDatesIdx.at(i)) ;
        accumulation[i] += espot * df * edf * fx;
        accumulationDF[i] += df * edf * fx;
    }
}

void SVQmcExpectedEQMM::setMMCorrection()
{
    for (size_t i = 0; i < futureDatesIdx.size(); ++i)
    {
        correction[i] = origExpectedFwdPrice.at(i) * accumulationDF.at(i) / accumulation.at(i);
    }

    if (doLog)
    {
        for (size_t i = 0; i < correction.size(); ++i)
        {
            correction[i] = log(correction.at(i));
        }
    }
}


//////////////////////////////////////////////////////////////////////////
SVQmcExpEnergyFuture::SVQmcExpEnergyFuture(
    IQMCDiffusibleEnergy * _pAsset,
    const DateTime & _measurementDate,
    const DateTimeArray & _futureDates,
    bool _doLog )
	:
    IQMCStateVariableExpected(_measurementDate, _futureDates, _doLog),
    pAsset(_pAsset),
    prepared(false)
{
	if (pAsset != NULL) {
        pAsset->addForwardDates(_measurementDate, _futureDates);
	}
}


void SVQmcExpEnergyFuture::prepare(bool mm)
{
    assert(!prepared);
    prepared = true;

    measurementDateIdx = (int)pAsset->getForwardForwardIndex(measurementDate);
    futureDatesIdx.resize(futureDates.size());

    for(int i=0; i<futureDates.size(); ++i) {
		futureDatesIdx[i] = pAsset->getForwardForwardIndex(futureDates[i]);
    }
    initialize(0, futureDatesIdx.size());

    //// possibly -- futureDates.clear();
}

SVQmcExpEnergyFutureMM::SVQmcExpEnergyFutureMM(
    IQMCDiffusibleEnergy * _pAsset,
    const DateTime & _measurementDate,
    const DateTimeArray & _futureDates,
    bool _doLog )
	:
    SVQmcExpEnergyFuture(_pAsset, _measurementDate, _futureDates, _doLog), nAcc(0)
{
	if (pAsset != NULL) {
        pAsset->addForwardDates(_measurementDate, _futureDates);
	}
}


void SVQmcExpEnergyFutureMM::prepare(bool mm)
{
    QLIB_VERIFY(false, "Moment Matched energy is not operational at the present moment.");

    if (doLog)
        correction = vector<double>(futureDates.size(), 0.0);
    else
        correction = vector<double>(futureDates.size(), 1.0);

    SVExpEnergyFuture::prepare(true);
}

void SVQmcExpEnergyFutureMM::resetMMCorrection()
{
    nAcc = 0;
    origLnFuturePrice.resize(futureDates.size());

    for (int i = 0; i < futureDates.size(); ++i)
    {
        origLnFuturePrice[i] = 0;  // DanielNg ???
    }

    accumulation = vector<double>(futureDates.size(), 0.0);
    if (doLog)
        correction = vector<double>(futureDates.size(), 0.0);
    else
        correction = vector<double>(futureDates.size(), 1.0);

}
void SVQmcExpEnergyFutureMM::accumulateMMCorrection()
{
    nAcc++;

    for (size_t i = 0; i < futureDatesIdx.size(); ++i)
    {
        accumulation[i] += pAsset->getExpectedPrice(measurementDateIdx, futureDatesIdx[i]);
    }

}

void SVQmcExpEnergyFutureMM::setMMCorrection()
{
    if (doLog)
    {
        for (int i = 0; i < futureDates.size(); ++i)
        {
            correction[i] = log(nAcc / accumulation.at(i)) +  origLnFuturePrice.at(i);
        }
    }
    else
    {
        for (int i = 0; i < futureDates.size(); ++i)
        {
            correction[i] = nAcc / accumulation.at(i) * exp(origLnFuturePrice.at(i));
        }
    }
}

/////////////////////////////////////////////////////////////////////////

SVQmcExpectedBasisFwdSpread::SVQmcExpectedBasisFwdSpread(
    IQMCDiffusibleBasisIndex* _pAsset,
    const DateTime& _measurementDate,
    const DateTimeArray& _resetDate):
        IQMCStateVariableExpected(_measurementDate, _resetDate, false),
        pAsset(_pAsset),
        prepared(false)
{
    if (pAsset != NULL)
    {
        pAsset->addForwardDates(measurementDate, futureDates);
    }

}

///FIXME: double check indices
void SVQmcExpectedBasisFwdSpread::prepare(bool mm)
{
    assert(!prepared);
    measurementDateIdx =
            pAsset->getForwardForwardIndex(measurementDate);
    futureDatesIdx.resize(futureDates.size());

    for(int i=0; i<futureDates.size(); ++i) {
        futureDatesIdx[i] =
                pAsset->getForwardForwardIndex(futureDates[i]);
    }
    initialize(0, futureDatesIdx.size());
    prepared = true;

    //// possibly -- futureDates.clear();
}


/////////////////////////////////////////////////////////////////////////
/*

SVExpectedBasisCoupon::SVExpectedBasisCoupon(
    SVGenIRSwap::IStateVarSP _swapYieldSV, // To do: only defined in cpp!!
    SVExpectedBasisFwdSpreadSP _fwdSpreadSV ) :
        swapYieldSV( _swapYieldSV ),
        fwdSpreadSV( _fwdSpreadSV )
        prepared(false)
{
}

void SVExpectedBasisCoupon::prepare()
{
    assert(!prepared);
    swapYieldSV->prepare();
    fwdSpreadSV->prepare();
    prepared = true;
}

double SVExpectedBasisCoupon::element(int cfIndex) const
{
    assert(prepared);
    double parYield;
    double annuity;
    swapYieldSV->parYield( parYield, annuity );
    double fwdSpread = fwdSpreadSV->getFwdSpread( cfIndex );
    return parYield + fwdSpread;
}
*/


DRLIB_END_NAMESPACE

