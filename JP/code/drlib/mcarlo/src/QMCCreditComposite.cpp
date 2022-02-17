//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : QMCCreditComposite.cpp
//
//   Description : Credit Asset generalizing weighted composition of other credit assets
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCCreditComposite.hpp"
#include "edginc/SRMUtil.hpp"
DRLIB_BEGIN_NAMESPACE

void QMCCreditComposite::ElementHolderBase::setWholeLogSDFRequest()
{
    getCreditAsset()->setWholeTimelineSurvProbRequest();
}

void QMCCreditComposite::ElementHolderBase::setTimelineMapping(
    IQMCHelperTimeLogicSP compositeTimelogic,
    IQMCHelperBoundedDiffusion* diffusionBounds)
{
    const DateTime maxCurveMat      = diffusionBounds->getMaxCurveMat();
    DateTimeArray esdfForwardDates= compositeTimelogic->getFwdEDFDates();
    esdfForwardDates= maxCurveMat.getPastDates(esdfForwardDates);

    bool neededFwdMap = false;
    fwdIdxMap.resize(esdfForwardDates.size());
    for(int idx=0; idx<esdfForwardDates.size(); ++idx)
    {
        fwdIdxMap[idx] = getCreditAsset()->getForwardForwardIndex(esdfForwardDates[idx]);
        QLIB_VERIFY(fwdIdxMap[idx]>=0, "Internal error: composite asset's fwd date is not in compoment's requested list");
        if (idx != fwdIdxMap[idx]) 
            neededFwdMap = true;
    }
    if (!neededFwdMap)
        fwdIdxMap.clear();

    bool neededSpotMap = false;
    const DateTime maxDiffusionDate = diffusionBounds->getMaxDiffDate();
    DateTimeArray sdfRequestedDates = compositeTimelogic->getDFDates();
    sdfRequestedDates = maxDiffusionDate.getPastDates(sdfRequestedDates);

    spotIdxMap.resize(sdfRequestedDates.size());
    for(int idx=0; idx<sdfRequestedDates.size(); ++idx)
    {
        spotIdxMap[idx] = getCreditAsset()->getSpotIndex(sdfRequestedDates[idx]);
        QLIB_VERIFY(spotIdxMap[idx]>=0, "Internal error: composite asset's spot date is not in compoment's requested list");
        if (idx != spotIdxMap[idx]) 
            neededSpotMap = true;
    }
    if (!neededSpotMap)
        spotIdxMap.clear();
}

class ElementHolderFull : public QMCCreditComposite::ElementHolderBase {
public:
    ElementHolderFull(IQMCDiffusibleCreditSpreadBaseSP _pCreditAsset)
        : pCreditAsset(_pCreditAsset) {}
    double getLnFwdValue(FwdIdx i, FwdIdx j)
    {
        return fwdIdxMap.empty() ?
            pCreditAsset->getLnExpectedSurvivalDiscFactor(i,j) :
            pCreditAsset->getLnExpectedSurvivalDiscFactor(fwdIdxMap[i], fwdIdxMap[j]);
    }
    double getSpotValue(SpotIdx i)
    {
        return spotIdxMap.empty() ?
            pCreditAsset->getSurvivalDiscFactor(i) :
            pCreditAsset->getSurvivalDiscFactor(spotIdxMap[i]);
    }
    IQMCDiffusibleCreditSpreadBaseSP getCreditAsset() { return pCreditAsset; }
    double getWholeLogSDF(size_t idx)
    {
        return pCreditAsset->getWholeTimelineLogSurvProb(idx);
    }

private:
    IQMCDiffusibleCreditSpreadBaseSP pCreditAsset;
};

class ElementHolderFractional : public QMCCreditComposite::ElementHolderBase {
public:
    ElementHolderFractional(IQMCCreditSpreadAccessingFractionalSP _pCreditAsset, double _weight)
        : pCreditAsset(_pCreditAsset), weight(_weight) {}
    double getLnFwdValue(FwdIdx i, FwdIdx j)
    {
        return fwdIdxMap.empty() ?
            pCreditAsset->getLnExpectedFractionalSurvivalDiscFactor(i,j,weight) :
            pCreditAsset->getLnExpectedFractionalSurvivalDiscFactor(fwdIdxMap[i], fwdIdxMap[j], weight);
    }
    double getSpotValue(SpotIdx i)
    {
        return std::pow( spotIdxMap.empty() ?
                            pCreditAsset->getSurvivalDiscFactor(i) :
                            pCreditAsset->getSurvivalDiscFactor(spotIdxMap[i]),
                            weight );
    }
    IQMCDiffusibleCreditSpreadBaseSP getCreditAsset() { return pCreditAsset; }
    double getWholeLogSDF(size_t idx)
    {
        return weight * pCreditAsset->getWholeTimelineLogSurvProb(idx);
    }

private:
    IQMCCreditSpreadAccessingFractionalSP pCreditAsset;
    double weight;
};




/** constructor */
QMCCreditComposite::QMCCreditComposite(
        const vector<IQMCDiffusibleCreditSpreadBaseSP>& components,
        const vector<double>&  weights) 
{
    QLIB_VERIFY(components.size() == weights.size(), "sizes of components and weights arrays should be the same");

    for (size_t i=0; i<components.size(); ++i)
    {
        if (weights[i] == 0.0) // TODO: use double-precision epsilon-equality here
            continue;  // do not add any component with weight 0;

        getDiffusionBound()->add(components[i]->getDiffusionBound());

        if (weights[i] == 1.0)
            vComponents.push_back(ElementHolderBaseSP(new ElementHolderFull(components[i])));
        else
        {
            IQMCCreditSpreadAccessingFractionalSP pFractional =
                IQMCCreditSpreadAccessingFractionalSP(
                    dynamic_cast<IQMCCreditSpreadAccessingFractional*>(components[i].get()));
            if (!pFractional.get())
                throw ModelException("QMCCreditComposite::setQMCCreditComposite",
                            "The credit composite component was requested with fractional weight "
                            "but its implementation does not support the fractional access interface. ");
            vComponents.push_back(ElementHolderBaseSP(new ElementHolderFractional(pFractional, weights[i])));
        }
    }
}

/** initialization */
void QMCCreditComposite::setQMCCreditComposite(const DateTime& _today, double _recovery)
{
    today = _today;
    setRecoveryRate(_recovery);
}


/** destructor */
QMCCreditComposite::~QMCCreditComposite()
{}


/** generates path across all dates */
void QMCCreditComposite::generatePathAndSurvivalRates(      
    IQMCRNGManagerSP rngMgr) 
{

    // for each point i on the timeline --
    //   sdf[i] =
    //             pow(pCommonMarket->sdf[i], coeffCommonMarket)
    //          *  pIdiosyncraticSpread->sdf[i]
    //          *  vJumps[0]->sdf[i] * ... * vJumps[k]->sdf[i]
    // this is aggregation for sdf
    for (size_t i=0; i<prob.size(); ++i)
    {
        prob[i] = 1.0;
        for(size_t j=0; j<vComponents.size(); ++j)
            prob[i] *= vComponents[j]->getSpotValue(i);
    }
}

/** finalize the timelines, allocate necessary memory */
void QMCCreditComposite::finalizePathGenerator(DateTimeArrayConstSP allDatesSP)
{
    // need to populate todayIdx, lastDiffusionIdx and mapping arrays in the components
//    timelineSP = allDatesSP;


    if (isWholeTimelineSurvivalRequested())
    {
        calcFirstAndLastDiffusionIdx(*allDatesSP, *allDatesSP);
    }

    DateTimeArray sdfRequestedDates = getSpotDates();

    prob.resize(sdfRequestedDates.size()); //the historic data are passed in, resize to accommodate future ones

//    numESDFDates = 0; // not needed in this model

    // creates necessary mappings
    for(size_t i=0; i<vComponents.size(); ++i)
    {
        vComponents[i]->setTimelineMapping(getTimeLogic(), getDiffusionBound());
    }
}

// setting the flags necessary for DoD calculation
void QMCCreditComposite::setDateOfDefaultEnquiry()
{
    QMCCreditDiffuse::setDateOfDefaultEnquiry();
    setWholeTimelineSurvProbRequest(); // request whole timeline from all the components
}

// requesting the timeline from all the components
void QMCCreditComposite::setWholeTimelineSurvProbRequest()
{
    QMCCreditDiffuse::setWholeTimelineSurvProbRequest();
    for(size_t i=0; i<vComponents.size(); ++i)
        vComponents[i]->setWholeLogSDFRequest();
}

double QMCCreditComposite::getWholeTimelineLogSurvProb(size_t idx)
{
    double val = 0.0;
    for(size_t i=0; i<vComponents.size(); ++i)
        val += vComponents[i]->getWholeLogSDF(idx);
    return val;
}


/** Accessing the natural log of the expected value ExpSDF(md, fd)
    where md is a simulated measurement date and fd is some future
    date after the measurement is made. */
double QMCCreditComposite::getLnExpectedSurvivalDiscFactor(
                                    FwdIdx measurementDateIdx,
                                    FwdIdx futureDateIdx)
{
    // composition for ExpectedSurvivalProb should be done like that:

    //
    //   lnESDF(i,j) =
    //             vFractional[0]->lnESDF(i,j,weight[0]) + ... + vFractional[n]->lnSDF(i,j, weight[n])
    //          +  vFull[0]->lnESDF(i,j)                 + ... + vFull[k]->lnSDF(i,j)

    double lnESDF = 0.0;
    for(size_t i=0; i<vComponents.size(); ++i)
        lnESDF += vComponents[i]->getLnFwdValue(measurementDateIdx, futureDateIdx);
    return lnESDF;
}



/** add all dates that will be needed later */
void QMCCreditComposite::addAggregatedDates(
    const DateTimeArray& _sdfDates,
    const DateTimeArray& _esdfRequestedDates,
    const DateTimeArray& _esdfForwardDates)
{
    // need to propagate the request to the components - to ensure that the components have all the
    // requested quantities available
    QMCCreditDiffuse::addAggregatedDates(_sdfDates, _esdfRequestedDates, _esdfForwardDates);
    for(size_t i=0; i<vComponents.size(); ++i)
        vComponents[i]->getCreditAsset()->addAggregatedDates(_sdfDates, _esdfRequestedDates, _esdfForwardDates);
//    irAsset->addAggregatedDates(_sdfDates, _esdfRequestedDates, _esdfForwardDates);
}


DRLIB_END_NAMESPACE

