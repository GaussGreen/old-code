//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : QMCCreditComposite.hpp
//
//   Description : Base class for SRM credit path generation
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_SRMCREDITCIDCOMPOSITE_HPP
#define QLIB_SRMCREDITCIDCOMPOSITE_HPP

#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/QMCCreditDiffuse.hpp"
#include "edginc/SRMRatesHJMDiffuse.hpp"
#include "edginc/SRMRatesHJMUtil.hpp"
#include "edginc/DECLARE.hpp"
#include <set>

#include "edginc/SRMCreditUtil.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/QMCHelperBoundedDiffusion.hpp"

//#include "edginc/QMCCreditCIDJumps.hpp"
//#include "edginc/QMCCreditCIDDiffuse.hpp"
DRLIB_BEGIN_NAMESPACE

// Work In Progress

/** base class for SRM credit path generation */
class QMCCreditComposite : public QMCCreditDiffuse {
public:

    /** provides remapping of date indices if necessary and encapsulates Full/Fractional function calls */
    class ElementHolderBase : public VirtualDestructorBase {
    public:
        virtual ~ElementHolderBase() {}
        virtual double getLnFwdValue(FwdIdx i, FwdIdx j) = 0;
        virtual double getSpotValue(SpotIdx i) = 0;
        virtual IQMCDiffusibleCreditSpreadBaseSP getCreditAsset() = 0;
        virtual double getWholeLogSDF(size_t idx) = 0;
        virtual void setWholeLogSDFRequest();
        void setTimelineMapping(
            IQMCHelperTimeLogicSP compositeTimelogic,
            IQMCHelperBoundedDiffusion* diffusionBounds);

    protected:
        vector<FwdIdx>  fwdIdxMap;
        vector<SpotIdx> spotIdxMap;
    };
    DECLARE(ElementHolderBase);

    /** constructor */
    QMCCreditComposite(
        const vector<IQMCDiffusibleCreditSpreadBaseSP>& components,
        const vector<double>&  weights);

    /** initialization */
    void setQMCCreditComposite(const DateTime& _today, double recovery);


    /** destructor */
    virtual ~QMCCreditComposite();


    /** add all dates that will be needed later */
    virtual void addAggregatedDates(
        const DateTimeArray& _sdfDates,
        const DateTimeArray& _esdfRequestedDates,
        const DateTimeArray& _esdfForwardDates);

    /** return underlying IR asset */
    virtual IQMCDiffusibleInterestRateSP getUnderlyingIRAsset() { return IQMCDiffusibleInterestRateSP(); }

    /** generates path across all dates */
    virtual void generatePathAndSurvivalRates(      
        IQMCRNGManagerSP rngMgr);

    /** finalize the timelines, allocate necessary memory */
    virtual void finalizePathGenerator(DateTimeArrayConstSP allDates);

    /** returns pointer to internal storage for MaxDiff(usion)Date / MaxCurveMaturity */
    virtual QMCHelperBoundedDiffusion * getDiffusionBound() { return &diffBound; }

// TO be activated later
//    virtual double getRecoveryRate(SpotIdx idx)
//    {
//        return recoveryRate; // TODO: implement the proper calculation with "catastophic recovery"
//    }
//    virtual double getExpRecoveryRate(
//                                        FwdIdx measurementDateIdx,
//                                        FwdIdx futureDateIdx) 
//    {
//        return recoveryRate; // base model does not have time-dep RR at the moment
//    }


    /** Accessing the expected value ExpSDF(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double  getExpectedSurvivalDiscFactor(
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx)
    { return std::exp(QMCCreditComposite::getLnExpectedSurvivalDiscFactor(measurementDateIdx, futureDateIdx));}

    /** Accessing the natural log of the expected value ExpSDF(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedSurvivalDiscFactor(
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx);


    /** Informs asset that the date of default will be asked for a path */
    /** using this method shall make it illegal to call SDF for this asset */
    virtual void setDateOfDefaultEnquiry();

    virtual void setWholeTimelineSurvProbRequest();

    virtual double getWholeTimelineLogSurvProb(size_t idx);

protected:

    QMCCreditComposite();
    vector<ElementHolderBaseSP>        vComponents;
    QMCCRBoundedDiffusion diffBound;
};

/** declare smartPtr versions */
DECLARE(QMCCreditComposite);

DRLIB_END_NAMESPACE

#endif // QLIB_SRMCREDITCIDCOMPOSITE_HPP
