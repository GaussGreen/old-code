//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : QMCCreditDiffuse.hpp
//
//   Description : Base class for credit path generation
//
//
//----------------------------------------------------------------------------

#ifndef QMCCREDITDIFFUSE_HPP
#define QMCCREDITDIFFUSE_HPP

#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/DECLARE.hpp"
#include <set>

#include "edginc/QMCHelperBoundedDiffusion.hpp"

DRLIB_BEGIN_NAMESPACE

/** base class for SRM credit path generation */
class QMCCreditDiffuse : public IQMCDiffusibleDefaultableCreditSpread {
public:

    /** constructor */
    QMCCreditDiffuse();

    /** destructor */
    virtual ~QMCCreditDiffuse();

    virtual void setRecoveryRate(double _recoveryRate) { recoveryRate = _recoveryRate; }

    /** get the simulation start date */
    virtual DateTime getBaseDate();

    /** accessing the diffused SDF, i.e. a non-default probability by a given date */
    virtual double getSurvivalDiscFactor(SpotIdx measurementDateIdx);

    /** this is for efficient access to SDF values stored for all sdf requested dates*/
    virtual double* getInternalSDFArray();

    /** returns IQMCHelperTimeLogic object */
    virtual IQMCHelperTimeLogicSP getTimeLogic() const;

    /** returns date on the timeline */
    DateTime getTimeLineDate(int idx);

    // DanielNg: need to check the refererence date
    /** returns log of deterministic surv prob */
    virtual double getOriginalLnSurvivalDiscFactor(SpotIdx measurementDateIdx);

    /** returns log of deterministic forward surv prob */
    virtual double getOriginalLnExpectedSurvivalDiscFactor(FwdIdx measurementDateIdx, FwdIdx futureDateIdx);

    /** Informs asset that the date of default will be asked for a path */
    /** using this method shall make it illegal to call SDF for this asset */
    virtual void setDateOfDefaultEnquiry() { isDateOfDefaultRequested = true; isWholeTimelineRequested=true; }

    /** Retrieves the simulated date of default */
    virtual DateTime getDateOfDefault()
    {
        if (!isDateOfDefaultRequested)
            throw ModelException("QMCCreditDiffuse::getDateOfDefault", "DateOfDefault not prepared.");

        return dateOfDefault;
    }

    /** Retrieves the recovery rate as of date of default */
    virtual double getRecoveryRateAtDefault() { return 0.0; }

    /** getting the simulated recovery rates */
    virtual double getRecoveryRate(SpotIdx idx)
    {
        return recoveryRate; // base model does not have time-dep RR at the moment
    }

    virtual double getExpRecoveryRate(
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx) 
    {
        return recoveryRate; // base model does not have time-dep RR at the moment
    }

    virtual void generatePathAndSurvivalRates(
        IQMCRNGManagerSP rngMgr) = 0;

    virtual void finalizePathGenerator(DateTimeArrayConstSP allDates) = 0;

    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)  // do not extend this function in derived classes!
                                                        // work with generatePathAndSurvivalRates instead
    { 
        generatePathAndSurvivalRates(rngMgr);
        if (isDateOfDefaultRequested)
            generateDefaultDate(rngMgr);
    }

    virtual void finalize(DateTimeArrayConstSP allDates) // do not extend this function in derived classes!
                                          // work with finalizePathGenerator instead
    {
        finalizePathGenerator(allDates);
        if (isWholeTimelineRequested)
            finalizeWholeTimelineRequest(allDates);
    }



    virtual void generateDefaultDate(IQMCRNGManagerSP rngMgr);
    virtual void finalizeWholeTimelineRequest(DateTimeArrayConstSP simDates);

    virtual void setWholeTimelineSurvProbRequest() { isWholeTimelineRequested = true; }

    virtual double getWholeTimelineLogSurvProb(size_t i);

protected:
    // purposefully non-virtual -- the whole information should stay at this level
    bool isWholeTimelineSurvivalRequested() const { return isWholeTimelineRequested; }

    DateTime        dateOfDefault;
    double          recoveryRate;

protected:
    // At this level we get committed to a timeline: this includes notions of allDates, simDates, today and various indexes/
    DateTime        today;
    int             todayIdx;
    int             lastDiffusionIdx;   //last index of simDates we do diffusion for, set in processSimDates() (called from finalize()) based on maxMatDate
    // cached variables
    vector<int>     probIndexes;        // indexes for when we save survival probabilities [Nsdf] set in processAllDates()
    vector<int>     expProbIndexes;     // indexes for when we compute expected survival probabilities [Nesdf] set in processAllDates()

    vector<double>  originalProbs;         // [Ndf], used only by moment matching as reference values
    vector<double>  originalExpProbs; 
    // computed discount factors at probIndexes
    vector<double>  prob;               // at the beginning historic lnQ, later historic lnQ + simulated lnQ
    int             probStart;          /* where to start in prob vector, where ...
                                            at the beginning:   prob vector = historic lnQ
                                            at the end:         prob vector = historic lnQ + simulated lnQ */

    void   calcFirstAndLastDiffusionIdx(const DateTimeArray& simDates, const DateTimeArray& allDates); // sets todayIdx and last required diffusion index
    void   processAllDates(DateTimeArrayConstSP allDates); // during finalize() initialize things that depend solely on allDates

private:
    IQMCHelperTimeLogicSP   timeLogic;      // aux object to manage DateTimeArray calculations
    DateTimeArrayConstSP    timelineSP;     // remember our timeline


private:
    bool            isDateOfDefaultRequested;
    bool            isWholeTimelineRequested;

    vector<double>    m_wholeTimeLineInYears;

//    // something to store the implied recovery rates array
//    // they are calculated only once
//    std::vector<double>  regularRecoveryRates;
//    std::vector<double>  catastrophicRecoveryRates;
//
//    DateTime             switchRecoveryRateDate;   // path-dependent

};

/** declare smartPtr versions */
DECLARE(QMCCreditDiffuse);

DRLIB_END_NAMESPACE

#endif // QMCCREDITDIFFUSE_HPP
