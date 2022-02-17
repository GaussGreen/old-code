//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : QMCCreditDiffuse.cpp
//
//   Description : Base class for SRM credit diffusion
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCCreditDiffuse.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/MemoryProfiler.hpp"
#include "edginc/Format.hpp"
#include "edginc/QMCHelperCachingTimeLogic.hpp"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/SmallLinearInterpolator.hpp"
#include "edginc/QMCHelperDateTimeCache.hpp"
#include "edginc/IQMCRNGManager.hpp"

#include <limits>
#include <fstream>

DRLIB_BEGIN_NAMESPACE

/** constructor */
QMCCreditDiffuse::QMCCreditDiffuse() :
    // default initialization with 'bad' values
    isDateOfDefaultRequested(false),
    isWholeTimelineRequested(false),
    probStart(-999),
    todayIdx(-999),
    lastDiffusionIdx(-999),
    recoveryRate(0.0),
    timeLogic(new QMCHelperCachingTimeLogic(IQMCHelperDateTimeCacheSP(new QMCHelperDateTimeCache),
              IQMCHelperDateTimeCacheSP(new QMCHelperDateTimeCache),
              IQMCHelperDateTimeCacheSP(new QMCHelperDateTimeCache))
             )
{}

/** destructor */
QMCCreditDiffuse::~QMCCreditDiffuse() {
}

/** get the simulation start date */
DateTime QMCCreditDiffuse::getBaseDate() {
    return today;
}

/** accessing the diffused SDF, i.e. a non-default probability by a given date */
double QMCCreditDiffuse::getSurvivalDiscFactor(SpotIdx measurementDateIdx) {
    if (isDateOfDefaultRequested)
       throw ModelException("QMCCreditDiffuse::getSurvivalDiscFactor",
            "After DateOfDefault requested - no Survival Probability shall be returned.");

    return prob[measurementDateIdx];
}

/** this is for efficient access to SDF values*/
double* QMCCreditDiffuse::getInternalSDFArray() {
    return (prob.empty() || isDateOfDefaultRequested) ? NULL : &prob[0];
}


/** returns IQMCHelperTimeLogic object */
IQMCHelperTimeLogicSP QMCCreditDiffuse::getTimeLogic() const {
    return timeLogic;
}

/** returns date on the timeline */
DateTime QMCCreditDiffuse::getTimeLineDate(int idx) {
    return (*timelineSP)[idx];
}

/** returns log of deterministic surv prob */
double QMCCreditDiffuse::getOriginalLnSurvivalDiscFactor(SpotIdx measurementDateIdx)
{
    if (originalProbs.empty())
        throw ModelException("QMCCreditDiffuse::getOriginalLnSurvivalDiscFactor", "originalProbs is not initialized");
    return originalProbs.at(measurementDateIdx);
}

/** returns log of deterministic forward surv prob */
double QMCCreditDiffuse::getOriginalLnExpectedSurvivalDiscFactor(FwdIdx measurementDateIdx, FwdIdx futureDateIdx)
{
    if (originalExpProbs.empty())
        throw ModelException("QMCCreditDiffuse::getOriginalLnExpectedSurvivalDiscFactor", "originalExpProbs is not initialized");
    return originalExpProbs.at(futureDateIdx) - originalExpProbs.at(measurementDateIdx);
}

void QMCCreditDiffuse::generateDefaultDate(IQMCRNGManagerSP rngMgr)
{
    IUniformRNGGenSP uniform = rngMgr->getSharedGen();
    double u = uniform->fetch();

//    ofstream ftest("e:/temp/defdatedump.txt", ios::app);
//
//    ftest<<" ------------------------- "<<endl;
//    ftest<<" years         logsurvprob "<<endl;
//    for(size_t z=0; z+1<m_wholeTimeLineInYears.size(); ++z)
//    {
//        ftest<<m_wholeTimeLineInYears[z]<<"       "<<getWholeTimelineLogSurvProb(z)<<endl;
//
//    }
//    ftest<<m_wholeTimeLineInYears.back()<<endl;
//    ftest.close();

    double _crdEvTime = smallLinearInterpolation(
                make_begin_iterator(this, &QMCCreditDiffuse::getWholeTimelineLogSurvProb),
                make_end_iterator(this, &QMCCreditDiffuse::getWholeTimelineLogSurvProb, m_wholeTimeLineInYears.size()-1),
                m_wholeTimeLineInYears.begin(),
                m_wholeTimeLineInYears.end(),
                -log(u));

    dateOfDefault = getBaseDate().rollDate(int(0.5 + 365 * _crdEvTime));
}

/* FIXME: normally we're passed here allDates (i.e. with past dates), but here we use allDates as simDates */

void QMCCreditDiffuse::finalizeWholeTimelineRequest(DateTimeArrayConstSP simDates)
{

    QLIB_VERIFY(todayIdx >= 0, "todayIdx has to be initialized to the position of today in the full set of simDates.");
    QLIB_VERIFY(lastDiffusionIdx >= todayIdx, "lastDiffusionIdx has to be initialized to the position of the max diffusion maturity in the full set of simDates.");
    m_wholeTimeLineInYears.resize(lastDiffusionIdx-todayIdx+1, 0.0);

    for(int i=todayIdx; i<=lastDiffusionIdx; ++i)
        m_wholeTimeLineInYears[i-todayIdx] = today.yearFrac((*simDates)[i]);

    m_wholeTimeLineInYears.push_back(m_wholeTimeLineInYears.back() + 1.0); 
    // value to be returned if no default happened in the observed time
}

void QMCCreditDiffuse::calcFirstAndLastDiffusionIdx(const DateTimeArray& simDates, const DateTimeArray& allDates)
{
    todayIdx = today.find(allDates);

    int idx = -1;
    const int numSimDates = simDates.size();
    if (! getDiffusionBound()->empty() && getDiffusionBound()->getMaxDiffDate() != DateTime()) {
        const DateTime maxDiffusionDate = getDiffusionBound()->getMaxDiffDate();

        idx    = maxDiffusionDate.findLower(simDates);
        QLIB_VERIFY(idx >= 0 && idx < numSimDates, "Unexpected value for the lastDiffusionIdx");
    }
    lastDiffusionIdx = idx;
}


void QMCCreditDiffuse::processAllDates(DateTimeArrayConstSP allDatesSP)
{

    timelineSP = allDatesSP;
//    todayIdx = today.find((*allDatesSP));

    DateTimeArray sdfRequestedDates =   getSpotDates();
    DateTimeArray esdfRequestedDates =  getForwardDates();

    // Consistency checks
    ASSERT(DateTime::isSubset((*allDatesSP), sdfRequestedDates));
    ASSERT(DateTime::isSubset((*allDatesSP), esdfRequestedDates));

    // this->numESDFDates = esdfRequestedDates.size();

     // translation of SDF dates
    probIndexes = DateTime::getIndexes((*allDatesSP),  sdfRequestedDates); // [Ndf]

    // translation of ESDF dates
    expProbIndexes = DateTime::getIndexes((*allDatesSP), esdfRequestedDates); // [Nedf]

    // make life easier - add request for index off the end
    // the index has to be positive or stepIdx will not work.
    const int IMPOSSIBLE_POSITIVE_INDEX = numeric_limits<int>::max();
    probIndexes.push_back(IMPOSSIBLE_POSITIVE_INDEX);
    expProbIndexes.push_back(IMPOSSIBLE_POSITIVE_INDEX);

    // determine probStart; NOTE today is a part of the past
    probStart = today.findUpper(sdfRequestedDates);
    if (probStart >= 0 && probStart != sdfRequestedDates.size() && sdfRequestedDates[probStart] == today)
        ++probStart; // skip today as diffusion results are always in futureDates

    prob.resize(sdfRequestedDates.size()); //the historic data are passed in, resize to accommodate future ones
}

double QMCCreditDiffuse::getWholeTimelineLogSurvProb(size_t i)
{
    throw ModelException(
        "QMCCreditDiffuse::getWholeTimelineLogSurvProb",
        "this method was not properly defined in the derived class");
}







DRLIB_END_NAMESPACE
