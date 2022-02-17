//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditUtil.hpp
//
//   Description : Base class for SRM Credit util classes- calibration etc
//
//
//----------------------------------------------------------------------------

#ifndef SRMCREDITUTIL_HPP
#define SRMCREDITUTIL_HPP

#include "edginc/CDSParSpreads.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DECLARE.hpp"

#include <cassert>

DRLIB_BEGIN_NAMESPACE

class SRMCreditUtil : public virtual VirtualDestructorBase {
public:

    /** constructor */
    SRMCreditUtil(
        const DateTime& baseDate,
        const string& smileParams,
        ICDSParSpreadsConstSP stochCDSCurve); 

    /** destructor */
    virtual ~SRMCreditUtil();

    /** returns log of the deterministic survival probability between two dates
    i.e. returns log [S(0,T) / S(0,t)] = - int_t^T lambda(0,u) du */
    double logSurvProbRatio(const DateTime& loDate, const DateTime& hiDate) const;

    /** returns log of the deterministic survival probability from today to a given date, 
    i.e. returns log S(0,t) = -int_0^t lambda(0,u) du */
    double logSurvProbRatio(const DateTime& date) const;

    /** computes and stores log of deterministic ratio of survival probabilities from first date in a 
    vector of given dates to each date (populating logFwdProbSimple and setting it to right length)
    result_i = log [S(0,T_i) / S(0,T_0)] = - int_{T_0}^{T_i} lambda(0,u) du */
    void computeLogFwdProbSimple(const DateTimeArray& myDates,
                                 vector<double>& fwdLnProbRatios) const;

    /** computes and stores log of the deterministic ratio of consecutive survival probabilities
    on extended time line (populating logFwdProbSimple and setting it to right length)
    result_i = log [S(0,T_{i+1}) / S(0,T_{i})] = - int_{T_{i}}^{T_{i+1}} lambda(0,u) du */
    void computeLogFwdProbSimple(vector<double>& logFwdProbSimple) const;

    /** returns the CDS curve that is diffused */
    ICDSParSpreadsConstSP getCdsCurve() const;

    /** returns all the simulation dates excluding today */
    const DateTimeArray& getSimDates() const;

    /** returns CDS dates merged with sim date */
    const DateTimeArray& getExtendedTimeLine() const;

    /** returns today */
    const DateTime& getBaseDate() const;

    /** finishes initialization */
    virtual void setTimeLine(DateTimeArrayConstSP simDates) = 0;
    
    void setMomentMatchingFlag(bool mm) {momentMatching = mm;}
    bool getMomentMatchingFlag() const {return momentMatching;}

protected:

    /** populates extendedTimeLine field 
    dates (date of simStart plus dates after simStart) + extra dates from input 
    along the lines of SRM3:ExtendFullTimeLine */
    void calcExtendedTimeLine();

    DateTime                baseDate;           // today
    DateTimeArrayConstSP    dates;              // sim start + dates after sim start // shared between CR
    DateTimeArray           extendedTimeLine;   // dates + extra dates from diffuse
    ICDSParSpreadsConstSP   stochCDSCurve;      // CDS curve
    bool                    momentMatching;
    bool                    initialized;        // true after setTimeLine() is called
};

/** declare smartPtr versions */
DECLARE(SRMCreditUtil);

DRLIB_END_NAMESPACE
#endif // SRMCREDITUTIL_HPP
