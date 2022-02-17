//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : QMCRatesUtil.hpp
//
//   Description : base class for IR util classes (calibration etc)
//
//----------------------------------------------------------------------------

#ifndef QMCRATESUTIL_HPP
#define QMCRATESUTIL_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DECLARE.hpp"

#include "edginc/YieldCurve.hpp"
#include "edginc/VolProcessedBSIR.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

// This class knows nothing about the volatilities or factors -- its prime purpose is to 
// declare the date logic, the yield curve logic, and some other minor things

/** Helper class for models to get data out of market cache */        

class MCARLO_DLL QMCRatesUtil : public virtual VirtualDestructorBase
{

public:
    
    ///// constructs and populates SRMRatesUtil object
    QMCRatesUtil(
            const DateTime&      baseDate,
            IYieldCurveConstSP   discYC,
            IYieldCurveConstSP   diffYC); // eg 6M
           
    //// computes log of DF [determinstic] from discount zero curve
    //// between dates on extended time line and populates
    //// logDiscFactor (and sets it to right length)
    void computeLogDiscFactor(vector<double>& logDiscFactor) const;
    void computeLogDiscFactor(
    						const DateTimeArray& myDates,
    						vector<double>& logDiscFactor,
    						IYieldCurveConstSP yc) const;

    /** Calculates a deterministic discount factor from date to today. Uses
        diffused curve if useDiffusedCurve is true else uses discount curve */
    double pv(const DateTime& date, bool useDiffusedCurve) const;

    /** Calculates a deterministic discount factor between the 2 dates. Uses
        diffused curve if useDiffusedCurve is true else uses discount curve */
    double pv(const DateTime& loDate, const DateTime& hiDate,
              bool useDiffusedCurve) const;


    /// Accessors (inline)
    
    /** returns the discount YC */
    IYieldCurveConstSP getDiscYC() const{ return discYC; }
    /** returns the YC that is diffused */
    IYieldCurveConstSP getDiffYC() const{ return diffYC; }
    /** Get the ISO code for this IR factor */
    string getCcy() const{ return discYC->getCcy(); }
    /** returns all the simulation dates excluding today */
    const DateTimeArray& getSimDates() const{ return (*dates);}
    /** returns zero dates merged with sim date */
    const DateTimeArray& getExtendedTimeLine() const{assertInitialized(); return extendedTimeLine; }
    /** returns today */
    const DateTime& getBaseDate() const{ return baseDate; }


	/** number of factors in model */
    // called when we know simDates
    virtual void setTimeLine(DateTimeArrayConstSP simDates) = 0;

protected:

    void calcExtendedTimeLine();
    
    DateTime  baseDate; // today

    DateTimeArrayConstSP dates; // sim start + dates after sim start
    
    DateTimeArray   extendedTimeLine; // dates + extra dates from diffuse
    vector<double>  FwdRate; // on extendedTimeLine
    IYieldCurveConstSP   discYC;
    IYieldCurveConstSP   diffYC;
  
    bool 				 initialized; // set when allDates are known
    void 				 assertInitialized(void) const;

};

DECLARE(QMCRatesUtil);

DRLIB_END_NAMESPACE
#endif
