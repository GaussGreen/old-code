//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditCIRUtil.hpp
//
//   Description : Derived SRMCreditUtil class for CIR model
//
//
//----------------------------------------------------------------------------

#ifndef SRMCREDITCIRUTIL_HPP
#define SRMCREDITCIRUTIL_HPP

#include "edginc/SRMCreditUtil.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DECLARE.hpp"

#include <cassert>

DRLIB_BEGIN_NAMESPACE

class SRMCreditCIRUtil : public SRMCreditUtil
{
public:

    /** constructor */
    SRMCreditCIRUtil(const DateTime&       baseDate,
                     const string&         smileParams,
                     ICDSParSpreadsConstSP stochCDSCurve); 
           
    /** destructor */
    virtual ~SRMCreditCIRUtil();

    /** returns spot vols */
    double getSpotVol() const;

    /** returns beta */
    double getBeta() const;

    /** returns qRight */
    double getQRight() const;

    /** returns qLeft */
    double getQLeft() const;

    /** finishes initialization */
    virtual void setTimeLine(DateTimeArrayConstSP simDates);


    bool getMomentMatchedFlag() {return false;}

private:

    double                  beta;       // mean reversion speed
    double                  SpotVol;    // spot vols
    double                  qLeft;      // smile parameter
    double                  qRight;     // smile parameter
};

/** declare smartPtr versions */
DECLARE(SRMCreditCIRUtil);

DRLIB_END_NAMESPACE

#endif // SRMCREDITCIRUTIL_HPP
