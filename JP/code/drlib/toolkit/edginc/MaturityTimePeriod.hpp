//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MaturityTimePeriod.hpp
//
//   Description : Defines floating expiries used to define yield curve & 
//                 vol surface points e.g. 1M, 5Y with a fixed time of day
//
//   Author      : Andrew J Swain
//
//   Date        : 15 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef MATURITYTIMEPERIOD_HPP
#define MATURITYTIMEPERIOD_HPP

#include "edginc/Expiry.hpp"
#include "edginc/MaturityPeriod.hpp"

DRLIB_BEGIN_NAMESPACE

/** Defines floating expiries used to define yield curve & 
    vol surface points e.g. 1M, 5Y with a fixed time of day
    e.g. roll forward 6M but use the time in the period, not
    the time in the day we're advancing from
    MaturityTimePeriods are immutable - see explanation
    in Expiry.hpp.
    DO NOT ADD ANY NON CONST METHODS */

class TOOLKIT_DLL  MaturityTimePeriod: public Expiry {
public:
    static CClassConstSP const TYPE;

    MaturityTimePeriod(const string& period, int time);

    virtual ~MaturityTimePeriod();

    /** Explicit implementation for performance */
    virtual IObject* clone() const;

    /** overrides default */
    virtual void validatePop2Object();

    /** Returns true if given object matches this */
    bool equalTo(const IObject* expiry) const;

    /** returns a hashcode for the object based upon the period and time.
        Overridden for performance */
    int hashCode() const;

    /** return this expiry as a string */
    virtual string toString() const;

    /** return this expiry as an absolute date */
    virtual DateTime toDate(const DateTime& date) const;

    /** write object out to writer */
    virtual void write(const string& tag, Writer* writer) const;
    /** populate an empty object in from reader  */
    virtual void import(Reader::Node* elem, Reader* reader); 

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** Returns true if given expiry matches this */
    virtual bool equals(const Expiry* expiry) const;

    /* specific functions to MaturityTimePeriod */
    /** Returns the MaturityPeriod portion of the MaturityTimePeriod */
    MaturityPeriodConstSP getMaturityPeriod() const;

    /** Returns the maturity portion of the MaturityTimePeriod */
    string  getMaturity() const;

    /** Returns the time portion of the MaturityTimePeriod */
    int  getTime() const;
private:
    friend class MaturityTimePeriodHelper;
    MaturityTimePeriod();
    MaturityTimePeriod(const MaturityTimePeriod& matTime);
    // fields
    const MaturityPeriodConstSP period;
    const int                   time;
};

typedef smartConstPtr<MaturityTimePeriod> MaturityTimePeriodConstSP;
typedef smartPtr<MaturityTimePeriod> MaturityTimePeriodSP;

DRLIB_END_NAMESPACE
#endif
