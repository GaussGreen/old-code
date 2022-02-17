//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenExpectedDiscFactor.hpp
//
//   Description : A Generator of MC EXPECTED Discount Factor State Variables
//                 aka Zero Coupon Bond Prices
//
//   Date        : 6 Aug 2004
//
//
//----------------------------------------------------------------------------

#ifndef SVGenExpectedDiscFactor_HPP
#define SVGenExpectedDiscFactor_HPP

//#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/ISVBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** A type-specific interface to the MC Expected Discount Factor State Variable */
class MCARLO_DLL SVExpectedDiscFactor : public virtual ISVBase
{
public:
    virtual ~SVExpectedDiscFactor() {}
    virtual double getFwdDF( int idx ) const  = 0;
    virtual double firstDF() const  = 0;
};
DECLARE(SVExpectedDiscFactor);

/** A Generator of MC Expected Discount Factor State Variables.
    It's important to understand the difference between regular discount
    factors (as per SVGenDiscFactor) and expected discount factors.
    Here we are saying, at a given date in the simulation, what is the 
    expected discount factor at that date between two future dates.
    
    For expected discount factors we (for the moment) fail if the calc date 
    is not the pvDate.
 */
class MCARLO_DLL SVGenExpectedDiscFactor: virtual public IElemStateVariableGen,
                            virtual public VirtualDestructorBase
{
public:

/** Note: The SVGenExpectedDiscFactor::IStateVar has been obsoleted by the 
    SVExpectedDiscFactor class.   */

    /** Constructor - for computing expected discount values on calcDate between
        pvDate and each date in dates. Note that currently past values for this
        are not supported. In particular a value of 0 is returned when calcDate
        is in the past (even if computeLog is true) */
    SVGenExpectedDiscFactor(const DateTime&      calcDate, // when to compute
                         const DateTime&      pvDate, // when to discount to
                         YieldCurveConstSP    yieldCurve,
                         const DateTimeArray& dates,
                         bool                 computeLog); // true: do log

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface). The
        previous IStateVariableSP (may be null) should be passed in.  The
        return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP              /*oldStateVar*/,
        IStateVariableGen::IStateGen*  pathGen) const 
    {
        return getSVExpectedDiscFactor(pathGen);
    }

    /** Returns a state variable which then provides access to the values etc. This
        is the method that products should call to get an SVExpectedDiscFactor. */
    SVExpectedDiscFactorSP getSVExpectedDiscFactor(
        IStateVariableGen::IStateGen* pathGen) const;

    /** Retrieve the yield curve associated with this SVGenExpectedDiscFactor */
    YieldCurveConstSP getYieldCurve() const;

    /** Retrieve the dates for which a discount factor is required */
    const DateTimeArray& getDates() const;
    /** Returns the date to which we pv to */
    const DateTime&   getPVDate() const;
    /** Returns the date on which the expected value should be computed */
    const DateTime&   getCalcDate() const;
    /** is the log of the expected discount factor wanted */
    bool logRequired() const;
    /** For use by Path Generators (past or future) that want to use
        deterministic rates */
   SVExpectedDiscFactor*  determinsticSV(const DateTime& today,
                                        bool            doingPast) const;
   void attachSVGen(IElemStateVariableGenVisitor*) const;
private:
    class DeterminsticSV;
    friend class DeterminsticSV;
    DateTime          calcDate; // when to compute expected value
    DateTime          pvDate; // when to discount to
    DateTimeArray     dates;
    YieldCurveConstSP yieldCurve;
    bool              computeLog;
};

typedef smartPtr<SVGenExpectedDiscFactor> SVGenExpectedDiscFactorSP;
//DECLARE(SVGenExpectedDiscFactor); 

DRLIB_END_NAMESPACE

#endif
