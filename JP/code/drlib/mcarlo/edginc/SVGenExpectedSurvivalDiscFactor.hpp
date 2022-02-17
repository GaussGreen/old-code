//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenExpectedSurvivalDiscFactor.hpp
//
//   Description : A Generator of MC EXPECTED Survival Discount Factor State Variables
//
//   Author      : Eva X Strasser
//
//
//----------------------------------------------------------------------------

#ifndef MCExpectedSurvivalDiscFactor_HPP
#define MCExpectedSurvivalDiscFactor_HPP

#include "edginc/CDSParSpreads.hpp"
#include "edginc/ISVBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** A type-specific interface to the MC Discount Factor State Variable */
class MCARLO_DLL SVExpSurvDiscFactor : public virtual ISVBase
{
public:
    virtual ~SVExpSurvDiscFactor() {}
    virtual double firstExpSDF() const = 0;
    virtual double getExpSDF(int idx) const = 0;
    virtual double getExpRecoveryRate(int idx) const = 0;
};
DECLARE(SVExpSurvDiscFactor);

/** A Generator of MC Expected Survival Discount Factor State Variables.
    It's important to understand the difference between regular discount
    factors (as per SVGenSurvivalDiscFactor) and expected survival discount factors.
    Here we are saying, at a given date in the simulation, what is the 
    expected survival discount factor at that date between two future dates.
    
    For expected survival discount factors we (for the moment) fail if the calc date is not the pvDate.
 */
class MCARLO_DLL SVGenExpectedSurvivalDiscFactor: 
    virtual public IElemStateVariableGen,
    virtual public VirtualDestructorBase {
public:
    /** Note: The SVGenExpectedSurvivalDiscFactor::IStateVar has been obsoleted by the 
        class MCARLO_DLL IQSVGenExpectedSurvivalDiscFactorSV.   */

    /** Constructor - for computing expected survival discount values on calcDate between
        pvDate and each date in dates. Note that currently past values for this
        are not supported. In particular a value of 0 is returned when calcDate
        is in the past (even if computeLog is true) */
    SVGenExpectedSurvivalDiscFactor(const DateTime&        calcDate,   // when to compute
                                 const DateTime&        pvDate,     // when to discount to
                                 ICDSParSpreadsConstSP  cdsParSpreadCurve,
                                 const DateTimeArray&   dates,
                                 bool                   computeLog); // true: do log

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface). The
        previous IStateVariableSP (may be null) should be passed in.  
        The return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP                oldStateVar,
                                    IStateVariableGen::IStateGen*   pathGen) const;

    /** Returns a Expected Survival Discount Factor state variable which then
        provides access to the values etc. This is the method that
        products should call to get an SVExpSurvDiscFactor. */
    
    SVExpSurvDiscFactorSP getSVExpectedSurvivalDiscFactor(
        IStateVariableGen::IStateGen* pathGen) const;

    /** Retrieve the CDS par spread curve associated with this SVGenExpectedSurvivalDiscFactor */
    ICDSParSpreadsConstSP getCDSParSpreadCurve() const;

    /** Retrieve the dates for which a discount factor is required */
    const DateTimeArray& getDates() const;
    
    /** Returns the date to which we pv to */
    const DateTime& getPVDate() const;
    
    /** Returns the date on which the expected value should be computed */
    const DateTime& getCalcDate() const;
    
    /** is the log of the expected discount factor wanted */
    bool logRequired() const;
    
    /** name of underlying credit */ 
    const string getName() const;
    
    /** For use by Path Generators (past or future) that want to use
        determinstic rates */
    SVExpSurvDiscFactor* determinsticSV(
        const DateTime& today,
        bool            doingPast) const;
    virtual void attachSVGen(IElemStateVariableGenVisitor*) const;    // support visitor pattern
private:
    class DeterminsticSV;
    friend class DeterminsticSV;
    DateTime                calcDate;           // when to compute expected value
    DateTime                pvDate;             // when to discount to
    DateTimeArray           dates;
    ICDSParSpreadsConstSP   cdsParSpreadCurve;
    bool                    computeLog;
};

DECLARE(SVGenExpectedSurvivalDiscFactor);

DRLIB_END_NAMESPACE

#endif
