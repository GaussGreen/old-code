//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenSurvivalDiscFactor.hpp
//
//   Description : A Generator of MC Survival Discount Factor State Variables
//
//   Author      : Eva X Strasser
//
//
//----------------------------------------------------------------------------

#ifndef SVGenSurvivalDiscFactor_HPP
#define SVGenSurvivalDiscFactor_HPP

#include "edginc/CDSParSpreads.hpp"
#include "edginc/ISVBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** A type-specific interface to the MC Discount Factor State Variable */
class MCARLO_DLL SVSurvivalDiscFactor : public virtual ISVBase
{
public:
    virtual ~SVSurvivalDiscFactor() {}
    virtual double firstSDF() const = 0;
    virtual double getSDF(int idx) const = 0;
    virtual double getRecoveryRate(int idx) const = 0;
};
DECLARE(SVSurvivalDiscFactor);


/** A Generator of MC Survival Discount Factor State Variables. Can return simulated
    survival discount factors between today and a set of dates. The set of dates can
    be specified in a variety of manners */

class MCARLO_DLL SVGenSurvivalDiscFactor:  virtual public IElemStateVariableGen,
                                        public virtual VirtualDestructorBase
{
public:
 
    /** Note: The SVGenSurvivalDiscFactor::IStateVar has been obsoleted by the
        class MCARLO_DLL SVSurvivalDiscFactor.   */

    /** Constructor - from  an array of dates.
        For computing discount factors between today and each date in dates */
    SVGenSurvivalDiscFactor(const DateTime&        today,
                         ICDSParSpreadsConstSP  cdsParSpreadCurve,
                         const DateTimeArray&   maturityDates);

    /** Constructor - from a single date */
    SVGenSurvivalDiscFactor(const DateTime&        today,
                         ICDSParSpreadsConstSP  cdsParSpreadCurve,
                         const DateTime&        maturityDate);

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface).
        The previous IStateVariableSP (may be null) should be passed in.
        The return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP                oldStateVar,
                                    IStateVariableGen::IStateGen*   pathGen) const;

    /** Returns a SVGenSurvivalDiscFactor state variable which then
        provides access to the path etc. This is the method that
        products should call to get an SVSurvivalDiscFactorSP. */
    SVSurvivalDiscFactorSP getSVSurvivalDiscFactor(IStateVariableGen::IStateGen* pathGen) const;

    /** Retrieve the CDS par spread curve associated with this
        SVGenSurvivalDiscFactor */
    ICDSParSpreadsConstSP getCDSParSpreadCurve() const;

    /** Retrieve the dates for which a survival discount factor is required */
    const DateTimeArray&   getDates() const;

    /** For use by Path Generators (past or future) that want to use determinstic rates */
    SVSurvivalDiscFactor* determinsticSV(bool doingPast) const;
    void attachSVGen(IElemStateVariableGenVisitor*) const;

private:
    /** Basic validation. Should be called by all constructors after
        population of all fields */
    void validate();
    class DeterminsticSV;
    DateTime                today;              // when to discount to
    ICDSParSpreadsConstSP   cdsParSpreadCurve;
    DateTimeArray           maturityDates;      // dates at which survival discount factors are requested
};

typedef smartPtr<SVGenSurvivalDiscFactor> SVGenSurvivalDiscFactorSP;
typedef vector<const SVGenSurvivalDiscFactor*> SVGenSurvivalDiscFactorArray;
// DECLARE(SVGenSurvivalDiscFactor);



DRLIB_END_NAMESPACE

#endif
