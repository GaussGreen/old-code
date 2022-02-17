//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenDateOfDefault.hpp
//
//   Description : A Generator of MC Default Date State Variables
//
//   Author      : Lawrence Siu
//
//
//----------------------------------------------------------------------------

#ifndef SVGenDateOfDefault_HPP
#define SVGenDateOfDefault_HPP

#include "edginc/CDSParSpreads.hpp"
#include "edginc/ISVBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** A type-specific interface to the MC Discount Factor State Variable */
class MCARLO_DLL SVDateOfDefault : public virtual ISVBase
{
public:
    virtual ~SVDateOfDefault() {}
    virtual DateTime getDateOfDefault() const = 0;
    virtual double   getRecoveryRateAtDefault() const = 0;
};
DECLARE(SVDateOfDefault);

/** A Generator of MC Date Of Default State Variables. Can return simulated date of default. */

class MCARLO_DLL SVGenDateOfDefault:  virtual public IElemStateVariableGen,
                                   public virtual VirtualDestructorBase
{
public:
 
    /** Constructor - no observation dates are necessary, 
        as the returned value of this SV is the date of default */
    SVGenDateOfDefault(const DateTime&     today,
                    ICDSParSpreadsConstSP  cdsParSpreadCurve,
                    const DateTime&        maxDate);

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface).
        The previous IStateVariableSP (may be null) should be passed in.
        The return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP                oldStateVar,
                                    IStateVariableGen::IStateGen*   pathGen) const;

    /** Returns a SVGenDateOfDefault state variable which then
        provides access to the path etc. This is the method that
        products should call to get an SVDateOfDefaultSP. */
    SVDateOfDefaultSP getSVDateOfDefault(IStateVariableGen::IStateGen* pathGen) const;

    /** Retrieve the CDS par spread curve associated with this
        SVGenDateOfDefault */
    //ICDSParSpreadsConstSP getCDSParSpreadCurve() const { return cdsParSpreadCurve; }

    /** Retrieve the credit name associated with this SV */
    string getCDSCurveName() const { return cdsParSpreadCurve->getName(); }

    const DateTime& getMaxDate() const  { return maxDate; }

    /** For use by Path Generators (past or future) that want to use determinstic rates */
    SVDateOfDefault* determinsticSV(bool doingPast) const;

    void attachSVGen(IElemStateVariableGenVisitor*) const;

private:
    /** Basic validation. Should be called by all constructors after
        population of all fields */
    void validate() {}
    class DeterminsticSV;
    DateTime                today;              // when to discount to
    DateTime                maxDate;            // the time horizon of our observation
    ICDSParSpreadsConstSP   cdsParSpreadCurve;
};

typedef smartPtr<SVGenDateOfDefault> SVGenDateOfDefaultSP;


DRLIB_END_NAMESPACE

#endif
