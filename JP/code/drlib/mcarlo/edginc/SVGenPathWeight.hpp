//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenPathWeight.hpp
//
//   Description : A Generator of MC Path Weight State Variables
//                 (SV as they can develop to become instrument specific)
//
//   Date        : 22 Nov 2006
//
//
//----------------------------------------------------------------------------

#ifndef SVGenPathWeight_HPP
#define SVGenPathWeight_HPP

#include "edginc/ISVBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** A type-specific interface to the MC Discount Factor State Variable */
class MCARLO_DLL SVPathWeight : public virtual ISVBase
{
public:
    virtual double getWeight(int i) const = 0;
};
DECLARE(SVPathWeight);

/** A Generator of MC Discount Factor State Variables. Can return simulated
    discount factors between today and a set of dates. The set of dates can
    be specified in a variety of manners */
class MCARLO_DLL SVGenPathWeight: virtual public IElemStateVariableGen,
                                  public virtual VirtualDestructorBase
{
public:

    /** Note: The SVGenDiscFactor::IStateVar has been obsoleted by the
        class MCARLO_DLL SVQmcDiscFactor.   */

    /** Constructor - from  an array of dates. For computing discount
        factors between today and each date in dates */
    SVGenPathWeight(const DateTimeArray& dates);

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface). The
        previous IStateVariableSP (may be null) should be passed in.  The
        return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP             oldStateVar,
                                 IStateVariableGen::IStateGen* pathGen) const;

    /** Returns a MC Discount Factor state variable which then
        provides access to the path etc. This is the method that
        products should call to get an SVQmcDiscFactor. */
    SVPathWeightSP getSVPathWeight(IStateVariableGen::IStateGen* pathGen) const;

    /** Retrieve the dates for which a discount factor is required */
    const DateTimeArray& getDates() const
        { return dates; }

    /** For use by Path Generators (past or future) that want to use
        deterministic rates */
    SVPathWeight* determinsticSV(bool doingPast) const;
    void attachSVGen(IElemStateVariableGenVisitor*) const;
private:

    class DeterminsticSV;
    DateTimeArray     dates;    //!< Dates at which discount factors are requested
};

DECLARE(SVGenPathWeight);

DRLIB_END_NAMESPACE

#endif // SVGenPathWeight
