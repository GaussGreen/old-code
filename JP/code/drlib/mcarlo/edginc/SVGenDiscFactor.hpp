//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenDiscFactor.hpp
//
//   Description : A Generator of MC Discount Factor State Variables
//
//   Date        : 20 May 2004
//
//
//----------------------------------------------------------------------------

#ifndef SVGenDiscFactor_HPP
#define SVGenDiscFactor_HPP

#include "edginc/YieldCurve.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/ISVBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** A type-specific interface to the MC Discount Factor State Variable */
class MCARLO_DLL SVDiscFactor : public virtual ISVBase
{
public:
    virtual ~SVDiscFactor() {}
    virtual double getSpotDF() const = 0;
    virtual double firstDF() const = 0;
    virtual double getDF(int i) const = 0;
};
DECLARE(SVDiscFactor);

/** A Generator of MC Discount Factor State Variables. Can return simulated
    discount factors between today and a set of dates. The set of dates can
    be specified in a variety of manners */
class MCARLO_DLL SVGenDiscFactor: virtual public IElemStateVariableGen,
                               public virtual VirtualDestructorBase
{
public:

    /** Note: The SVGenDiscFactor::IStateVar has been obsoleted by the
        class MCARLO_DLL SVQmcDiscFactor.   */

    /** Constructor - from  an array of dates. For computing discount
        factors between today and each date in dates */
    SVGenDiscFactor(const DateTime&      today,
                 YieldCurveConstSP    yieldCurve,
                 const DateTimeArray& dates);

    /** Constructor - from  an array of dates and a settlement */
    SVGenDiscFactor(const DateTime&        today,
                 YieldCurveConstSP      yieldCurve,
                 InstrumentSettlementSP instSettle,
                 const DateTimeArray&   dates);

    /** Constructor - from a single date */
    SVGenDiscFactor(const DateTime&     today,
                 YieldCurveConstSP   yieldCurve,
                 const DateTime&     maturityDate);

    /** Constructor - from a single date and a settlement */
    SVGenDiscFactor(const DateTime&        today,
                 YieldCurveConstSP      yieldCurve,
                 InstrumentSettlementSP instSettle,
                 const DateTime&        maturityDate);

    /** Constructor - from a single date and a settlement. Does support
        physical settlement */
    SVGenDiscFactor(const DateTime&        today,
                 YieldCurveConstSP      yieldCurve,
                 InstrumentSettlementSP instSettle,
                 const DateTime&        maturityDate,
                 CAssetSP               asset);

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface). The
        previous IStateVariableSP (may be null) should be passed in.  The
        return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP             oldStateVar,
                                 IStateVariableGen::IStateGen* pathGen) const;

    /** Returns a MC Discount Factor state variable which then
        provides access to the path etc. This is the method that
        products should call to get an SVQmcDiscFactor. */
    SVDiscFactorSP getSVDiscFactor(IStateVariableGen::IStateGen* pathGen) const;

    /** Retrieve the yield curve associated with this SVGenDiscFactor */
    YieldCurveConstSP getYieldCurve() const;
    /** Retrieve the dates for which a discount factor is required */
    const DateTimeArray&   getDates() const;
    /** For use by Path Generators (past or future) that want to use
        determinstic rates */
    SVDiscFactor* determinsticSV(bool doingPast) const;
    void attachSVGen(IElemStateVariableGenVisitor*) const;
private:
    void initDatesFromInstSettlement(
        InstrumentSettlementSP instSettle,
        const DateTimeArray&   theDates,
        CAsset*                asset);

    /** Basic validation. Should be called by all
        constructors after population of all fields */
    void validate();

    class DeterminsticSV;
    DateTime          today; // when to discount to
    bool              isMargin;
    YieldCurveConstSP yieldCurve;
    DateTimeArray     originalDates;    //!< Dates at which discount factors are requested
    DateTimeArray     adjustedDates;    //!< Settlement adjusted discount factor dates
};

typedef smartPtr<SVGenDiscFactor> SVGenDiscFactorSP;
typedef smartConstPtr<SVGenDiscFactor> SVGenDiscFactorConstSP;
typedef vector<const SVGenDiscFactor*> SVGenDiscFactorArray;
//DECLARE(SVGenDiscFactor);

DRLIB_END_NAMESPACE

#endif // SVGenDiscFactor
