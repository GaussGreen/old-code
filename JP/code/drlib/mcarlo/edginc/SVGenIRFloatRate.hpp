//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenIRFloatRate.hpp
//
//   Description : A Generator of MC IR Floater State Variables
//
//   Author      : Bruno O Melka
//
//   Date        : 23 Feb 2005
//
//
//----------------------------------------------------------------------------

#ifndef SVGenIRFloatRate_HPP
#define SVGenIRFloatRate_HPP

#include "edginc/StateVariableClient.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/SVGenIRSwap.hpp"


DRLIB_BEGIN_NAMESPACE
class MCARLO_DLL SVGenIRFloatRate: public virtual IStateVariableGen,
                public virtual IStateVariableClient,
                public virtual VirtualDestructorBase
{
public:

    /** Interface for the state variable that SVGenIRFloatRate produces. This is
        the type that products deal with in the payoff. The payoff obtains
        it by calling the getSVDiscFactor() method below. Note support here
        for only one yield curve per state variable (unlike SVGenSpot). May
        need to reconsider at some point */
    class MCARLO_DLL IStateVar: public virtual IStateVariable{
    public:
        virtual ~IStateVar();

		/** returns the Array of swap rates to apply in the leg. */
		virtual DoubleArraySP getYields() const = 0;

		/** returns the Array of year fractions to apply in the leg. */
		virtual DoubleArraySP getYearFracs() const = 0;
    };
    typedef smartPtr<IStateVar> IStateVarSP;

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector. Implementations typically call
        IStateVariableCollector::append */
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;

    /** Constructor for precomputed swap dates.
	If a discountCurve is supplied then the formula
    sum(FwdLibor*Z)/A is used otherwise (Z0-Zn)/A is used */
    SVGenIRFloatRate(
			YieldCurveConstSP		couponCurve,		//	for estimating coupons etc
			YieldCurveConstSP		discountCurve,		//	in domestic ccy
			const DateTimeArray&	refixDates,
			const DateTimeArray&	payDates,			//	coupon dates of swap
			string					rateInterval,		//	rate interval
			string					fixedPayInterval,   //	payment interval of fixed leg
			string					rateDCC,			//	day count convention of rate
			string					swapDCC,			//	day count convention of swap
			string					stubType,           //	type of front or back stub
			bool					stubAtEnd,          //	only matters if stubType != NONE
			string					accrueBadDayConv,
			string					payBadDayConv,
			HolidaySP				hols,
			bool					isCashSettled);		//	T=PV wrt ytm, F=PV wrt zcurve

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface). The
        previous IStateVariableSP (may be null) should be passed in.  The
        return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP              oldStateVar,
                                 IStateVariableGen::IStateGen*  pathGen) const;

    /** Returns a MC IRFloat state variable which then
        provides access to simulated values etc. This is the method that
        products should call to get an SVGenIRFloatRate::IStateVar. */
    IStateVarSP getIRFloatSV(IStateVariableGen::IStateGen* pathGen) const;

private:
	SVGenIRSwapArraySP		mcSwapRates;	// this is our StateVar generator
	DoubleArray			coupons;		// for a coupon rate of 1%. just to get the coverage
	class StateVar;
};

typedef smartPtr<SVGenIRFloatRate> SVGenIRFloatRateSP;

DRLIB_END_NAMESPACE

#endif
