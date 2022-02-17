//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenIRFloatRate.cpp
//
//   Description : A Generator of MC IR Floater State Variables
//
//   Author      : Bruno O Melka
//
//   Date        : 23 Feb 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenIRFloatRate.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/DayCountConventionFactory.hpp"


DRLIB_BEGIN_NAMESPACE
SVGenIRFloatRate::IStateVar::~IStateVar(){}

/** Appends 'true' (ie non derived) state variable generators
    required to the supplied collector. Implementations typically call
    IStateVariableCollector::append */
void SVGenIRFloatRate::collectStateVars(IStateVariableCollectorSP svCollector) const {
	for (size_t j = 0; j < mcSwapRates->size(); j++)
	{
		svCollector->append((*mcSwapRates)[j].get());
	}
}

/** Constructor for precomputed refix dates.
	If a discountCurve is supplied then the formula
    sum(FwdLibor*Z)/A is used otherwise (Z0-Zn)/A is used */
SVGenIRFloatRate::SVGenIRFloatRate( 
			YieldCurveConstSP		couponCurve,		//	for estimating coupons etc
			YieldCurveConstSP		discountCurve,		//	in domestic ccy
			const DateTimeArray&	refixDates,
			const DateTimeArray&	payDates,			//	coupon dates of swap
			string					rateInterval,		//	rate interval
			string					fixedPayInterval,   //	payment interval of fixed leg
			string					rateDCC,	        //	day count convention of Rate
			string					swapDCC,			//	day count convention of swap
			string					stubType,           //	type of front or back stub
			bool					stubAtEnd,          //	only matters if stubType != NONE
			string					accrueBadDayConv,    
			string					payBadDayConv,       
			HolidaySP				hols,
			bool					isCashSettled):		//	T=PV wrt ytm, F=PV wrt zcurve
			coupons(refixDates.size()),
			mcSwapRates(new SVGenIRSwapArray(refixDates.size())) {

	try {
		
		BadDayConventionSP		accrueBDC(BadDayConventionFactory::make(accrueBadDayConv));
		DayCountConventionSP	dccSwap(DayCountConventionFactory::make(swapDCC));
		MaturityPeriod			interval(rateInterval);


		for (int i = 0; i < refixDates.size(); i++){
			(*mcSwapRates)[i] = SVGenIRSwapSP(new SVGenIRSwap(couponCurve, 
															discountCurve, 
															refixDates[i],
															accrueBDC->adjust(interval.toDate(refixDates[i]), hols.get()),
															fixedPayInterval,
															rateDCC,
															stubType,
															stubAtEnd,
															accrueBadDayConv,
															payBadDayConv,
															hols,
															isCashSettled));
			coupons[i] = dccSwap->years(refixDates[i],payDates[i]); // rate = 1.0
		}
	}
	catch (exception& e){ 
        throw ModelException(e, "SVGenIRFloatRate::SVGenIRFloatRate","Failed to build swap rates");
    }
}


/** Create the corresponding State Variable for this State
    Variable Generator (from IStateVariableGen interface). The
    previous IStateVariableSP (may be null) should be passed in.  The
    return object may or may not be the same as oldStateVar. */
IStateVariableSP SVGenIRFloatRate::create(IStateVariableSP              oldStateVar,
											IStateVariableGen::IStateGen*  pathGen) const {
    return getIRFloatSV(pathGen);
}

class SVGenIRFloatRate::StateVar: public IStateVar{
public:
    virtual ~StateVar(){}

	/** Returns true if this state variable is being used to 'simulate'
        the past. This is a useful method for users of state variables - 
        need to see how hard it is to implement */
    bool doingPast() const{
        return false; // todo swapRates->doingPast();
    }
	
	/** returns the Array of swap rates to apply in the leg. */
	virtual DoubleArraySP getYields() const {
		double parYield, annuity;
		DoubleArraySP legYields(new DoubleArray(swapRates->size()));
		for (size_t j = 0; j < swapRates->size(); j++) {
			(*swapRates)[j].get()->parYield(parYield, annuity);
			(*legYields)[j] = parYield;
		}
		return(legYields);
	}

	/** returns the Array of year fractions to apply in the leg. */
	virtual DoubleArraySP getYearFracs() const {
		DoubleArraySP legYearFracs(new DoubleArray(swapRates->size()));
		for (size_t j = 0; j < swapRates->size(); j++) {
			(*legYearFracs)[j] = yearFracs[j];
		}
		return(legYearFracs);
	}

    // constructor
    StateVar(SVGenIRSwap::StateVarArraySP swapRates,
				const DoubleArray& coupons):
				swapRates(swapRates),
				yearFracs(coupons.size()) {
		for (size_t i = 0; i < swapRates->size(); i++) {
			yearFracs[i] = coupons[i];
		}
	}

protected:
	SVGenIRSwap::StateVarArraySP	swapRates;
	vector<double>				yearFracs;

};

/** Returns a MC IRSwap state variable which then
    provides access SVQmcDiscFactorSP discFactoto simulated values etc. This is the method that
    products should call to get an SVGenIRFloatRate::IStateVar. */
SVGenIRFloatRate::IStateVarSP SVGenIRFloatRate::getIRFloatSV(IStateVariableGen::IStateGen* pathGen) const{

	SVGenIRSwap::StateVarArraySP swapRates(new SVGenIRSwap::StateVarArray(mcSwapRates->size()));
	for (size_t i =0; i < mcSwapRates->size(); i++)
	{
		(*swapRates)[i] = SVGenIRSwap::IStateVarSP(((*mcSwapRates)[i])->getIRSwapSV(pathGen));
	}
	return SVGenIRFloatRate::IStateVarSP(new StateVar(swapRates, coupons)); 
}

DRLIB_END_NAMESPACE
