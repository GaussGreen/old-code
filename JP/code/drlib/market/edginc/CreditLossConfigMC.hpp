//----------------------------------------------------------------------------
//
// Group       : CH Quantitative Research
//
// Filename    : CreditLossConfigMC.hpp
//
// Description : Representing interfaces that every ICreditLossConfig should implement to support MC pricing
//
// Date        : Aug 2006
//
//----------------------------------------------------------------------------


#ifndef CREDITLOSSCONFIGMC_HPP
#define CREDITLOSSCONFIGMC_HPP

#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/IProtectionProvider.hpp"
#include "edginc/IRebateCalculator.hpp"
#include "edginc/DateTimeLite.hpp"


DRLIB_BEGIN_NAMESPACE

// ################################################################################################################
/**
	Interface that a CreditLossConfig should implement to construct a SV generator that would be used in turn for MC pricing
*/
class MARKET_DLL ICreditLossConfigSVGenMC
{
public:
	/** constructor */
	ICreditLossConfigSVGenMC() {}

	/** virtual destructor */
	virtual ~ICreditLossConfigSVGenMC() {}

	/** method to create State Variable Generator */
	virtual ICreditLossConfigSVGenConstSP createSVGen(
		const DateTimeLiteVectorConstSP& timeline,
		const CIntConstSP& triggerDelay,
		const CIntConstSP& defaultToCalculationDelay,
		double temporaryLossAmount,
		const DateTime& lastTriggerDate,
		const AccrualPeriodArrayConstSP& accrualPeriods,
		const IBadDayAdjusterConstSP& bda,
		const IProtectionProviderConstSP& protect,
		const IRebateCalculatorConstSP& rebateCalc,
        const bool recoverNotional
		) const = 0;

private:
    ICreditLossConfigSVGenMC(const ICreditLossConfigSVGenMC&);
	ICreditLossConfigSVGenMC& operator=(const ICreditLossConfigSVGenMC&);
};

// ################################################################################################################
/**
	Interface that a CreditLossConfig should implement to construct a SV generator (of the indexed type) that would
	be used in turn for MC pricing.
*/
class MARKET_DLL ICreditLossConfigIndexedSVGenMC
{
public:
	/** constructor */
	ICreditLossConfigIndexedSVGenMC() {}

	/** virtual destructor */
	virtual ~ICreditLossConfigIndexedSVGenMC() {}

	/** method to create Indexed State Variable Generator */
	virtual ICreditLossConfigIndexedSVGenConstSP createIndexedSVGen(
		const DateTimeArrayConstSP& timeline,
		const CIntConstSP& triggerDelay,
		const CIntConstSP& defaultToCalculationDelay,
		double temporaryLossAmount,
		const DateTime& lastTriggerDate,
		const AccrualPeriodArrayConstSP& accrualPeriods,
		const IBadDayAdjusterConstSP& bda,
		const IProtectionProviderConstSP& protect,
		const IRebateCalculatorConstSP& rebateCalc,
        const bool recoverNotional
		) const = 0;

private:
    ICreditLossConfigIndexedSVGenMC(const ICreditLossConfigIndexedSVGenMC&);
    ICreditLossConfigIndexedSVGenMC& operator=(const ICreditLossConfigIndexedSVGenMC&);
};

// ################################################################################################################

DRLIB_END_NAMESPACE

#endif







