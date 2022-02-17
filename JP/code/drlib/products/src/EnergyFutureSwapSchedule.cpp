//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyFutureSwapSchedule.cpp
//
//   Description : Energy Swap instrument
//
//   Author      : Sean Chen
//
//   Date        : Aug. 29, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/EnergyFutureSwapSchedule.hpp"
#include "edginc/EnergyStreamSchedule.hpp"

#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE


void EnergyFutureSwapSchedule::validatePop2Object()
{
    static const string method("EnergyFutureSwapSchedule::validatePop2Object");

	EnergySwapScheduleBase::validatePop2Object();
	
	payFrequencyFixed = EnergyDWMYPeriod("1M");
	payFrequencyFloating = EnergyDWMYPeriod("1M");

	int aContractType = (scheduleType=="INDEX")?EnergyUnderlyer::INDEX:EnergyUnderlyer::OPTION;

	// Overwriting for contract period
    payFrequencyFixed.setToContractPeriod(aContractType);
	payFrequencyFloating.setToContractPeriod(aContractType);

	EnergyContractLabel aStartLabel = dealPeriod->getStartContractLabel();
	EnergyContractLabel aEndLabel = dealPeriod->getEndContractLabel();

    if (dealPeriod->isPrompt())
	{
		aStartLabel = energyUnderlyer.getSP()->calculateContractLabel(valueDate,aContractType);
		aEndLabel = aStartLabel;
	}
	else if (!(aStartLabel.isValid()) || !(aEndLabel.isValid()) )
	{
		// Deal period is wrong
		throw ModelException(method, "Deal period must be PROMPT or start/end labels");
	}

	EnergyContractLabel tmpLabel(aStartLabel);
    tmpLabel.addMonths(-1);
	startDate = energyUnderlyer.getSP()->expiryDate(tmpLabel,aContractType);
	endDate = energyUnderlyer.getSP()->expiryDate(aEndLabel,aContractType);

	string notionalSetting = "MATCHDELIVERY";


	// rollDay not used
	int rollDay = 0;

	fixedEnergyStreamScheduleSP = EnergyStreamScheduleFixedSP( 
		              new EnergyStreamScheduleFixed(
                          energyUnderlyer.getSP(), valueDate, startDate, endDate, settleDays,
						  payFrequencyFixed, avgDaysB4End, avgPeriod,
						  avgFrequency, rollDay, notionalSetting)
					  );

	fixedEnergyStreamScheduleSP->buildStreamSchedule();

	floatingEnergyStreamScheduleSP = EnergyStreamScheduleFloatingSP( 
		              new EnergyStreamScheduleFloating(
                          energyUnderlyer.getSP(), valueDate, startDate,endDate, settleDays,
						  payFrequencyFixed, avgDaysB4End, avgPeriod,
						  avgFrequency, rollDay, nearbyRel, nearbyAbsLabel, notionalSetting)
					  );

	floatingEnergyStreamScheduleSP->buildStreamSchedule();

}

EnergyFutureSwapSchedule::EnergyFutureSwapSchedule(): EnergySwapScheduleBase(TYPE),scheduleType("INDEX")
{    
}

EnergyFutureSwapSchedule::~EnergyFutureSwapSchedule()
{    
}

class EnergyFutureSwapScheduleHelper
{

public:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); 
        REGISTER(EnergyFutureSwapSchedule, clazz);
        SUPERCLASS(EnergySwapScheduleBase);
        EMPTY_SHELL_METHOD(defaultFutureSwapSchedule);

        FIELD(scheduleType,     "Schedule Type"); // INDEX or OPTION, rule to follow
		FIELD_MAKE_OPTIONAL(scheduleType);

    }

    static IObject* defaultFutureSwapSchedule()
    {
        return new EnergyFutureSwapSchedule();
    }
};
CClassConstSP const EnergyFutureSwapSchedule::TYPE = CClass::registerClassLoadMethod(
    "EnergyFutureSwapSchedule", typeid(EnergyFutureSwapSchedule), EnergyFutureSwapScheduleHelper::load);

bool  EnergyFutureSwapScheduleLoad() { return (EnergyFutureSwapSchedule::TYPE != 0);   }


DRLIB_END_NAMESPACE

