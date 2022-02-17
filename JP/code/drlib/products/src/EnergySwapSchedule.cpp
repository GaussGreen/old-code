//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergySwapSchedule.cpp
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
#include "edginc/EnergySwapSchedule.hpp"

#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE


void EnergySwapSchedule::validatePop2Object()
{
    static const string method = "EnergySwapSchedule::validatePop2Object";

	EnergySwapScheduleBase::validatePop2Object();
	
	startDate = dealPeriod->getStartDate();
	endDate = dealPeriod->getEndDate();
	if ( startDate.empty() || endDate.empty() )
		throw ModelException(method, "Deal period must be a benchmark or start/end dates");

	if ( !payFrequencyFixed.getPeriod().size() )
		throw ModelException(method,"Payment Frequency is mandatory");
	if ( !payFrequencyFloating.getPeriod().size() )
		throw ModelException(method,"Payment Frequency is mandatory");


	fixedEnergyStreamScheduleSP = EnergyStreamScheduleFixedSP( 
		              new EnergyStreamScheduleFixed(
                          energyUnderlyer.getSP(), valueDate, startDate, endDate, settleDays,
						  payFrequencyFixed, avgDaysB4End, avgPeriod,
						  avgFrequency,rollDay, "MATCHCOUPON")
					  );

	fixedEnergyStreamScheduleSP->buildStreamSchedule();
    

	floatingEnergyStreamScheduleSP = EnergyStreamScheduleFloatingSP( 
		              new EnergyStreamScheduleFloating(
                          energyUnderlyer.getSP(), valueDate, startDate,endDate, settleDays,
						  payFrequencyFloating, avgDaysB4End, avgPeriod,
						  avgFrequency, rollDay, nearbyRel, nearbyAbsLabel,"MATCHCOUPON" )
					  );
	floatingEnergyStreamScheduleSP->buildStreamSchedule();
}


EnergySwapSchedule::EnergySwapSchedule(): EnergySwapScheduleBase(TYPE), rollDay(0)
{    
}

EnergySwapSchedule::~EnergySwapSchedule()
{    
}


class EnergySwapScheduleHelper
{

public:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); 
        REGISTER(EnergySwapSchedule, clazz);
        SUPERCLASS(EnergySwapScheduleBase);
        EMPTY_SHELL_METHOD(defaultSwapSchedule);

        FIELD(payFrequencyFixed,     "Fixed Payment Frequency"); // nD,M,W,Y
        FIELD(payFrequencyFloating,  "Floating Payment Frequency");

        FIELD(rollDay,         "Roll Day Overriding");
        FIELD_MAKE_OPTIONAL(rollDay);


    }

    static IObject* defaultSwapSchedule()
    {
        return new EnergySwapSchedule();
    }
};

CClassConstSP const EnergySwapSchedule::TYPE = CClass::registerClassLoadMethod(
    "EnergySwapSchedule", typeid(EnergySwapSchedule), EnergySwapScheduleHelper::load);

bool  EnergySwapScheduleLoad() { return (EnergySwapSchedule::TYPE != 0);   }

DRLIB_END_NAMESPACE

