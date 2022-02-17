//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyTermPeriodDates.cpp
//
//   Description : Defines a deal period with start/end dates for energy. 
//
//   Author      : Sean Chen
//
//   Date        : 26 Sept. 2005
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/EnergyTermPeriodDates.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE


void EnergyTermPeriodDates::validatePop2Object()
{
    static const string routine = "EnergyTermPeriodDates::validatePop2Object"; 
	
	if (dateStart>= dateEnd)
		throw ModelException(routine, "Start date must be before End Date");

	startDate = dateStart;
	endDate = dateEnd;
}

EnergyTermPeriodDates::EnergyTermPeriodDates() : EnergyTermPeriod(TYPE)
{
}

EnergyTermPeriodDates::~EnergyTermPeriodDates()
{
}


/* for reflection */
class EnergyTermPeriodDatesHelper
{

public:

   /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EnergyTermPeriodDates, clazz);
        SUPERCLASS(EnergyTermPeriod);
        EMPTY_SHELL_METHOD(defaultEnergyTermPeriodDates);
        FIELD(dateStart, "Term Start Date");
		FIELD(dateEnd, "Term End Date");
		
    }
    
    static IObject* defaultEnergyTermPeriodDates(){
        return new EnergyTermPeriodDates();
    }
};

CClassConstSP const EnergyTermPeriodDates::TYPE = CClass::registerClassLoadMethod(
    "EnergyTermPeriodDates", typeid(EnergyTermPeriodDates), EnergyTermPeriodDatesHelper::load);


DRLIB_END_NAMESPACE
