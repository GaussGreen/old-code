//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyTermPeriod.cpp
//
//   Description : Defines a deal period for energy. 
//
//   Author      : Sean Chen
//
//   Date        : 26 Sept. 2005
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Format.hpp"

#include "edginc/EnergyTermPeriod.hpp"


DRLIB_BEGIN_NAMESPACE


EnergyTermPeriod::EnergyTermPeriod(const EnergyTermPeriod& thePeriod):
    CObject(TYPE), startDate(thePeriod.startDate), endDate(thePeriod.endDate),
	startContract(thePeriod.startContract), endContract(thePeriod.endContract), prompt(thePeriod.prompt)
{}

void EnergyTermPeriod::validatePop2Object()
{
    static const string routine = "EnergyTermPeriod::validatePop2Object";   
}

EnergyTermPeriod::EnergyTermPeriod(CClassConstSP clazz) : CObject(clazz), prompt(false)
{
}

EnergyTermPeriod::~EnergyTermPeriod()
{
}


DateTime EnergyTermPeriod::getStartDate() const
{
    return startDate;
}

DateTime EnergyTermPeriod::getEndDate() const
{
    return endDate;
}

EnergyContractLabel EnergyTermPeriod::getStartContractLabel() const
{
	return startContract;
}

EnergyContractLabel EnergyTermPeriod::getEndContractLabel() const
{
	return endContract;
}

bool EnergyTermPeriod::isPrompt() const
{
	return prompt;
}


/* for reflection */
class EnergyTermPeriodHelper
{

public:

   /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EnergyTermPeriod, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultEnergyTermPeriod);
        FIELD(startDate, "Term Start Date");
		FIELD_MAKE_TRANSIENT(startDate);
		FIELD(endDate, "Term End Date");
		FIELD_MAKE_TRANSIENT(endDate);
		FIELD(startContract, "Start Contract");
		FIELD_MAKE_TRANSIENT(startContract);
		FIELD(endContract, "End Contract");
		FIELD_MAKE_TRANSIENT(endContract);
    }
    
    static IObject* defaultEnergyTermPeriod(){
        return new EnergyTermPeriod();
    }
};

CClassConstSP const EnergyTermPeriod::TYPE = CClass::registerClassLoadMethod(
    "EnergyTermPeriod", typeid(EnergyTermPeriod), EnergyTermPeriodHelper::load);


DRLIB_END_NAMESPACE
