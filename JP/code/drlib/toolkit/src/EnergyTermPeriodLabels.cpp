//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyTermPeriodLabels.cpp
//
//   Description : Defines a deal period with a benchmark for energy. 
//
//   Author      : Sean Chen
//
//   Date        : 26 Sept. 2005
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Format.hpp"

#include "edginc/EnergyTermPeriodLabels.hpp"

DRLIB_BEGIN_NAMESPACE


void EnergyTermPeriodLabels::validatePop2Object()
{
    static const string routine = "EnergyTermPeriodLabels::validatePop2Object"; 
	
	try
	{
	    startContract = EnergyContractLabel(startLabel);
	    endContract   = EnergyContractLabel(endLabel);
	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}
}

EnergyTermPeriodLabels::EnergyTermPeriodLabels() : EnergyTermPeriod(TYPE)
{
}

EnergyTermPeriodLabels::~EnergyTermPeriodLabels()
{
}

/* for reflection */
class EnergyTermPeriodLabelsHelper
{

public:

   /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EnergyTermPeriodLabels, clazz);
        SUPERCLASS(EnergyTermPeriod);
        EMPTY_SHELL_METHOD(defaultEnergyTermPeriodLabels);
        FIELD(startLabel, "Start Label");
		FIELD(endLabel, "End Label");

    }
    
    static IObject* defaultEnergyTermPeriodLabels(){
        return new EnergyTermPeriodLabels();
    }
};

CClassConstSP const EnergyTermPeriodLabels::TYPE = CClass::registerClassLoadMethod(
    "EnergyTermPeriodLabels", typeid(EnergyTermPeriodLabels), EnergyTermPeriodLabelsHelper::load);


DRLIB_END_NAMESPACE
