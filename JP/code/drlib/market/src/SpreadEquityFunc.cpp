//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpreadEquityFunc.cpp
//
//   Description : 
//
//   Author      : Qing Hou
//
//   Date        : 18 Dec 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/SpreadEquityFunc.hpp"
#include "edginc/TimeMetric.hpp"

DRLIB_BEGIN_NAMESPACE

SpreadEquityFunc::SpreadEquityFunc(const CClassConstSP& clazz) 
: CObject(clazz)
{}

void SpreadEquityFunc::setTimeMetric(const TimeMetric *metric) 
{
	this->metric = TimeMetricSP(copy(metric));
}

void SpreadEquityFunc::load(CClassSP& clazz)
{
	clazz->setPrivate(); // make invisible to EAS/spreadsheet
	REGISTER(SpreadEquityFunc, clazz);
	SUPERCLASS(CObject);
	FIELD(metric,"");
}

CClassConstSP const SpreadEquityFunc::TYPE = CClass::registerClassLoadMethod(
              "SpreadEquityFunc", typeid(SpreadEquityFunc), SpreadEquityFunc::load);

class SpreadEquityFuncAddin: public CObject{
    static CClassConstSP const TYPE;

    SpreadEquityFuncSP      spreadFunc;
    DateTime				startDate;   
    DateTime				endDate;
	double					strike;

    static double getSpread(SpreadEquityFuncAddin* params){
        static const string routine = "SpreadEquityFuncAddin::getSpread";
        try {

            return params->spreadFunc->getSpreadCC(params->strike,	params->startDate, params->endDate);

        } catch (exception& e){
            throw ModelException(e, routine);
        }
	}

    /** for reflection */
    SpreadEquityFuncAddin():  CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(SpreadEquityFuncAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSpreadEquityFuncAddin);
        FIELD(spreadFunc,		"spread function");
        FIELD(startDate,	"start date");
        FIELD(endDate,	"end date");
        FIELD(strike,	"equity price");

        Addin::registerInstanceDoubleMethod("GET_SPREAD",
                                         Addin::MARKET,
                                         "get equity price dependent spread",
                                         TYPE,
                                         (Addin::DoubleMethod*)getSpread);
    }

    static IObject* defaultSpreadEquityFuncAddin(){
        return new SpreadEquityFuncAddin();
    }   
};

CClassConstSP const SpreadEquityFuncAddin::TYPE = CClass::registerClassLoadMethod(
    "SpreadEquityFuncAddin", typeid(SpreadEquityFuncAddin), load);

DRLIB_END_NAMESPACE
