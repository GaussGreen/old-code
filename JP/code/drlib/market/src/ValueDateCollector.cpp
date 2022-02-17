//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CValueDateCollector.cpp
//
//   Description : Value date collector class
//
//   Author      : Andre Segger
//
//   Date        : 30 Apr 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Equity.hpp"
#include "edginc/ObjectIteration.hpp"


DRLIB_BEGIN_NAMESPACE
/** Invoked when TestCollect is 'loaded' */
void CValueDateCollector::load(CClassSP& clazz){
    REGISTER(CValueDateCollector, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICollector);
}

CClassConstSP const CValueDateCollector::TYPE = CClass::registerClassLoadMethod(
    "CValueDateCollector", typeid(CValueDateCollector), load);

CValueDateCollector::CValueDateCollector(const DateTime& valueDate,
                                         const string&   source): CObject(TYPE),
    valueDate(valueDate), source(source) {}

void CValueDateCollector::valueDateValidate(const DateTime& valueDateToCompare,
                                            const string&   text)
{
    /* ensure option does not start before basket */
    if ( ! valueDate.equals(valueDateToCompare) ) {
        throw ModelException("ValueDateCollector::valueDateValidate", 
            text + " value date (" + valueDateToCompare.toString() + 
            ") is different from " + source + " value date (" +
            valueDate.toString() +")");
    }
}

void CValueDateCollector::validateAll(IObjectSP       obj,
                                      const DateTime& valueDate,
                                      const string&   text)
{
    // class for handling call back
    class Action: public ObjectIteration::IActionConst{
    public:
        CValueDateCollector* collector;
        // called by ObjectIteration
        bool invoke(const ObjectIteration::State& state, IObjectConstSP obj) {
            obj->accept(collector);
            return true;
        }
    };

    CValueDateCollector collector(valueDate, text);
    // create the action object
    Action action;
    action.collector = &collector;

    // build an instance of the class that drives the recursion
    ObjectIteration iteration1(CAsset::TYPE); 
    // go
    iteration1.recurse(action, obj);


    // build an instance of the class that drives the recursion
    ObjectIteration iteration2(Equity::TYPE); 
    // go
    iteration2.recurse(action, obj);

    // build an instance of the class that drives the recursion
    ObjectIteration iteration3(CVolBase::TYPE);
    // go
    iteration3.recurse(action, obj);
}

DRLIB_END_NAMESPACE

