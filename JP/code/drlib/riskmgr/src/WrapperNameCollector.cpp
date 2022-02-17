//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : WrapperNameCollector.cpp
//
//   Description : Collects all wrapper names 
//
//   Author      : André Segger
//
//   Date        : 04 December 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/WrapperNameCollector.hpp"
#include "edginc/ObjectIteration.hpp"


DRLIB_BEGIN_NAMESPACE

void WrapperNameCollector::load(CClassSP& clazz) {
    REGISTER(WrapperNameCollector, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICollector);
}

CClassConstSP const WrapperNameCollector::TYPE = 
CClass::registerClassLoadMethod(
    "WrapperNameCollector", typeid(WrapperNameCollector), load);

WrapperNameCollector::WrapperNameCollector():CObject(TYPE) {
    wrapperNames = CStringArraySP(new CStringArray(0));
}

CStringArraySP WrapperNameCollector::getWrapperNames(IObjectSP obj)
{
    // class for handling call back
    class Action: public ObjectIteration::IActionConst{
    public:
        WrapperNameCollector collector;
        // called by ObjectIteration
        bool invoke(const ObjectIteration::State& state, IObjectConstSP obj) {
            obj->accept(&collector);
            return true;
        }
    };
    // create the action object
    Action action;
    // build an instance of the class that drives the recursion
    ObjectIteration iteration(CObject::TYPE);
    // go
    iteration.recurse(action, obj);

    return action.collector.wrapperNames;
}

/** Implementation of IAction interface */
void WrapperNameCollector::addName(const string& assetName){
    wrapperNames->push_back(assetName);
}


DRLIB_END_NAMESPACE
