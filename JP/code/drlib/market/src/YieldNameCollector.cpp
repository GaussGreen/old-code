//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : YieldNameCollector.cpp
//
//   Description : Collects all yield names 
//
//   Author      : Steohen
//
//   Date        : 23 Jan 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/YieldNameCollector.hpp"
#include "edginc/ObjectIteration.hpp"

DRLIB_BEGIN_NAMESPACE

void YieldNameCollector::load(CClassSP& clazz) {
    REGISTER(YieldNameCollector, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICollector);
}

CClassConstSP const YieldNameCollector::TYPE = CClass::registerClassLoadMethod(
    "YieldNameCollector", typeid(YieldNameCollector), load);

YieldNameCollector::YieldNameCollector(): 
    CObject(TYPE), yieldNames(new StringArray()){}

CStringArraySP YieldNameCollector::getYieldNames(IObjectSP obj)
{
    // class for handling call back
    class Action: public ObjectIteration::IActionConst{
    public:
        YieldNameCollector collector;
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

    return action.collector.yieldNames;
}

void YieldNameCollector::addName(const string& ycName){
    yieldNames->push_back(ycName);
}


DRLIB_END_NAMESPACE
