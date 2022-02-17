//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ScenarioShift.cpp
//
//   Description : Defines a scenario shift to market data
//
//   Author      : Andrew J Swain
//
//   Date        : 24 April 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/RiskMgrInterface.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Control.hpp"

DRLIB_BEGIN_NAMESPACE

IScenarioShift::~IScenarioShift(){}

/** Invoked when class is 'loaded' */
void IScenarioShift::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IScenarioShift, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IScenarioShift::TYPE = CClass::registerInterfaceLoadMethod(
    "IScenarioShift", typeid(IScenarioShift), load);

DEFINE_TEMPLATE_TYPE(IScenarioShiftArray);


ScenarioShift::~ScenarioShift(){}

ScenarioShift::ScenarioShift(IPerturbationSP  sensCtrl, 
                             const string&    marketDataName) :
    CObject(TYPE), sensCtrl(sensCtrl), marketDataName(marketDataName) {
}

/** apply this scenario shift to the supplied object. The return
        value indicates if anything was actually shifted (true => yes) */
bool ScenarioShift::applyScenario(IObjectSP object) {
    TRACE_METHOD;
    static const string method = "ScenarioShift::applyScenario";
    try {
        if (!sensCtrl->applyBeforeGetMarket()) {
            TRACE("The IPerturbation wants to be applied now, after the getMarket phase");
            OutputNameSP outputName; // leave as null if string empty
            if (!marketDataName.empty()){
                outputName.reset(new OutputName(marketDataName));
            }
            return sensCtrl->findAndShift(object, outputName);
        }
        else {
            TRACE("(The IPerturbation doesn't want to be applied now, after the getMarket phase.)");
            return false;
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

// Nothing to do before market data is retrieved.
bool ScenarioShift::preapplyScenario(IObjectSP object){
    TRACE_METHOD;
    static const string method = "ScenarioShift::applyScenario";
    try {
        if (sensCtrl->applyBeforeGetMarket()) {
            TRACE("The IPerturbation wants to be applied now, before the getMarket phase");
            OutputNameSP outputName; // leave as null if string empty
            if (!marketDataName.empty()){
                outputName.reset(new OutputName(marketDataName));
            }
            return sensCtrl->findAndShift(object, outputName);
        }
        else {
            TRACE("(The IPerturbation doesn't want to be applied now, before the getMarket phase.)");
            return false;
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

class ScenarioShiftHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ScenarioShift, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IScenarioShift);
        EMPTY_SHELL_METHOD(defaultScenarioShift);
        FIELD(sensCtrl, "scenario to use in conjunction with marketDataName");
        FIELD(marketDataName, "market data name");
        FIELD_MAKE_OPTIONAL(marketDataName);
    }

    static IObject* defaultScenarioShift(){
        return new ScenarioShift();
    }
};

CClassConstSP const ScenarioShift::TYPE = CClass::registerClassLoadMethod(
    "ScenarioShift", typeid(ScenarioShift), ScenarioShiftHelper::load);

DEFINE_TEMPLATE_TYPE(ScenarioShiftArray);

// for reflection
ScenarioShift::ScenarioShift(): CObject(TYPE){}

/** For COGS: apply a native scenario to supplied RiskMgr object */
class ApplyScenario: public CObject,
                     public virtual ClientRunnable{
    RiskMgrInterfaceSP root; // market data must have already been retrieved
    ScenarioShiftSP    scenario;
public:
    ~ApplyScenario(){}

    static CClassConstSP const TYPE;

    virtual IObjectSP run(){
        // create TweakGroup (holds what we want to tweak)
        TweakGroupSP tweakGroup(new TweakGroup(root->inst, root->model));
        // apply shift
        scenario->applyScenario(tweakGroup);
        // return new RiskMgrInterface object
        return IObjectSP(new RiskMgrInterface(tweakGroup->getModelSP(),
                                              tweakGroup->getInstrumentSP(),
                                              root->ctrl,
                                              root->market));
    }
private:
    ApplyScenario(): CObject(TYPE){}

    static IObject* defaultConstructor(){
        return new ApplyScenario();
    }

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(ApplyScenario, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        IMPLEMENTS(ClientRunnable);
        FIELD(root, "Instrument, market, model data etc");
        FIELD(scenario, "scenario to apply");
    }
};

CClassConstSP const ApplyScenario::TYPE =
CClass::registerClassLoadMethod("ApplyScenario", typeid(ApplyScenario), load);

DRLIB_END_NAMESPACE

