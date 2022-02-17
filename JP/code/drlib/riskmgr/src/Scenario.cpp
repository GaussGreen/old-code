//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Scenario.cpp
//
//   Description : Defines a collection of scenario shifts
//
//   Author      : Andrew J Swain
//
//   Date        : 24 April 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/ScenarioInterface.hpp"
#include "edginc/MultiScenarioInterface.hpp"
#include "edginc/SpreadSheetMode.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/Format.hpp"
#include "edginc/XMLWriter.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE


void Scenario::validatePop2Object()
{
    static const string method = "Scenario::validatePop2Object";
    for (int i = 0; i < shifts->size(); i++) {
            if (!(*shifts)[i]) {
                throw ModelException(method, "scenario shift " +
                                     Format::toString(i) + " is null");
            }
    }
}

/** Takes copy of sens array */
Scenario::Scenario(IScenarioShiftArrayConstSP shifts): CObject(TYPE){
    try{
        this->shifts = IScenarioShiftArraySP(shifts.clone());
    } catch (exception& e){
        throw ModelException(e, "Scenario::Scenario");
    }
}

/** Apply this scenario to supplied instrument/market */
void Scenario::apply(IModelSP                model, 
                     IInstrumentCollectionSP instruments) const { 
    static const string method = "Scenario::apply";
    try {
        // create TweakGroup (holds what we want to tweak)
        MultiTweakGroupSP tweakGroup(new MultiTweakGroup(instruments, model));
        for (int i = 0; i < shifts->size(); i++) {
            (*shifts)[i]->applyScenario(tweakGroup);
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Apply this scenario to supplied instrument/market, before market data is retrieved */
void Scenario::preapply(IModelSP                model, // (M)
                        IInstrumentCollectionSP instruments) const { // (M)
    static const string method = "Scenario::preapply";
    try {
        // create TweakGroup (holds what we want to tweak)
        MultiTweakGroupSP tweakGroup(new MultiTweakGroup(instruments, model));
        for (int i = 0; i < shifts->size(); i++) {
            if (!(*shifts)[i]) {
                throw ModelException(method, "scenario shift " +
                                     Format::toString(i) + " is null");
            }
            (*shifts)[i]->preapplyScenario(tweakGroup);
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Writes a ScenarioInterface object to file which consists of this and the
    supplied parameters. A ScenarioInterface object is capable of being run
    through the regression tester */
void Scenario::writeInputs(IModelSP      model,
                           CInstrumentSP instrument,
                           CControlSP    control,
                           MarketDataSP  market) 
{
    ScenarioInterface si(ScenarioSP::attachToRef(this),
                         model, instrument, control, market);
    XMLWriter xml(control->getFileName(true));
    si.write("SCENARIO", &xml); 
}

/** Writes a ScenarioInterface object to file which consists of this and the
    supplied parameters. A ScenarioInterface object is capable of being run
    through the regression tester */
void Scenario::writeInputs(IModelSP                model,
                           IInstrumentCollectionSP instruments,
                           CControlSP              control,
                           MarketDataSP            market) 
{
    MultiScenarioInterface msi(ScenarioSP::attachToRef(this),
                              model, instruments, control, market);
    XMLWriter xml(control->getFileName(true));
    msi.write("SCENARIO", &xml); 
}

/** Invoked when class is 'loaded' */
void Scenario::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Scenario, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultScenario);
    FIELD(shifts, "List of shifts to apply");
}

// for reflection
Scenario::Scenario(): CObject(TYPE){}

IObject* Scenario::defaultScenario(){
    return new Scenario();
}

CClassConstSP const Scenario::TYPE = CClass::registerClassLoadMethod(
    "Scenario", typeid(Scenario), load);

DRLIB_END_NAMESPACE

