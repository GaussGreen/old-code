//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DRWrapperInterface.cpp
//
//   Description : Defines DR Wrapper interface to RiskMgr
//
//   Author      : Andrew J Swain
//
//   Date        : 20 April 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DRWrapperInterface.hpp"
#include "edginc/DRWrapper.hpp"
#include "edginc/Addin.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE

DRWrapperInterface::DRWrapperInterface(ScenarioSP    scenario,
                                       IObjectSP     model, 
                                       IObjectSP     inst,
                                       IObjectSP     ctrl,
                                       CMarketDataSP market) :
    CObject(TYPE), scenario(scenario), model(model), 
    inst(inst), ctrl(ctrl), market(market) {
    // empty
}

// for reflection
DRWrapperInterface::DRWrapperInterface() : 
    CObject(TYPE), scenario(0), model(0), inst(0), ctrl(0), market(0){
    // empty
}

/** Runs 'regression test' for model, inst and control in this class */
IObjectSP DRWrapperInterface::runTest() const{
    static const string method = "DRWrapperInterface::runTest";
    IObjectSP results;
    try {
        CMarketDataSP marketData(market);
        if (!marketData){
            // review in conjunction with avoiding writing all market data out
            marketData = CMarketDataSP(new MarketData());
        }
        // if any of the 3 Objects are actually wrappers, convert them
        // to real objects

        IModelSP      modelSP;
        CInstrumentSP instSP;
        CControlSP    ctrlSP;
        DRWrapperSP   wrapperSP;

        if (DRWrapper::TYPE->isInstance(model)) {
            wrapperSP = DRWrapperSP::dynamicCast(model);        
            modelSP   = IModelSP::dynamicCast(wrapperSP->pop2Object());
        }
        else {
            modelSP = IModelSP::dynamicCast(model);
        }

        if (!modelSP) {
            throw ModelException(method, "couldn't interpret model parameter");
        }

        if (DRWrapper::TYPE->isInstance(inst)) {
            wrapperSP = DRWrapperSP::dynamicCast(inst);   
            wrapperSP->validateExeName();
        
            instSP    = CInstrumentSP::dynamicCast(wrapperSP->pop2Object());
        }
        else {
            instSP = CInstrumentSP::dynamicCast(inst);
        }

        if (!instSP) {
            throw ModelException(method, 
                                 "couldn't interpret instrument parameter");
        }

        if (DRWrapper::TYPE->isInstance(ctrl)) {
            wrapperSP = DRWrapperSP::dynamicCast(ctrl);        
            ctrlSP    = CControlSP::dynamicCast(wrapperSP->pop2Object());
        }
        else {
            ctrlSP = CControlSP::dynamicCast(ctrl);
        }

        if (!ctrlSP) {
            throw ModelException(method, "couldn't interpret ctrl parameter");
        }
        
        results = modelSP->go(instSP, scenario, ctrlSP, marketData);
    }
    catch (exception& e){
        throw ModelException(&e, method);
    }

    return results;
}

/** addin function wrapper for DRWrapper */
static IObjectSP addinDRWrapper(DRWrapperInterface* params){
    // copy instrument before we run - avoids problems of whether market
    // data has been obtained or not
    params->inst = IObjectSP(params->inst.clone());
    return params->runTest();
}


class DRWrapperInterfaceHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DRWrapperInterface, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultInterface);
        FIELD(scenario, "scenario (optional)");
        FIELD_MAKE_OPTIONAL(scenario);  
        FIELD(model, "model");
        FIELD(inst, "instrument");
        FIELD(ctrl, "ctrl");
        FIELD(market, "market");
        FIELD_MAKE_OPTIONAL(market);  

        // registration for addin function
        Addin::registerClassObjectMethod("RUN_DR_WRAPPER",
                                         Addin::RISK,
                                         "Calculates price and greeks",
                                         DRWrapperInterface::TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)addinDRWrapper);
              
    }

    static IObject* defaultInterface(){
        return new DRWrapperInterface();
    }
};

CClassConstSP const DRWrapperInterface::TYPE = CClass::registerClassLoadMethod(
    "DRWrapperInterface", typeid(DRWrapperInterface), DRWrapperInterfaceHelper::load);
   

DRLIB_END_NAMESPACE

