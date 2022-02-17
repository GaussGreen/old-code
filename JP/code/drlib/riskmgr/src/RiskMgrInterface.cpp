//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RiskMgrInterface.cpp
//
//   Description : Defines interface to RiskMgr
//
//   Author      : Andrew J Swain
//
//   Date        : 19 January 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RiskMgrInterface.hpp"
#include "edginc/MultiRiskMgrInterface.hpp"
#include "edginc/Addin.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/Scenario.hpp"

DRLIB_BEGIN_NAMESPACE

RiskMgrInterface::RiskMgrInterface(IModelSP      model, 
                                   CInstrumentSP inst,
                                   CControlSP    ctrl,
                                   CMarketDataSP market) :
    CObject(TYPE), model(model), inst(inst), ctrl(ctrl), market(market) {}

// for reflection
RiskMgrInterface::RiskMgrInterface(): CObject(TYPE){}

/** Runs 'regression test' for model, inst and control in this class */
IObjectSP RiskMgrInterface::runTest() const{
    // don't want to write out an input file
    ctrl->switchOffWriteToFile();
    return go();
}

// EdrAction 
IObjectSP RiskMgrInterface::run(){
    return go();
}

IObjectSP RiskMgrInterface::go() const{
    CMarketDataSP marketData(market);
    if (!marketData){
        // review in conjunction with avoiding writing all market data out
        marketData.reset(new MarketData());
    }
    CControlSP control = ctrl->applyCommandLineOptions();
    return model->go(inst, ScenarioSP(), control, marketData);
}

class RiskMgrInterfaceHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RiskMgrInterface, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultInterface);
        FIELD(model, "model");
        FIELD(inst, "instrument");
        FIELD(ctrl, "ctrl");
        FIELD(market, "market");
        FIELD_MAKE_OPTIONAL(market);
        // registration for addin function
        Addin::registerObjectMethod("RISK_MGR",
                                    Addin::RISK,
                                    "Calculates price and greeks",
                                    true,
                                    Addin::returnHandle,
                                    &RiskMgrInterface::run);
    }

    static IObject* defaultInterface(){
        return new RiskMgrInterface();
    }
};

CClassConstSP const RiskMgrInterface::TYPE = CClass::registerClassLoadMethod(
    "RiskMgrInterface", typeid(RiskMgrInterface), RiskMgrInterfaceHelper::load);

DRLIB_END_NAMESPACE

