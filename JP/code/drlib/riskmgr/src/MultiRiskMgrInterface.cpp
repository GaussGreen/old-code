/**
 * @file MultiRiskMgrInterface.cpp
 */

#include "edginc/config.hpp"
#include "edginc/MultiRiskMgrInterface.hpp"
#include "edginc/Addin.hpp"
#include "edginc/RiskMgr.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/Scenario.hpp"

DRLIB_BEGIN_NAMESPACE

MultiRiskMgrInterface::MultiRiskMgrInterface(IModelSP model, 
                                             IInstrumentCollectionSP insts,
                                             CControlSP ctrl,
                                             CMarketDataSP market) :
    CObject(TYPE), model(model), insts(insts), ctrl(ctrl), market(market) {
    // empty
}

// for reflection
MultiRiskMgrInterface::MultiRiskMgrInterface(): CObject(TYPE){}

/** Runs 'regression test' for model, insts and control in this class */
IObjectSP MultiRiskMgrInterface::runTest() const{
    // don't want to write out an input file
    ctrl->switchOffWriteToFile();
    return go();
}

// EdrAction 
IObjectSP MultiRiskMgrInterface::run(){
    return go();
}

IObjectSP MultiRiskMgrInterface::go() const{
    CMarketDataSP marketData(market);
    if (!marketData){
        // review in conjunction with avoiding writing all market data out
        marketData.reset(new MarketData());
    }
    CControlSP control = ctrl->applyCommandLineOptions();
    return model->go(insts, ScenarioSP(), control, market);
}

IObject *MultiRiskMgrInterface::defaultMultiRiskMgrInterface() {
    return new MultiRiskMgrInterface();
}

void MultiRiskMgrInterface::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(MultiRiskMgrInterface, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultMultiRiskMgrInterface);
    IMPLEMENTS(ClientRunnable);
    FIELD(model, "model");
    FIELD(insts, "instruments");
    FIELD(ctrl, "ctrl");
    FIELD(market, "market");
    FIELD_MAKE_OPTIONAL(market);

    // registration for addin function
    Addin::registerObjectMethod("MULTI_RISK_MGR",
                                Addin::RISK,
                                "Calculates price and greeks",
                                true,
                                Addin::returnHandle,
                                &MultiRiskMgrInterface::run);
}

CClassConstSP const MultiRiskMgrInterface::TYPE =
CClass::registerClassLoadMethod(
    "MultiRiskMgrInterface", typeid(MultiRiskMgrInterface), load);

bool MultiRiskMgrInterfaceLinkIn() {
    return MultiRiskMgrInterface::TYPE != NULL;
}

DRLIB_END_NAMESPACE
