//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : ExposureReporter.cpp
//
//   Description : Given an imnt/model/market, return a set of Results
//                 where a non-zero entry in each element indicates exposure 
//                 exists - essentially figures out where exposure is without
//                 a full blown pricing
//
//   Author      : Andrew J Swain
//
//   Date        : 19 December 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/Model.hpp"
#include "edginc/ModelFilter.hpp"
#include "edginc/Scenario.hpp"
#include "edginc/RiskMgr.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE

class ExposureReporter: public CObject,
                        virtual public ClientRunnable {
public:
    static CClassConstSP const TYPE;
    
    // EdrAction 
    IObjectSP run() {
        // create a highlighter from the "real" model
        IModelSP highlighter(model->exposureHighlighter());

        // disable write to file
        CControlSP control(ctrl.clone());
        control->switchOffWriteToFile();

        return highlighter->go(inst, ScenarioSP(), control, market);
    }

private:
    /** for addin */
    static IObjectSP addinFunc(ExposureReporter* params) {
        return params->run();
    }

    ExposureReporter():CObject(TYPE), model(0), inst(0), ctrl(0), market(0){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ExposureReporter, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultInterface);
        FIELD(model, "model");
        FIELD(inst, "instrument");
        FIELD(ctrl, "ctrl");
        FIELD(market, "market");
        // registration for addin function
        Addin::registerClassObjectMethod("EXPOSURE_REPORTER",
                                         Addin::RISK,
                                         "Highlights sensitive regions of market data",
                                         ExposureReporter::TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)addinFunc);
        
    }

    static IObject* defaultInterface(){
        return new ExposureReporter();
    }

    // fields
    IModelSP      model;
    CInstrumentSP inst;
    CControlSP    ctrl;
    CMarketDataSP market;
};

bool ExposureReporterLinkIn() {
    return ExposureReporter::TYPE != NULL;
}

CClassConstSP const ExposureReporter::TYPE = CClass::registerClassLoadMethod(
    "ExposureReporter", typeid(ExposureReporter), ExposureReporter::load);
   



DRLIB_END_NAMESPACE
