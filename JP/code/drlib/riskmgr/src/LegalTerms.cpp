//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LegalTerms.cpp
//
//   Description : Transform object into one which matches term sheet
//                 rather than one relevant for pricing/greeks
//                 e.g. use economic barrier instead of risk barrier
//
//   Date        : May 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/LegalTerms.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/Results.hpp"
#include "edginc/Model.hpp"

DRLIB_BEGIN_NAMESPACE
LegalTerms::Shift::~Shift(){} // empty

/** for reflection */
LegalTerms::LegalTerms(): CObject(TYPE){}


//// LEGAL_TERMS_FV output request is not handled by the instrument
//// but at the top level through the scenario. This calculator class
//// sits inside the OutputRequestHelper
class LegalTerms::Calculator : public OutputRequestCalculator{
public:
    //// implementation to handle top level LEGAL_TERMS_FV out put request
    void calculate(OutputRequest*     request, 
                   const IModel*      model, 
                   const CInstrument* instrument,
                   Results*           results) {
        // build the scenario on the fly
        IScenarioShiftSP legalTermsScenario(new LegalTerms());
        // preserve original instrument
        CInstrumentSP copyInst(copy(instrument));
        // we don't want to calculate the FV again if there are no legal terms
        if (legalTermsScenario->applyScenario(copyInst)){
            IModelSP copyModel(copy(model));
            CControlSP ctrl(new Control(SensitivityArrayConstSP(   ), 
                                        OutputRequestArrayConstSP(   ),false,""));
            ResultsSP tempResults(new Results());
            copyModel->Price(copyInst.get(), ctrl.get(), tempResults.get());
            results->storeRequestResult(request, tempResults->retrievePrice());
        } else {
            results->storeRequestResult(request, results->retrievePrice());
        }
    }
    Calculator() {};
};

/** Apply 'legal terms' */
bool LegalTerms::applyScenario(IObjectSP object){
    SensMgr sensMgr(object);
    sensMgr.shift(this);
    return sensMgr.getShiftStatus();
}

// Nothing to do before market data is retrieved.
bool LegalTerms::preapplyScenario(IObjectSP object){
    return false;
}

/** IPerturbation implementation - for backwards compatibility only.
    Equivalent here to applyScenario(objectToShift) */
bool LegalTerms::findAndShift(IObjectSP         objectToShift, 
                              OutputNameConstSP name){
    SensMgr sensMgr(objectToShift);
    sensMgr.shift(this);
    return sensMgr.getShiftStatus();
}

  
//// implementation of ITweakID
void LegalTerms::reset(){}

//// implementation of ITweakNameResolution - return null
ITweakNameResolver* LegalTerms::nameResolver(){
    return 0;
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP LegalTerms::shiftInterface() const{
    return Shift::TYPE;
}
 
bool LegalTerms::shift(IObjectSP obj) {
    // cast obj to LegalTerms::Shift and then invoke shift method
    Shift& LegalTermsObj = dynamic_cast<Shift&>(*obj);
    return LegalTermsObj.sensShift(this);
}
  
/** Creates an OutputRequestCalculator which can be used for a default
    implementation of the LegalTerms output request */
OutputRequestCalculator* LegalTerms::createCalculator(){
    return new Calculator();
}
  
class LegalTermsHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(LegalTerms, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IScenarioShift);
        IMPLEMENTS(IPerturbation); // for backwards compat
        FIELD(toTweak, "ignored - do not use");
        FIELD_MAKE_OPTIONAL(toTweak);
        EMPTY_SHELL_METHOD(defaultLegalTerms);
    }

    static IObject* defaultLegalTerms(){
        return new LegalTerms();
    }
};

CClassConstSP const LegalTerms::TYPE = CClass::registerClassLoadMethod(
    "LegalTerms", typeid(LegalTerms), LegalTermsHelper::load);

CClassConstSP const LegalTerms::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "LegalTerms::Shift", typeid(LegalTerms::Shift), 0);



DRLIB_END_NAMESPACE
