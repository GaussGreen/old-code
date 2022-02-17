//////////////////////////////////////////////////////////////////////////

#include "edginc/config.hpp"
#include "edginc/SensControlAllNames.hpp"
#include "edginc/SensMgr.hpp"

DRLIB_BEGIN_NAMESPACE
SensControlAllNames::~SensControlAllNames(){}

/** Returns null (=> all names match)  */
ITweakNameResolver* SensControlAllNames::nameResolver(){
    return 0;
}


/** IPerturbation implementation - for backwards compatibility only.
    Equivalent here to applyScenario(objectToShift) */
bool SensControlAllNames::findAndShift(IObjectSP         objectToShift, 
                                       OutputNameConstSP name)
{
    return applyScenario(objectToShift);
}


/** apply this scenario shift to the supplied object */
bool SensControlAllNames::applyScenario(IObjectSP object) {
    SensMgr sensMgr(object);
    sensMgr.shift(this);
    return sensMgr.getShiftStatus();
}

// Nothing to do before market data is retrieved.
bool SensControlAllNames::preapplyScenario(IObjectSP object){
    return false;
}

/** Returns null as this tweak shifts all instances rather than
    specific named types */
OutputNameConstSP SensControlAllNames::getMarketDataName() const{
    return OutputNameConstSP();
}


/** Note SensControl is abstract. Create a sens control of type clazz
    and which uses outputName (eg VEGA_PARALLEL) to identify results */
SensControlAllNames::SensControlAllNames(const CClassConstSP& clazz,
                                         const string&        outputName):
    SensControl(clazz, outputName){}


void SensControlAllNames::load(CClassSP& clazz){
    REGISTER(SensControlAllNames, clazz);
    SUPERCLASS(SensControl);
    IMPLEMENTS(IScenarioShift);
    IMPLEMENTS(IPerturbation); // for backwards compatibility
}    

CClassConstSP const SensControlAllNames::TYPE = CClass::registerClassLoadMethod(
    "SensControlAllNames", typeid(SensControlAllNames), load);

DRLIB_END_NAMESPACE
