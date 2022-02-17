//////////////////////////////////////////////////////////////////////////

#include "edginc/config.hpp"
#include "edginc/SensControlPerName.hpp"
#include "edginc/Results.hpp"
#include "edginc/NotApplicable.hpp"
#include "edginc/SensMgr.hpp"

DRLIB_BEGIN_NAMESPACE

SensControlPerName::~SensControlPerName(){}

/** sets the name identifying the market data to be shifted */
void SensControlPerName::setMarketDataName(OutputNameConstSP name){
    marketDataName = name;
}

/** does the given name match the name identifying the market data
    to be shifted */
bool SensControlPerName::marketDataNameMatches(OutputNameConstSP name) const{
    return name->equals(marketDataName.get());
}

/** Returns array of output names which need to be tweaked for this
    sensitivity. In particular, if hasOverrideNames is true then returns
    overrideNames() else generates a list based on object supplied (which
    is typically, but needn't be, a TweakGroup */
OutputNameArrayConstSP SensControlPerName::names(const IObject* tweakGroup) const{
    return OutputName::trim(hasOverrideNames()? 
                            overrideNames():
                            SensMgrConst(tweakGroup).allNames(
                                const_cast<SensControlPerName*>(this)));
}

OutputNameArrayConstSP SensControlPerName::allNames(const IObject* tweakGroup) const{
    SensMgrConst sensMgr(tweakGroup);
    SensControlPerName* mthis = const_cast<SensControlPerName*>(this);
    return sensMgr.allNames(mthis);
}

/** The same as the above method except that for the case where
    this has its own names then for each of those names which is
    not present NotApplicable is stored in the Results. */
OutputNameArrayConstSP SensControlPerName::names(const IObject* tweakGroup,
                                                 Results*       results){
    if (!hasOverrideNames()) {
        return names(tweakGroup);
    }
    OutputNameArrayConstSP explicitNames(OutputName::trim(overrideNames()));
    // get all names actually present
    OutputNameArrayConstSP names(OutputName::trim(
                                     SensMgrConst(tweakGroup).allNames(this)));
    // now identify any extra names
    OutputNameArraySP extraNames(OutputName::difference(explicitNames,
                                                        names));
    for (int i = 0; i < extraNames->size(); i++){
        results->storeGreek(IObjectSP(new NotApplicable()),
                            getPacketName(),
                            (*extraNames)[i]);
    }
    return explicitNames;
}

/** apply this Perturbation to the object with the specified name contained
    within the supplied object */
bool SensControlPerName::findAndShift(IObjectSP         objectToShift, 
                                      OutputNameConstSP name){
    OutputNameConstSP tmpName = marketDataName;
    try{
        setMarketDataName(name);
        SensMgr sensMgr(objectToShift);
        sensMgr.shift(this);
        setMarketDataName(tmpName);
        return sensMgr.getShiftStatus();
    } catch (exception&){
        setMarketDataName(tmpName);
        throw;
    } 
}

/** Returns this => tweaks defined by classes derived from this type are
    associated with individual names */
ITweakNameResolver* SensControlPerName::nameResolver(){
    return this;
}

/** returns the name identifying the market data to be shifted. Returns
    null if not set */
OutputNameConstSP SensControlPerName::getMarketDataName() const{
    return marketDataName;
}

/** Note SensControl is abstract. Create a sens control of type clazz
    and which uses outputName (eg VEGA_PARALLEL) to identify results */
SensControlPerName::SensControlPerName(const CClassConstSP& clazz,
                                       const string&        outputName):
    SensControl(clazz, outputName){}

void SensControlPerName::load(CClassSP& clazz){
    REGISTER(SensControlPerName, clazz);
    SUPERCLASS(SensControl);
    IMPLEMENTS(IPerturbation);
    IMPLEMENTS(IPerNameSensitivity);
    FIELD(marketDataName, "Market data name");
    FIELD_MAKE_TRANSIENT(marketDataName);
}

CClassConstSP const SensControlPerName::TYPE = CClass::registerClassLoadMethod(
    "SensControlPerName", typeid(SensControlPerName), load);


DRLIB_END_NAMESPACE
