#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/Addin.hpp"
#include "edginc/RadarRepAddin.hpp"

DRLIB_BEGIN_NAMESPACE

RadarRepAddin::RadarRepAddin(RadarRepDealWrapperSP p) :
CObject(TYPE), radarDealW(p) 
{
}

//double RadarRepAddin::getValue(DoubleArraySP dArr) const {
double RadarRepAddin::getValue() 
{
    return radarDealW->getValue(radarDate,fittingVals);
}

void RadarRepAddin::validatePop2Object() {
}

void RadarRepAddin::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("RadarRep object");
    REGISTER(RadarRepAddin, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(radarDealW, "RadarRep wrapper field");
    FIELD(fittingVals, "Values of the fitting variables");
    FIELD(radarDate, "Date on which radar values will be calculated");

    Addin::registerConstructor(Addin::UTILITIES, RadarRepAddin::TYPE);

    Addin::registerDoubleMethod("GET_RADAR_REP_VALUE",
        Addin::RISK,
        "Returns a specific radar value",
        &RadarRepAddin::getValue);
}

CClassConstSP const RadarRepAddin::TYPE = CClass::registerClassLoadMethod(
    "RadarRepAddin", typeid(RadarRepAddin), RadarRepAddin::load);

/******************************/
// for type linking
bool RadarRepAddinLoad(void){
    return (RadarRepAddin::TYPE != 0);
}

DRLIB_END_NAMESPACE
