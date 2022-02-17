#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/RadarRepDealWrapper.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE


RadarRepDealWrapper::RadarRepDealWrapper(RadarRepDealSP newRadarRep) : CObject(TYPE), underlyingRadarRep(newRadarRep) {}

double RadarRepDealWrapper::getValue(DateTime date, DoubleArraySP fittingVals) {
    vector<double> fittingVec = vector<double>(fittingVals->begin(), fittingVals->end());
    return underlyingRadarRep->getRadarRepValue(date, fittingVec);
}

void RadarRepDealWrapper::validatePop2Object() {
}

void RadarRepDealWrapper::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("RadarRep object");
    REGISTER(RadarRepDealWrapper, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    //FIELD(radarMap, "Mapping between datetime and radar");

    Addin::registerConstructor(Addin::UTILITIES, RadarRepDealWrapper::TYPE);
}

CClassConstSP const RadarRepDealWrapper::TYPE = CClass::registerClassLoadMethod(
    "RadarRepDealWrapper", typeid(RadarRepDealWrapper), RadarRepDealWrapper::load);

/******************************/
// for type linking
bool RadarRepDealWrapperLoad(void){
    return (RadarRepDealWrapper::TYPE != 0);
}

IObject* RadarRepDealWrapper::clone() const {
    IObjectSP newRadarRep = RadarRepDealWrapperSP(new RadarRepDealWrapper(underlyingRadarRep));
    return newRadarRep.release();
}

DRLIB_END_NAMESPACE
