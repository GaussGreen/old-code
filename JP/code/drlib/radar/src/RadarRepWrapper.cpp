#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/RadarRepWrapper.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

/*RadarRepWrapper::RadarRepWrapper(DoubleArraySP coef, IFuncBasisWrapperSP funcBasis, IFittingVarTransWrapperSP fittingTransform) :
    CObject(TYPE), coefArr(coef), funcBasis(funcBasis), fittingTransform(fittingTransform) {
}*/

RadarRepWrapper::RadarRepWrapper(RadarRepSP radar) : CObject(TYPE), underlyingRadarRep(radar) {
    
}

double RadarRepWrapper::getValue(DoubleArraySP fittingVals) {
    /*vector<double> coefVec = vector<double>(coefArr->begin(), coefArr->end());
    underlyingRadarRep = RadarRepSP(new RadarRep(coefVec, funcBasis->getBasis(), fittingTransform->getTransform()));*/
    FittingArray fArr = FittingArray(fittingVals->begin(), fittingVals->end());
    return underlyingRadarRep->getValue(fArr);
}

void RadarRepWrapper::validatePop2Object() {
}

void RadarRepWrapper::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("RadarRep object");
    REGISTER(RadarRepWrapper, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);

    /*FIELD(coefArr, "Fitting variable product vector");
    FIELD(funcBasis, "What basis functions are being used");
    FIELD_MAKE_OPTIONAL(funcBasis);
    FIELD(fittingTransform, "What fitting variable transformation is being used");*/

    Addin::registerConstructor(Addin::UTILITIES, RadarRepWrapper::TYPE);
}

CClassConstSP const RadarRepWrapper::TYPE = CClass::registerClassLoadMethod(
    "RadarRepWrapper", typeid(RadarRepWrapper), RadarRepWrapper::load);

IObject* RadarRepWrapper::clone() const {
    IObjectSP newRadarRep = RadarRepWrapperSP(new RadarRepWrapper(underlyingRadarRep));
    return newRadarRep.release();
}

/******************************/
// for type linking
bool RadarRepWrapperLoad(void){
    return (RadarRepWrapper::TYPE != 0);
}

DRLIB_END_NAMESPACE
