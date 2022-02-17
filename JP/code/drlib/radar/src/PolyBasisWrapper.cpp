#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/PolyBasisWrapper.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

IFunctionBasisSP PolyBasisWrapper::getBasis(void) const {
    return IFunctionBasisSP(new PolynomialBasis(numVars, maxDeg));
}

PolyBasisWrapper::PolyBasisWrapper(int numVars, int maxDeg) :
    CObject(TYPE), numVars(numVars), maxDeg(maxDeg)
{    
}

void PolyBasisWrapper::validatePop2Object() {
}

void PolyBasisWrapper::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("PolyBasis object");
    REGISTER(PolyBasisWrapper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IFuncBasisWrapper);

    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(numVars, "The no. of variables");
    FIELD(maxDeg, "Max degree of the polynomials in the basis");

    Addin::registerConstructor(Addin::UTILITIES, PolyBasisWrapper::TYPE);
}

CClassConstSP const PolyBasisWrapper::TYPE = CClass::registerClassLoadMethod(
    "PolynomialBasis", typeid(PolyBasisWrapper), PolyBasisWrapper::load);

/******************************/
// for type linking
bool PolyBasisWrapperLoad(void){
    return (PolyBasisWrapper::TYPE != 0);
}

DRLIB_END_NAMESPACE
