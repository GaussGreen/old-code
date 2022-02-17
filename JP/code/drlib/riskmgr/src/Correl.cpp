#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/BoxedEnum.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#define QLIB_CORREL_CPP
#include "edginc/Correl.hpp"

DRLIB_BEGIN_NAMESPACE

Correl::Correl(bool asset, bool fx, bool other, Operation op):
    CObject(TYPE),
    asset(asset),
    fx(fx),
    other(other),
    op(op)
{}

Correl::~Correl() {}

CorrelConstSP Correl::SP(bool asset, bool fx, bool other, Operation op) {
    return CorrelConstSP(new Correl(asset, fx, other, op));
}

static void Correl_load(CClassSP& clazz) {
    REGISTER(Correl, clazz);
    SUPERCLASS(CObject);
    FIELD(asset, "Whether to tweak asset-asset correlations");
    FIELD(fx, "Whether to tweak FX-FX and FX-other correlations");
    FIELD(other, "Whether to tweak other correlations");
    FIELD(op, "How to tweak");
    EMPTY_SHELL_METHOD(DefaultConstructor<Correl>::iObject);
}

CClassConstSP const Correl::TYPE = CClass::registerClassLoadMethod("Correl", typeid(Correl), Correl_load);

RiskProperty_TYPES(Correl)

INSTANTIATE_TEMPLATE(class RISKMGR_DLL ITweakableWithRespectTo<Correl>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL IRestorableWithRespectTo<Correl>);

START_PUBLIC_ENUM_DEFINITION(Correl::Operation, "How correlations are to be tweaked");
ENUM_VALUE_AND_NAME(Correl::ABSOLUTE_INWARD, "Absolute (inward)",
                    "shiftSize is amount by which correl is moved towards zero");
ENUM_VALUE_AND_NAME(Correl::ABSOLUTE, "Absolute",
                    "shiftSize is amount correl is moved");
ENUM_VALUE_AND_NAME(Correl::SQUEEZE, "Squeeze",
                    "shiftSize is the proportion of the way correl is moved towards one");
END_ENUM_DEFINITION(Correl::Operation);

DRLIB_END_NAMESPACE
