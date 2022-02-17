#include "edginc/config.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/VSCurveDeltaPointwise.hpp"
#include "edginc/fieldRiskProperty.hpp"
#include "edginc/TweakFunction.hpp"
#include "edginc/FieldPath.hpp"
#include "edginc/VolVSCurve.hpp"


DRLIB_BEGIN_NAMESPACE

const string VSCurveDeltaPointwise::NAME = "VSCURVE_DELTA_POINTWISE";
const double VSCurveDeltaPointwise::DEFAULT_SHIFT = 0.0001;
const double VSCurveDeltaPointwise::SENSITIVITY_UNIT = 0.01;

VSCurveDeltaPointwise::VSCurveDeltaPointwise(const string& name, double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, name)
{}

VSCurveDeltaPointwise::VSCurveDeltaPointwise(double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
{}



PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv VSCurveDeltaPointwise::deriv() const {
    return Deriv(IResultsFunction::price(),
                 fieldRiskProperty::pointwise(
                     VolVSCurve::TYPE,
                     FieldPath::SP("valueCurve"), FieldPathSP(),
                     FieldPath::SP("tenorCurve"), FieldPathSP(),
                     0,
                     IFieldTweak::IOperator::numeric(
                         TweakFunction::additive(),
                         InfiniteRange::InfiniteRange(),
                         true, false),
                     CDouble::SP(1.),
                     true),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}





VSCurveDeltaPointwise::~VSCurveDeltaPointwise() {}

void VSCurveDeltaPointwise::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(VSCurveDeltaPointwise, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
    EMPTY_SHELL_METHOD(&DefaultConstructor<VSCurveDeltaPointwise>::iObject);
    SensitivityFactory::addSens(VSCurveDeltaPointwise::NAME, 
                                new GenericSensitivityFactory<VSCurveDeltaPointwise>(), 
                                new VSCurveDeltaPointwise(VSCurveDeltaPointwise::DEFAULT_SHIFT),
                                VolVSCurve::TYPE);
}

CClassConstSP const VSCurveDeltaPointwise::TYPE = CClass::registerClassLoadMethod(
    "VSCurveDeltaPointwise", typeid(VSCurveDeltaPointwise), load);

bool VSCurveDeltaPointwiseLinkIn() {
    return VSCurveDeltaPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE
