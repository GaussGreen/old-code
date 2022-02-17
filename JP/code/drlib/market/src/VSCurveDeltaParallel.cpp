#include "edginc/config.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/VSCurveDeltaParallel.hpp"
#include "edginc/fieldRiskProperty.hpp"
#include "edginc/TweakFunction.hpp"
#include "edginc/FieldPath.hpp"
#include "edginc/VolVSCurve.hpp"


DRLIB_BEGIN_NAMESPACE

const string VSCurveDeltaParallel::NAME = "VSCURVE_DELTA_PARALLEL";
const double VSCurveDeltaParallel::DEFAULT_SHIFT = 0.0001;
const double VSCurveDeltaParallel::SENSITIVITY_UNIT = 0.01;

VSCurveDeltaParallel::VSCurveDeltaParallel(const string& name, double shiftSize):
    PerNameRiskPropertySensitivity<Void>(TYPE, shiftSize, name)
{}

VSCurveDeltaParallel::VSCurveDeltaParallel(double shiftSize):
    PerNameRiskPropertySensitivity<Void>(TYPE, shiftSize, NAME)
{}


PerNameRiskPropertySensitivity<Void>::Deriv VSCurveDeltaParallel::deriv() const {
    return Deriv(IResultsFunction::price(),
                 fieldRiskProperty::scalar(
                     VolVSCurve::TYPE,
                     IFieldTweak::bulk(
                         FieldPath::SP("valueCurve"), FieldPathSP(),
                         CDouble::SP(1.),
                         IFieldTweak::IOperator::numeric(
                             TweakFunction::additive(),
                             InfiniteRange::InfiniteRange(),
                             true, false)),
                     true),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}


VSCurveDeltaParallel::~VSCurveDeltaParallel() {}

void VSCurveDeltaParallel::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(VSCurveDeltaParallel, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<Void>);
    EMPTY_SHELL_METHOD(&DefaultConstructor<VSCurveDeltaParallel>::iObject);
    SensitivityFactory::addSens(VSCurveDeltaParallel::NAME, 
                                new GenericSensitivityFactory<VSCurveDeltaParallel>(), 
                                new VSCurveDeltaParallel(VSCurveDeltaParallel::DEFAULT_SHIFT),
                                VolVSCurve::TYPE);
}

CClassConstSP const VSCurveDeltaParallel::TYPE = CClass::registerClassLoadMethod(
    "VSCurveDeltaParallel", typeid(VSCurveDeltaParallel), load);

bool VSCurveDeltaParallelLinkIn() {
    return VSCurveDeltaParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
