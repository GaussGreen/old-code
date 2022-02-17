#include "edginc/config.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/BoxedInt.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/fieldRiskProperty.hpp"
#include "edginc/TweakFunction.hpp"
#include "edginc/FieldPath.hpp"
#include "edginc/VolVSCurve.hpp"

DRLIB_BEGIN_NAMESPACE

struct VSCurveMRRTweak:
        public PerNameRiskPropertySensitivity<BoxedInt>,
        virtual public Additive {

    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;

    VSCurveMRRTweak(double shiftSize = DEFAULT_SHIFT):
        PerNameRiskPropertySensitivity<BoxedInt>(
            TYPE, shiftSize, NAME)
    {}

    Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     fieldRiskProperty::elementwise(
                         VolVSCurve::TYPE,
                         FieldPath::SP("meanReversRate"), FieldPathSP(),
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

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(VSCurveMRRTweak, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(PerNameRiskPropertySensitivity<BoxedInt>);
        EMPTY_SHELL_METHOD(&DefaultConstructor<VSCurveMRRTweak>::iObject);
        SensitivityFactory::addSens(NAME, 
                                    new GenericSensitivityFactory<VSCurveMRRTweak>(), 
                                    new VSCurveMRRTweak(DEFAULT_SHIFT),
                                    VolVSCurve::TYPE);
    }
};

const string VSCurveMRRTweak::NAME = "VSCURVE_MRR_TWEAK";
const double VSCurveMRRTweak::DEFAULT_SHIFT = 0.001;
const double VSCurveMRRTweak::SENSITIVITY_UNIT = 0.01;

CClassConstSP const VSCurveMRRTweak::TYPE = CClass::registerClassLoadMethod(
    "VSCurveMRRTweak", typeid(VSCurveMRRTweak), load);

bool VSCurveMRRTweakLinkIn() {
    return VSCurveMRRTweak::TYPE != NULL;
}

DRLIB_END_NAMESPACE
