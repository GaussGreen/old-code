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

struct VSCurveVVega:
        public PerNameRiskPropertySensitivity<BoxedInt>,
        virtual public Additive {

    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;

    VSCurveVVega(double shiftSize = DEFAULT_SHIFT):
        PerNameRiskPropertySensitivity<BoxedInt>(
            TYPE, shiftSize, NAME)
    {}

    Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     fieldRiskProperty::elementwise(
                         VolVSCurve::TYPE,
                         FieldPath::SP("volVol"), FieldPathSP(),
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
        REGISTER(VSCurveVVega, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(PerNameRiskPropertySensitivity<BoxedInt>);
        EMPTY_SHELL_METHOD(&DefaultConstructor<VSCurveVVega>::iObject);
        SensitivityFactory::addSens(NAME, 
                                    new GenericSensitivityFactory<VSCurveVVega>(), 
                                    new VSCurveVVega(DEFAULT_SHIFT),
                                    VolVSCurve::TYPE);
    }
};

const string VSCurveVVega::NAME = "VSCURVE_VVEGA";
const double VSCurveVVega::DEFAULT_SHIFT = 0.001;
const double VSCurveVVega::SENSITIVITY_UNIT = 0.01;

CClassConstSP const VSCurveVVega::TYPE = CClass::registerClassLoadMethod(
    "VSCurveVVega", typeid(VSCurveVVega), load);

bool VSCurveVVegaLinkIn() {
    return VSCurveVVega::TYPE != NULL;
}

DRLIB_END_NAMESPACE
