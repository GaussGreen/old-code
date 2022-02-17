#include "edginc/config.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/ICrossDerivative.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/CrossRiskPropertySensitivity.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/VSCurveDeltaPointwise.hpp"
#include "edginc/fieldRiskProperty.hpp"
#include "edginc/TweakFunction.hpp"
#include "edginc/FieldPath.hpp"
#include "edginc/VolVSCurve.hpp"

DRLIB_BEGIN_NAMESPACE

struct VSCurveCrossGamma:
        public CrossRiskPropertySensitivity<ExpiryWindow, ExpiryWindow>,
        virtual public Additive {

    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;

    VSCurveCrossGamma(double shiftSize = DEFAULT_SHIFT):
        CrossRiskPropertySensitivity<ExpiryWindow, ExpiryWindow>(
            TYPE, shiftSize, shiftSize, NAME)
    {}

    Deriv deriv() const {
        IExpiryRiskPropertyConstSP prop = fieldRiskProperty::pointwise(
            VolVSCurve::TYPE,
            FieldPath::SP("valueCurve"), FieldPathSP(),
            FieldPath::SP("tenorCurve"), FieldPathSP(),
            0,
            IFieldTweak::IOperator::numeric(
                TweakFunction::additive(),
                InfiniteRange::InfiniteRange(),
                true, false),
            IObjectSP(CDouble::create(1.)),
            true);

        return Deriv(IResultsFunction::price(), prop, prop,
                     ICrossDerivative::cross(), SENSITIVITY_UNIT);
    }

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(VSCurveCrossGamma, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(CrossRiskPropertySensitivity<ExpiryWindow _COMMA_ ExpiryWindow>);
        EMPTY_SHELL_METHOD(&DefaultConstructor<VSCurveCrossGamma>::iObject);
        SensitivityFactory::addSens(NAME, 
                                    new GenericSensitivityFactory<VSCurveCrossGamma>(), 
                                    new VSCurveCrossGamma(DEFAULT_SHIFT),
                                    VolVSCurve::TYPE);
    }
};

const string VSCurveCrossGamma::NAME = "VSCURVE_CROSS_GAMMA";
const double VSCurveCrossGamma::DEFAULT_SHIFT = VSCurveDeltaPointwise::DEFAULT_SHIFT;
const double VSCurveCrossGamma::SENSITIVITY_UNIT =
    VSCurveDeltaPointwise::SENSITIVITY_UNIT * VSCurveDeltaPointwise::SENSITIVITY_UNIT;

CClassConstSP const VSCurveCrossGamma::TYPE = CClass::registerClassLoadMethod(
    "VSCurveCrossGamma", typeid(VSCurveCrossGamma), load);

bool VSCurveCrossGammaLinkIn() {
    return VSCurveCrossGamma::TYPE != NULL;
}

DRLIB_END_NAMESPACE
