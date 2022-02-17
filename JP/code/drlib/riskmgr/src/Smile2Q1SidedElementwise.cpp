/**
 * @file Smile2Q1SidedElementwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/BoxedInt.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/Smile2QElementwise.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/ParSpreadMaxTweakSize.hpp"
#include "edginc/IFieldRiskPropertyDefinition.hpp"
#include "edginc/FieldSensitivityDefinition.hpp"
#include "edginc/BoxedInt.hpp"

DRLIB_BEGIN_NAMESPACE

static const double ONE_BASIS_POINT = 0.0001;
static const double SENSITIVITY_UNIT = ONE_BASIS_POINT;

/**
 * A greek calculated by tweaking each market name's par spreads at all its
 * defined expiries.
 *
 * Smile2Q1SidedElementwise has been rewritten to use the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview.
 */

class Smile2Q1SidedElementwise: public PerNameRiskPropertySensitivity<BoxedInt>,
                             public virtual Additive {

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(Smile2Q1SidedElementwise, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(PerNameRiskPropertySensitivity<BoxedInt>);
        EMPTY_SHELL_METHOD(&DefaultConstructor<Smile2Q1SidedElementwise>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<Smile2Q1SidedElementwise>(), 
                                    new Smile2Q1SidedElementwise(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<Smile2QElementwise>::TYPE);
    }

    PerNameRiskPropertySensitivity<BoxedInt>::Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     RiskProperty<Smile2QElementwise>::SP(),
                     IScalarDerivative::oneSided(),
                     SENSITIVITY_UNIT);
    }

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const char* NAME;

    Smile2Q1SidedElementwise(double shiftSize = DEFAULT_SHIFT):
        PerNameRiskPropertySensitivity<BoxedInt>(TYPE, shiftSize, NAME)
    {}
};

const char* Smile2Q1SidedElementwise::NAME = "SMILE2Q_1SIDED_ELEMENTWISE";
const double Smile2Q1SidedElementwise::DEFAULT_SHIFT = 0.01;

PerNameRiskPropertySensitivity<BoxedInt>* newSmile2Q1SidedElementwise(double shiftSize) {
    return new Smile2Q1SidedElementwise(shiftSize);
}

CClassConstSP const Smile2Q1SidedElementwise::TYPE = CClass::registerClassLoadMethod(
    "Smile2Q1SidedElementwise", typeid(Smile2Q1SidedElementwise), load);

bool Smile2Q1SidedElementwiseLinkIn() {
    return Smile2Q1SidedElementwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
