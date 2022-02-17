/**
 * @file ParSpreadRhoPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/ParSpreadPointwise.hpp"
#include "edginc/ParSpreadRhoPointwise.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/ParSpreadMaxTweakSize.hpp"
#include "edginc/IFieldRiskPropertyDefinition.hpp"
#include "edginc/FieldSensitivityDefinition.hpp"

DRLIB_BEGIN_NAMESPACE

static const double ONE_BASIS_POINT = 0.0001;
static const double SENSITIVITY_UNIT = ONE_BASIS_POINT;

/**
 * A greek calculated by tweaking each market name's par spreads at all its
 * defined expiries.
 *
 * ParSpreadRhoPointwise has been rewritten to use the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview.
 */

class ParSpreadRhoPointwise: public PerNameRiskPropertySensitivity<ExpiryWindow>,
                             public virtual Additive,
                             public CreditTweak {

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(ParSpreadRhoPointwise, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
        EMPTY_SHELL_METHOD(&DefaultConstructor<ParSpreadRhoPointwise>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<ParSpreadRhoPointwise>(), 
                                    new ParSpreadRhoPointwise(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<ParSpreadPointwise>::TYPE);
    }

    PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     ParSpreadMaxTweakSize::adapted(
                         RiskProperty<ParSpreadPointwise>::SP()),
                     IScalarDerivative::oneSided(),
                     SENSITIVITY_UNIT);
    }

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const char* NAME;

    ParSpreadRhoPointwise(double shiftSize = DEFAULT_SHIFT):
        PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
    {}
};

const char* ParSpreadRhoPointwise::NAME = "PAR_SPREAD_RHO_POINTWISE";

PerNameRiskPropertySensitivity<ExpiryWindow>* newParSpreadRhoPointwise(double shiftSize) {
    return new ParSpreadRhoPointwise(shiftSize);
}

CClassConstSP const ParSpreadRhoPointwise::TYPE = CClass::registerClassLoadMethod(
    "ParSpreadRhoPointwise", typeid(ParSpreadRhoPointwise), load);

FORWARD_DECLARE(ParSpreadRhoPointwiseTwoSided)

class ParSpreadRhoPointwiseTwoSided: public PerNameRiskPropertySensitivity<ExpiryWindow>,
                                     public virtual Additive,
                                     public CreditTweak {

    static IExpiryRiskPropertyConstSP builtinProp(IObjectConstSP arg) {
        return RiskProperty<ParSpreadPointwise>::SP();
    }

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(ParSpreadRhoPointwiseTwoSided, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
        EMPTY_SHELL_METHOD(&DefaultConstructor<ParSpreadRhoPointwiseTwoSided>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<ParSpreadRhoPointwiseTwoSided>(), 
                                    new ParSpreadRhoPointwiseTwoSided(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<ParSpreadPointwise>::TYPE);

        // We use this as a test case for the "generic greeks" internal override mechanism

        FieldSensitivityDefinition::registerBuiltin(
            "2-sided ParSpreadRhoPointwise",
            RQFFactory_2packet<ParSpreadRhoPointwiseTwoSided>::newOne);

        IFieldRiskPropertyDefinition::registerPointwiseBuiltin(
            "par spread pointwise", builtinProp);
    }

    PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     ParSpreadMaxTweakSize::adapted(
                         RiskProperty<ParSpreadPointwise>::SP()),
                     IScalarDerivative::twoSided(),
                     SENSITIVITY_UNIT,
                     IScalarDerivative::second(),
                     SENSITIVITY_UNIT);
    }

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;

    static const char* NAME;
    static const char* NAME2;

    ParSpreadRhoPointwiseTwoSided(
            double shiftSize = DEFAULT_SHIFT,
            const string& packetName1 = NAME,
            const string& packetName2 = NAME2):
        PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, packetName1, packetName2)
    {}
};

const char* ParSpreadRhoPointwiseTwoSided::NAME = "PAR_SPREAD_RHO_POINTWISE_2SIDED";
const char* ParSpreadRhoPointwiseTwoSided::NAME2 = "PAR_SPREAD_RHO_POINTWISE_GAMMA";

CClassConstSP const ParSpreadRhoPointwiseTwoSided::TYPE = CClass::registerClassLoadMethod(
    "ParSpreadRhoPointwiseTwoSided", typeid(ParSpreadRhoPointwiseTwoSided), load);

const double ParSpreadRhoPointwise::DEFAULT_SHIFT = ONE_BASIS_POINT;
const double ParSpreadRhoPointwiseTwoSided::DEFAULT_SHIFT = ONE_BASIS_POINT;

bool ParSpreadRhoPointwiseLinkIn() {
    return ParSpreadRhoPointwise::TYPE != NULL &&
           ParSpreadRhoPointwiseTwoSided::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
