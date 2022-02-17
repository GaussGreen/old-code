/**
 * @file CreditDeltaPointwiseWithShift.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/PropertyRiskAxis.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/RiskQuantity.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/ParSpreadPointwise.hpp"
#include "edginc/ParSpreadParallelRelative.hpp"


DRLIB_BEGIN_NAMESPACE

/**
 * A greek calculated by shifting all par spreads and then calculating a
 * pointwise tweak at the defined expiries.
 *
 * CreditDeltaPointwiseWithShift has been written to use the "declarative" 
 * sensitivities framework: see IRiskQuantityFactory for an overview.
 */

class CreditDeltaPointwiseWithShift:
        public PerNameRiskPropertySensitivity<ExpiryWindow>,
        public virtual Additive,
        public virtual CreditTweak {

public:
    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;

    CreditDeltaPointwiseWithShift(double shiftSize = DEFAULT_SHIFT,
                                  double parallelShiftSize = 0.10 /* 10% */):
        PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME),
        parallelShiftSize(parallelShiftSize)
    {}

    ~CreditDeltaPointwiseWithShift() {}

private:

    double parallelShiftSize;

    PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     RiskProperty<ParSpreadPointwise>::SP(),
                     IScalarDerivative::oneSided()->underScenario(
                         RiskProperty<ParSpreadParallelRelative>::SP(),
                         parallelShiftSize),
                     0.0001); // expressed in dPrice per BP
    }

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(CreditDeltaPointwiseWithShift, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
        EMPTY_SHELL_METHOD(&DefaultConstructor<CreditDeltaPointwiseWithShift>::iObject);
        FIELD(parallelShiftSize, "Shift size for parallel par spread scenario");
        FIELD_MAKE_OPTIONAL(parallelShiftSize);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<CreditDeltaPointwiseWithShift>(), 
                                    new CreditDeltaPointwiseWithShift(),
                                    ITweakableWithRespectTo<ParSpreadPointwise>::TYPE);
    }
};

const string CreditDeltaPointwiseWithShift::NAME = "CREDIT_DELTA_POINTWISE_WITH_SHIFT";
const double CreditDeltaPointwiseWithShift::DEFAULT_SHIFT = 0.0001;

CClassConstSP const CreditDeltaPointwiseWithShift::TYPE = 
    CClass::registerClassLoadMethod("CreditDeltaPointwiseWithShift", 
                                    typeid(CreditDeltaPointwiseWithShift), 
                                    load);

bool CreditDeltaPointwiseWithShiftLinkIn() {
    return CreditDeltaPointwiseWithShift::TYPE != NULL;
}

DRLIB_END_NAMESPACE
