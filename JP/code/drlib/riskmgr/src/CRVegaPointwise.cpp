/**
 * @file CRVegaPointwise.cpp
 * 
 * Group : Credit QR NY
 * 
 * File  : CRVegaPointwise.cpp
 * 
 * Description : the tweak tag that describes a pointwise tweak for spread volatility cubes. The greek is calculated by tweaking each market volatility smile at all its
 * defined expiries.
 * 
 */

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/PropertyRiskAxis.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/RiskQuantity.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/CRVolPointwise.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/CRVegaPointwise.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

const string CRVegaPointwise::NAME = "CR_VEGA_POINTWISE";
const double CRVegaPointwise::DEFAULT_SHIFT = 0.001;
const double CRVegaPointwise::SENSITIVITY_UNIT = 0.01;

CRVegaPointwise::CRVegaPointwise(double shiftSize):
    PerNameRiskPropertySensitivity<ExpiryPair>(TYPE, shiftSize, NAME)
{}

PerNameRiskPropertySensitivity<ExpiryPair>::Deriv CRVegaPointwise::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<CRVolPointwise>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

CRVegaPointwise::~CRVegaPointwise() {}

void CRVegaPointwise::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CRVegaPointwise, clazz);
    IMPLEMENTS(Additive);
    SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryPair>);
    EMPTY_SHELL_METHOD(&DefaultConstructor<CRVegaPointwise>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<CRVegaPointwise>(), 
                                new CRVegaPointwise(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<CRVolPointwise>::TYPE);
}

CClassConstSP const CRVegaPointwise::TYPE = CClass::registerClassLoadMethod(
    "CRVegaPointwise", typeid(CRVegaPointwise), load);

bool CRVegaPointwiseLinkIn() {
    return CRVegaPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

