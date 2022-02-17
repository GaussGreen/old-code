/**
 * @file RiskPropertySensitivity.cpp
 */

#include "edginc/config.hpp"
#define QLIB_RISKPROPERTYSENSITIVITY_CPP
#include "edginc/Atomic.hpp"
#include "edginc/BoxedInt.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/ExpiryPair.hpp"
#include "edginc/ExpiryAndStrike.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/RiskPropertySensitivity.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/RiskMapping.hpp"

DRLIB_BEGIN_NAMESPACE

template <class QUALIFIER>
RiskPropertySensitivity<QUALIFIER>::RiskPropertySensitivity(
        CClassConstSP type,
        double shiftSize,
        const string& name,
        const string& name2):
    RiskQuantityFactorySensitivity(type, name, name2),
    shiftSize(shiftSize),
    hasPredefinedQualifiers(false)
{}

#if defined(_MSC_VER) && (_MSC_VER <1300)
// VC6 work around
#define IHATEVC6(QUALIFIER)                                                   \
    template <>                                                               \
    RiskPropertySensitivity<QUALIFIER>::Deriv::Deriv(                         \
            IResultsFunctionConstSP derivand,                                 \
            smartConstPtr<IRiskProperty<QUALIFIER> > property,                \
            IScalarDerivativeConstSP derivative,                              \
            double sensitivityUnit,                                           \
            IScalarDerivativeConstSP derivative2,                             \
            double sensitivityUnit2):                                         \
        derivand(derivand),                                                   \
        property(property),                                                   \
        derivatives(new IScalarDerivativeArray(                               \
                        1, IScalarDerivativeSP::constCast(derivative))),      \
        sensitivityUnits(new DoubleArray(1, sensitivityUnit))                 \
    {                                                                         \
        if (!!derivative2) {                                                  \
            /* good grief, const array nastiness */                           \
            IScalarDerivativeArraySP::constCast(derivatives)->push_back(      \
                IScalarDerivativeSP::constCast(derivative2));                 \
            DoubleArraySP::constCast(sensitivityUnits)->push_back(            \
                sensitivityUnit2);                                            \
        }                                                                     \
    }

IHATEVC6(Void);
IHATEVC6(ExpiryWindow);
IHATEVC6(ExpiryPair);
IHATEVC6(ExpiryAndStrike);
IHATEVC6(BoxedInt)
#endif


template <class QUALIFIER>
RiskPropertySensitivity<QUALIFIER>::Deriv::~Deriv(){}

// these two methods are explicit to avoid the need for the compiler to generate
// them everytime it sees this class
template <class QUALIFIER>
RiskPropertySensitivity<QUALIFIER>::Deriv::Deriv(const Deriv& rhs){
    *this = rhs;
}

template <class QUALIFIER>
typename RiskPropertySensitivity<QUALIFIER>::Deriv& RiskPropertySensitivity<QUALIFIER>::Deriv::operator=(const Deriv& rhs){
    this->derivand = rhs.derivand;
    this->property = rhs.property;
    this->derivatives = rhs.derivatives;
    this->sensitivityUnits = rhs.sensitivityUnits;
    return *this;
}


template <class QUALIFIER>
RiskPropertySensitivity<QUALIFIER>::Deriv::Deriv(
        IResultsFunctionConstSP derivand,
        smartConstPtr<IRiskProperty<QUALIFIER> > property,
        IScalarDerivativeConstSP derivative,
        double sensitivityUnit,
        IScalarDerivativeConstSP derivative2,
        double sensitivityUnit2):
    derivand(derivand),
    property(property),
    derivatives(new IScalarDerivativeArray(
                    1, IScalarDerivativeSP::constCast(derivative))),
    sensitivityUnits(new DoubleArray(1, sensitivityUnit))
{
    if (!!derivative2) {
        // good grief
        IScalarDerivativeArraySP::constCast(derivatives)->push_back(
            IScalarDerivativeSP::constCast(derivative2));
        DoubleArraySP::constCast(sensitivityUnits)->push_back(
            sensitivityUnit2);
    }
}

template <class QUALIFIER>
RiskPropertySensitivity<QUALIFIER>::~RiskPropertySensitivity() {}

template <class QUALIFIER>
void RiskPropertySensitivity<QUALIFIER>::ensureDeriv() const {
    if (!_derivand) {
        Deriv it(deriv());
        _derivand = it.derivand;
        _property = it.property;
        _derivatives = it.derivatives;
        _sensitivityUnits = it.sensitivityUnits;
    }
}

template <class QUALIFIER>
IResultsFunctionConstSP RiskPropertySensitivity<QUALIFIER>::derivand() const {
    ensureDeriv();
    return _derivand;
}

template <class QUALIFIER>
smartConstPtr<IRiskProperty<QUALIFIER> > RiskPropertySensitivity<QUALIFIER>::property() const {
    ensureDeriv();
    return _property;
}

template <class QUALIFIER>
IScalarDerivativeArrayConstSP RiskPropertySensitivity<QUALIFIER>::derivatives() const {
    ensureDeriv();
    return _derivatives;
}

template <class QUALIFIER>
DoubleArrayConstSP RiskPropertySensitivity<QUALIFIER>::sensitivityUnits() const {
    ensureDeriv();
    return _sensitivityUnits;
}

template <class QUALIFIER>
bool RiskPropertySensitivity<QUALIFIER>::discreteShift() const {
    return property()->discrete();
}

template <class QUALIFIER>
void RiskPropertySensitivity<QUALIFIER>::load(CClassSP& clazz) {
    REGISTER(RiskPropertySensitivity, clazz);
    SUPERCLASS(RiskQuantityFactorySensitivity);
    FIELD(_derivand, "derivand");
    FIELD_MAKE_TRANSIENT(_derivand);
    FIELD(_property, "property");
    FIELD_MAKE_TRANSIENT(_property);
    FIELD(_derivatives, "derivatives");
    FIELD_MAKE_TRANSIENT(_derivatives);
    FIELD(_sensitivityUnits, "sensitivityUnits");
    FIELD_MAKE_TRANSIENT(_sensitivityUnits);
    FIELD(shiftSize, "shiftSize");
    FIELD(predefinedQualifiers, "predefinedQualifiers");
    FIELD_MAKE_TRANSIENT(predefinedQualifiers);
    FIELD(hasPredefinedQualifiers, "hasPredefinedQualifiers");
    FIELD_MAKE_TRANSIENT(hasPredefinedQualifiers);
}

#if defined (_MSC_VER)

    // "no suitable definition provided" for availableRiskQuantities: yeah, it's pure virtual
#   pragma warning( disable : 4661 )

#endif

template <>
CClassConstSP const RiskPropertySensitivity<Void>::TYPE = CClass::registerClassLoadMethod(
    "RiskPropertySensitivity<Void>", typeid(RiskPropertySensitivity<Void>), load);

template <>
CClassConstSP const RiskPropertySensitivity<ExpiryWindow>::TYPE = CClass::registerClassLoadMethod(
    "RiskPropertySensitivity<ExpiryWindow>", typeid(RiskPropertySensitivity<ExpiryWindow>), load);

template <>
CClassConstSP const RiskPropertySensitivity<ExpiryPair>::TYPE = CClass::registerClassLoadMethod(
    "RiskPropertySensitivity<ExpiryPair>", typeid(RiskPropertySensitivity<ExpiryPair>), load);

template <>
CClassConstSP const RiskPropertySensitivity<ExpiryAndStrike>::TYPE = CClass::registerClassLoadMethod(
    "RiskPropertySensitivity<ExpiryAndStrike>", typeid(RiskPropertySensitivity<ExpiryAndStrike>), load);

template <>
CClassConstSP const RiskPropertySensitivity<BoxedInt>::TYPE = CClass::registerClassLoadMethod(
    "RiskPropertySensitivity<BoxedInt>", typeid(RiskPropertySensitivity<BoxedInt>), load);

/* The order/positioning of these lines is unbelievably important. In
   particular, a) the instantiation of the Deriv inner class needs to appear
   before the outer class is instantiated (at least I think) and b) before
   any code appears which would cause this template to be instantiated - in
   particular if the VC6 specific code was included on windows */

// Must use normal template instantiation since we want these instantiated
// for both debug and optimised
template class RISKMGR_DLL RiskPropertySensitivity<Void>::Deriv;
template class RISKMGR_DLL RiskPropertySensitivity<ExpiryWindow>::Deriv;
template class RISKMGR_DLL RiskPropertySensitivity<BoxedInt>::Deriv;
template class RISKMGR_DLL RiskPropertySensitivity<ExpiryPair>::Deriv;
template class RISKMGR_DLL RiskPropertySensitivity<ExpiryAndStrike>::Deriv;
template class RISKMGR_DLL RiskPropertySensitivity<Void>;
template class RISKMGR_DLL RiskPropertySensitivity<ExpiryWindow>;
template class RISKMGR_DLL RiskPropertySensitivity<BoxedInt>;
template class RISKMGR_DLL RiskPropertySensitivity<ExpiryPair>;
template class RISKMGR_DLL RiskPropertySensitivity<ExpiryAndStrike>;

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
