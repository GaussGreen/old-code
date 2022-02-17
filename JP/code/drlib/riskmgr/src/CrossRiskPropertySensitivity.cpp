/**
 * @file CrossRiskPropertySensitivity.cpp
 */

#include "edginc/config.hpp"
#define QLIB_CROSSRISKPROPERTYSENSITIVITY_CPP
#include "edginc/TRACE.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/ExpiryPair.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/CrossRiskPropertySensitivity.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/RiskMapping.hpp"

DRLIB_BEGIN_NAMESPACE

template <class Q0, class Q1>
CrossRiskPropertySensitivity<Q0, Q1>::CrossRiskPropertySensitivity(
        CClassConstSP type,
        double shiftSize0,
        double shiftSize1,
        const string& name):
    RiskQuantityFactorySensitivity(type, name),
    shiftSize(shiftSize0),
    shiftSize1(shiftSize1),
    _sensitivityUnit(0) // to avoid UMRs
{}

template <class Q0, class Q1>
CrossRiskPropertySensitivity<Q0, Q1>::~CrossRiskPropertySensitivity() {}

template <class Q0, class Q1>
void CrossRiskPropertySensitivity<Q0, Q1>::ensureDeriv() const {
    if (!_derivand) {
        Deriv it = deriv();
        _derivand = it.derivand;
        _property0 = it.property0;
        _property1 = it.property1;
        _derivative = it.derivative;
        _sensitivityUnit = it.sensitivityUnit;
    }
}

template <class Q0, class Q1>
IResultsFunctionConstSP CrossRiskPropertySensitivity<Q0, Q1>::derivand() const {
    ensureDeriv();
    return _derivand;
}

template <class Q0, class Q1>
smartConstPtr<IRiskProperty<Q0> > CrossRiskPropertySensitivity<Q0, Q1>::property0() const {
    ensureDeriv();
    return _property0;
}

template <class Q0, class Q1>
smartConstPtr<IRiskProperty<Q1> > CrossRiskPropertySensitivity<Q0, Q1>::property1() const {
    ensureDeriv();
    return _property1;
}

template <class Q0, class Q1>
ICrossDerivativeConstSP CrossRiskPropertySensitivity<Q0, Q1>::derivative() const {
    ensureDeriv();
    return _derivative;
}

template <class Q0, class Q1>
double CrossRiskPropertySensitivity<Q0, Q1>::sensitivityUnit() const {
    ensureDeriv();
    return _sensitivityUnit;
}

template <class Q0, class Q1>
bool CrossRiskPropertySensitivity<Q0, Q1>::discreteShift() const {
    return property0()->discrete() && property1()->discrete();
}

template <class Q0, class Q1>
CrossRiskPropertySensitivity<Q0, Q1>::Deriv::Deriv(
        IResultsFunctionConstSP derivand,
        smartConstPtr<IRiskProperty<Q0> > property0,
        smartConstPtr<IRiskProperty<Q1> > property1,
        ICrossDerivativeConstSP derivative,
        double sensitivityUnit):
    derivand(derivand),
    property0(property0),
    property1(property1),
    derivative(derivative),
    sensitivityUnit(sensitivityUnit)
{}

template <class Q0, class Q1>
NamedRiskQuantityArraySP CrossRiskPropertySensitivity<Q0, Q1>::nameRiskQuantities(
        MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping) const {

    TRACE_METHOD;

    OutputNameArrayConstSP inWorld0 =
        OutputName::trim(riskMapping->subjectNames(property0(), world));

    OutputNameArrayConstSP inWorld1 =
        OutputName::trim(riskMapping->subjectNames(property1(), world));

    OutputNameArrayConstSP tweakables0 = hasOverrideNames() ?
        OutputName::intersection(inWorld0, overrideNames()) : inWorld0;

    OutputNameArrayConstSP tweakables1 = hasOverrideNames() ?
        OutputName::intersection(inWorld1, overrideNames()) : inWorld1;

    OutputNameArrayConstSP tweakables =
        OutputName::intersection(tweakables0, tweakables1);

    NamedRiskQuantityArraySP rqs(new NamedRiskQuantityArray());

    for (int t = 0; t < tweakables->size(); ++t) {
        TRACE_BLOCK("Market name " << *(*tweakables)[t]);
        try {
            Qualifier0ArrayConstSP quals0 =
                riskMapping->subjectQualifiers(property0(), (*tweakables)[t],
                                               world);

            Qualifier1ArrayConstSP quals1 =
                riskMapping->subjectQualifiers(property1(), (*tweakables)[t],
                                               world);

/**/        this->marketDataName = (*tweakables)[t]; // so that CDSPS::endDate() works

            BoolArrayConstSP isSens0 = property0()->mayHaveEffect(
                world, (*tweakables)[t], quals0, this);
            ASSERT(!!isSens0 && isSens0->size() == quals0->size());

            DoubleArrayConstSP coeffs0 = property0()->adaptiveCoefficients(
                world, (*tweakables)[t], quals0, shiftSize);
            ASSERT(!!coeffs0 && coeffs0->size() == quals0->size());

            BoolArrayConstSP isSens1 = property1()->mayHaveEffect(
                world, (*tweakables)[t], quals1, this);
            ASSERT(!!isSens1 && isSens1->size() == quals1->size());

            DoubleArrayConstSP coeffs1 = property1()->adaptiveCoefficients(
                world, (*tweakables)[t], quals1, shiftSize1);
            ASSERT(!!coeffs1 && coeffs1->size() == quals1->size());

            for (int q0 = 0; q0 < quals0->size(); ++q0) {
                IRiskAxisConstSP axis0 = property0()->axisFor((*tweakables)[t],
                                                              (*quals0)[q0]);
                for (int q1 = 0; q1 <= q0; ++q1) {
                    IRiskAxisConstSP axis1 = property1()->axisFor((*tweakables)[t],
                                                                  (*quals1)[q1]);

                    rqs->push_back(NamedRiskQuantity::SP(
                        (*isSens0)[q0] && (*isSens1)[q1] ?
                            derivative()->riskQuantity(derivand(),
                                                       axis0, (*coeffs0)[q0],
                                                       axis1, (*coeffs1)[q1]) :
                            RiskQuantity::constant(0),
                        IResultsIdentifier::SP(
                            getSensOutputName(),
                            (*tweakables)[t], (*quals0)[q0], (*quals1)[q1]),
                        sensitivityUnit()));

                    if ((*isSens0)[q0] && (*isSens1)[q1]) {
                        TRACE("Including request to compute " <<
                              *rqs->back()->resultsName);
                    }
                    else {
                        TRACE("Including request to just set " <<
                              *rqs->back()->resultsName <<
                              " to zero (known to be insensitive)");
                    }
                }
            }
        }
        catch (exception& e) {
            rqs->push_back(NamedRiskQuantity::SP(
                               RiskQuantity::untweakable(e),
                               IResultsIdentifier::SP(getSensOutputName(),
                                                      (*tweakables)[t])));

            TRACE("Including request to mark " <<
                  *rqs->back()->resultsName << " as Untweakable");
        }
    }

/**/this->marketDataName.reset();

    if (hasOverrideNames()) {
        OutputNameArrayConstSP na = OutputName::trim(OutputName::difference(
            overrideNames(), tweakables));
        for (int n = 0; n < na->size(); ++n) {
            rqs->push_back(NamedRiskQuantity::SP(
                RiskQuantity::notApplicable(),
                // FIXME probably wrong in some sense
                IResultsIdentifier::SP(getPacketName(),
                                       (*na)[n])));
        }
    }

    return rqs;
}

template <class Q0, class Q1>
OutputNameConstSP CrossRiskPropertySensitivity<Q0, Q1>::getMarketDataName() const {
    return marketDataName;
}

// A bit impressionistic

template <class Q0, class Q1>
OutputNameArrayConstSP CrossRiskPropertySensitivity<Q0, Q1>::allNames(
        const IObject* object) const {
    IObjectConstSP o = IObjectConstSP::attachToRef(object);
    return OutputName::intersection(property0()->subjectNames(o),
                                    property1()->subjectNames(o));
}

// A bit impressionistic

template <class Q0, class Q1>
CClassConstSP CrossRiskPropertySensitivity<Q0, Q1>::shiftInterface() const {
    return property0()->subjectInterface();
}

template <class Q0, class Q1>
void CrossRiskPropertySensitivity<Q0, Q1>::load(CClassSP& clazz) {
    REGISTER(CrossRiskPropertySensitivity, clazz);
    SUPERCLASS(RiskQuantityFactorySensitivity);
    IMPLEMENTS(IPerNameSensitivity);
    FIELD(_derivand, "derivand");
    FIELD_MAKE_TRANSIENT(_derivand);
    FIELD(_property0, "property0");
    FIELD_MAKE_TRANSIENT(_property0);
    FIELD(_property1, "property1");
    FIELD_MAKE_TRANSIENT(_property1);
    FIELD(_derivative, "derivative");
    FIELD_MAKE_TRANSIENT(_derivative);
    FIELD(_sensitivityUnit, "sensitivityUnit");
    FIELD_MAKE_TRANSIENT(_sensitivityUnit);
    FIELD(shiftSize, "shiftSize");
    FIELD_MAKE_OPTIONAL(shiftSize);
    FIELD(shiftSize1, "shiftSize1");
    FIELD_MAKE_OPTIONAL(shiftSize1);
}

template <>
CClassConstSP const CrossRiskPropertySensitivity<Void, Void>::TYPE = CClass::registerClassLoadMethod(
    "CrossRiskPropertySensitivity<Void, Void>", typeid(CrossRiskPropertySensitivity<Void, Void>), load);

template <>
CClassConstSP const CrossRiskPropertySensitivity<ExpiryWindow, ExpiryWindow>::TYPE = CClass::registerClassLoadMethod(
    "CrossRiskPropertySensitivity<ExpiryWindow, ExpiryWindow>", typeid(CrossRiskPropertySensitivity<ExpiryWindow, ExpiryWindow>), load);
/* The order/positioning of these lines is unbelievably important. In
   particular, a) the instantiation of the Deriv inner class needs to appear
   before the outer class is instantiated (at least I think) and b) before
   any code appears which would cause this template to be instantiated - in
   particular the instantiation of the GeneralCrossRiskPropertySensitivity 
   template below. */
// Must use normal template instantiation since we want these instantiated
// for both debug and optimised
template struct RISKMGR_DLL CrossRiskPropertySensitivity<Void, Void>::Deriv;
template struct RISKMGR_DLL 
CrossRiskPropertySensitivity<ExpiryWindow, ExpiryWindow>::Deriv;

template <class Q0, class Q1>
class GeneralCrossRiskPropertySensitivity:
        public CrossRiskPropertySensitivity<Q0, Q1>,
        public virtual Additive {

    static CClassConstSP const TYPE;

private:

    GeneralCrossRiskPropertySensitivity():
        CrossRiskPropertySensitivity<Q0, Q1>(TYPE, 0., 0., "") {}

    static IObject* emptyShell() {
        return new GeneralCrossRiskPropertySensitivity();
    }        

    static void load(CClassSP& clazz) {
        REGISTER(GeneralCrossRiskPropertySensitivity, clazz);
        SUPERCLASS(CrossRiskPropertySensitivity<Q0 _COMMA_ Q1>);
        EMPTY_SHELL_METHOD(emptyShell);
        IMPLEMENTS(Additive);
    }

    typename CrossRiskPropertySensitivity<Q0, Q1>::Deriv deriv() const {
        return typename CrossRiskPropertySensitivity<Q0, Q1>::Deriv(
             IResultsFunctionSP(),
             smartPtr<IRiskProperty<Q0> >(),
             smartPtr<IRiskProperty<Q1> >(),
             ICrossDerivativeSP(), 0.);
    }

public:

    GeneralCrossRiskPropertySensitivity(
            const string& outputName,
            IResultsFunctionConstSP derivand,
            smartConstPtr<IRiskProperty<Q0> > property0,
            smartConstPtr<IRiskProperty<Q1> > property1,
            ICrossDerivativeConstSP derivative,
            double sensitivityUnit,
            double shiftSize0,
            double shiftSize1):
        CrossRiskPropertySensitivity<Q0, Q1>(
            TYPE, shiftSize0, shiftSize1, outputName)
    {
        this->_derivand = derivand;
        this->_property0 = property0;
        this->_property1 = property1;
        this->_derivative = derivative;
        this->_sensitivityUnit = sensitivityUnit;
    }
};

template <class Q0, class Q1>
smartPtr<CrossRiskPropertySensitivity<Q0, Q1> > CrossRiskPropertySensitivity<Q0, Q1>::SP(
        const string& outputName,
        IResultsFunctionConstSP derivand,
        smartConstPtr<IRiskProperty<Q0> > property0,
        smartConstPtr<IRiskProperty<Q1> > property1,
        ICrossDerivativeConstSP derivative,
        double sensitivityUnit,
        double shiftSize0, double shiftSize1) {
    return smartPtr<CrossRiskPropertySensitivity<Q0, Q1> >(
        new GeneralCrossRiskPropertySensitivity<Q0, Q1>(
            outputName, derivand, property0, property1,
            derivative, sensitivityUnit, shiftSize0, shiftSize1));
}

template class RISKMGR_DLL CrossRiskPropertySensitivity<Void, Void>;
template class RISKMGR_DLL CrossRiskPropertySensitivity<ExpiryWindow, ExpiryWindow>;

template class GeneralCrossRiskPropertySensitivity<Void, Void>;
template class GeneralCrossRiskPropertySensitivity<ExpiryWindow, ExpiryWindow>;

template <>
CClassConstSP const GeneralCrossRiskPropertySensitivity<Void, Void>::TYPE = CClass::registerClassLoadMethod(
    "GeneralCrossRiskPropertySensitivity<Void, Void>", typeid(GeneralCrossRiskPropertySensitivity<Void, Void>), load);

template <>
CClassConstSP const GeneralCrossRiskPropertySensitivity<ExpiryWindow, ExpiryWindow>::TYPE = CClass::registerClassLoadMethod(
    "GeneralCrossRiskPropertySensitivity<ExpiryWindow, ExpiryWindow>", typeid(GeneralCrossRiskPropertySensitivity<ExpiryWindow, ExpiryWindow>), load);

DRLIB_END_NAMESPACE
