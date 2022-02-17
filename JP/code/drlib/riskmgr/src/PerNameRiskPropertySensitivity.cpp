/**
 * @file PerNameRiskPropertySensitivity.cpp
 */

#include "edginc/config.hpp"
#define QLIB_PERNAMERISKPROPERTYSENSITIVITY_CPP
#include "edginc/TRACE.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/BoxedInt.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/ExpiryPair.hpp"
#include "edginc/ExpiryAndStrike.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/RiskMapping.hpp"

DRLIB_BEGIN_NAMESPACE


template <class QUALIFIER>
PerNameRiskPropertySensitivity<QUALIFIER>::PerNameRiskPropertySensitivity(
        CClassConstSP type,
        double shiftSize,
        const string& name,
        const string& name2):
    RiskPropertySensitivity<QUALIFIER>(type, shiftSize, name, name2)
{}

template <class QUALIFIER>
PerNameRiskPropertySensitivity<QUALIFIER>::~PerNameRiskPropertySensitivity() {}

template <class QUALIFIER>
NamedRiskQuantityArraySP PerNameRiskPropertySensitivity<QUALIFIER>::nameRiskQuantities(
        MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping) const {

    TRACE_METHOD;

    OutputNameArrayConstSP inWorld =
        OutputName::trim(riskMapping->subjectNames(this->property(), world));

    OutputNameArrayConstSP tweakables = this->hasOverrideNames() ?
        OutputName::intersection(inWorld, this->overrideNames()) :
        inWorld;

    NamedRiskQuantityArraySP rqs(new NamedRiskQuantityArray());

    for (int t = 0; t < tweakables->size(); ++t) {
        TRACE_BLOCK("Market name " << *(*tweakables)[t]);
        try {
            DECLARE(QUALIFIER)
            QUALIFIERArrayConstSP quals = this->hasPredefinedQualifiers ?
                this->predefinedQualifiers :
                riskMapping->subjectQualifiers(this->property(), (*tweakables)[t], world);

/**/        this->marketDataName = (*tweakables)[t]; // so that CDSPS::endDate() works

            BoolArrayConstSP isSens = this->property()->mayHaveEffect(
                                          world, this->marketDataName, quals, this);
            ASSERT(!!isSens && isSens->size() == quals->size());

            DoubleArrayConstSP coeffs = this->property()->adaptiveCoefficients(
                world, this->marketDataName, quals, this->shiftSize);
            ASSERT(!!coeffs && coeffs->size() == quals->size());

            for (int q = 0; q < quals->size(); ++q) {
                IRiskAxisConstSP axis = this->property()->axisFor(this->marketDataName,
                                                            (*quals)[q]);

                for (int d = 0; d < this->derivatives()->size(); ++d) {
                    rqs->push_back(NamedRiskQuantity::SP(
                        (*isSens)[q] ? (*(this->derivatives()))[d]->riskQuantity(
                                           this->derivand(), axis, (*coeffs)[q]) :
                                       RiskQuantity::constant(0),
                        IResultsIdentifier::SP(
                            d == 0 ? this->getSensOutputName() :
                                     (*this->outputNames)[d],
                            this->marketDataName,
                            (*quals)[q]),
                        pow((*(this->sensitivityUnits()))[d],
                            (*(this->derivatives()))[d]->order())));

                    if ((*isSens)[q]) {
                        TRACE("Including request to compute " <<
                              *rqs->back()->resultsName <<
                              " with shift size " << (*coeffs)[q]);
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
            for (int u = 0; u < this->derivatives()->size(); ++u) {
                rqs->push_back(NamedRiskQuantity::SP(
                                   RiskQuantity::untweakable(e),
                                   IResultsIdentifier::SP((*this->outputNames)[u],
                                                          (*tweakables)[t])));
            }

            TRACE("Including request to mark " <<
                  *rqs->back()->resultsName << " as Untweakable");
        }
    }

/**/this->marketDataName.reset();

    if (this->hasOverrideNames()) {
        OutputNameArrayConstSP na = OutputName::trim(OutputName::difference(
            this->overrideNames(), tweakables));
        for (int n = 0; n < na->size(); ++n) {
            for (int u = 0; u < this->outputNames->size(); ++u) {
                rqs->push_back(NamedRiskQuantity::SP(
                    RiskQuantity::notApplicable(),
                    // FIXME probably wrong in some sense
                    IResultsIdentifier::SP(
                        u == 0 ? this->getPacketName() : (*this->outputNames)[u],
                        (*na)[n])));
            }
        }
    }

    return rqs;
}

template <class QUALIFIER>
OutputNameConstSP PerNameRiskPropertySensitivity<QUALIFIER>::getMarketDataName() const {
    return this->marketDataName;
}

template <class QUALIFIER>
OutputNameArrayConstSP PerNameRiskPropertySensitivity<QUALIFIER>::allNames(
        const IObject* object) const {
    return this->property()->subjectNames(IObjectConstSP::attachToRef(object));
}

template <class QUALIFIER>
CClassConstSP PerNameRiskPropertySensitivity<QUALIFIER>::shiftInterface() const {
    return this->property()->subjectInterface();
}

template <class QUALIFIER>
class GeneralPerNameRiskPropertySensitivity:
        public PerNameRiskPropertySensitivity<QUALIFIER>,
        public virtual Additive {

    static CClassConstSP const TYPE;

private:

    GeneralPerNameRiskPropertySensitivity():
        PerNameRiskPropertySensitivity<QUALIFIER>(TYPE, 0., "", "") {}

    static IObject* emptyShell() {
        return new GeneralPerNameRiskPropertySensitivity();
    }        

    static void load(CClassSP& clazz) {
        REGISTER(GeneralPerNameRiskPropertySensitivity, clazz);
        SUPERCLASS(PerNameRiskPropertySensitivity<QUALIFIER>);
        EMPTY_SHELL_METHOD(emptyShell);
        IMPLEMENTS(Additive);
    }

    typename RiskPropertySensitivity<QUALIFIER>::Deriv deriv() const {
        return typename RiskPropertySensitivity<QUALIFIER>::Deriv
            (IResultsFunctionSP(),
             smartPtr<IRiskProperty<QUALIFIER> >(),
             IScalarDerivativeSP(),
             0.);
    }

public:

    GeneralPerNameRiskPropertySensitivity(
            const string& outputName1,
            const string& outputName2,
            IResultsFunctionConstSP derivand,
            smartConstPtr<IRiskProperty<QUALIFIER> > property,
            IScalarDerivativeConstSP derivative1,
            IScalarDerivativeConstSP derivative2,
            double sensitivityUnit1,
            double sensitivityUnit2,
            double shiftSize):
        PerNameRiskPropertySensitivity<QUALIFIER>(
            TYPE, shiftSize, outputName1, !derivative2 ? string() : outputName2)
    {
        this->_derivand = derivand;
        this->_property = property;
        this->_derivatives.reset(new IScalarDerivativeArray(
            1, IScalarDerivativeSP::constCast(derivative1)));
        this->_sensitivityUnits.reset(new DoubleArray(1, sensitivityUnit1));

        if (!!derivative2) {
            IScalarDerivativeArraySP::constCast(this->_derivatives)->push_back(
                IScalarDerivativeSP::constCast(derivative2));
            DoubleArraySP::constCast(this->_sensitivityUnits)->push_back(
                sensitivityUnit2);
        }
    }
};

template <class QUALIFIER>
smartPtr<PerNameRiskPropertySensitivity<QUALIFIER> > PerNameRiskPropertySensitivity<QUALIFIER>::SP(
        const string& outputName,
        IResultsFunctionConstSP derivand,
        smartConstPtr<IRiskProperty<QUALIFIER> > property,
        IScalarDerivativeConstSP derivative,
        double sensitivityUnit,
        double shiftSize) {
    return smartPtr<PerNameRiskPropertySensitivity<QUALIFIER> >(
        new GeneralPerNameRiskPropertySensitivity<QUALIFIER>(
             outputName, "", derivand, property,
             derivative, IScalarDerivativeConstSP(),
             sensitivityUnit, 0, shiftSize));
}

template <class QUALIFIER>
smartPtr<PerNameRiskPropertySensitivity<QUALIFIER> > PerNameRiskPropertySensitivity<QUALIFIER>::SP(
        const string& outputName1,
        const string& outputName2,
        IResultsFunctionConstSP derivand,
        smartConstPtr<IRiskProperty<QUALIFIER> > property,
        IScalarDerivativeConstSP derivative1,
        IScalarDerivativeConstSP derivative2,
        double sensitivityUnit1,
        double sensitivityUnit2,
        double shiftSize) {
    return smartPtr<PerNameRiskPropertySensitivity<QUALIFIER> >(
        new GeneralPerNameRiskPropertySensitivity<QUALIFIER>(
            outputName1, outputName2, derivand, property,
            derivative1, derivative2,
            sensitivityUnit1, sensitivityUnit2, shiftSize));
}

template <class QUALIFIER>
smartPtr<PerNameRiskPropertySensitivity<QUALIFIER> > PerNameRiskPropertySensitivity<QUALIFIER>::withQualifier(
        QualifierConstSP qualifier) const {
    smartPtr<PerNameRiskPropertySensitivity<QUALIFIER> > it(copy(this));

    it->predefinedQualifiers.reset(
        new QualifierArray(1, QualifierSP::constCast(qualifier)));

    it->hasPredefinedQualifiers = true;
    return it;
}

template <class QUALIFIER>
void PerNameRiskPropertySensitivity<QUALIFIER>::load(CClassSP& clazz) {
    REGISTER(PerNameRiskPropertySensitivity, clazz);
    SUPERCLASS(RiskPropertySensitivity<QUALIFIER>);
    IMPLEMENTS(IPerNameSensitivity);
}

template <>
CClassConstSP const PerNameRiskPropertySensitivity<Void>::TYPE = CClass::registerClassLoadMethod(
    "PerNameRiskPropertySensitivity<Void>", typeid(PerNameRiskPropertySensitivity<Void>), load);

template <>
CClassConstSP const PerNameRiskPropertySensitivity<ExpiryWindow>::TYPE = CClass::registerClassLoadMethod(
    "PerNameRiskPropertySensitivity<ExpiryWindow>", typeid(PerNameRiskPropertySensitivity<ExpiryWindow>), load);

template <>
CClassConstSP const PerNameRiskPropertySensitivity<ExpiryPair>::TYPE = CClass::registerClassLoadMethod(
    "PerNameRiskPropertySensitivity<ExpiryPair>", typeid(PerNameRiskPropertySensitivity<ExpiryPair>), load);

template <>
CClassConstSP const PerNameRiskPropertySensitivity<ExpiryAndStrike>::TYPE = CClass::registerClassLoadMethod(
    "PerNameRiskPropertySensitivity<ExpiryAndStrike>", typeid(PerNameRiskPropertySensitivity<ExpiryAndStrike>), load);

template <>
CClassConstSP const PerNameRiskPropertySensitivity<BoxedInt>::TYPE = CClass::registerClassLoadMethod(
    "PerNameRiskPropertySensitivity<BoxedInt>", typeid(PerNameRiskPropertySensitivity<BoxedInt>), load);

// Must use normal template instantiation since we want these instantiated
// for both debug and optimised
template class RISKMGR_DLL PerNameRiskPropertySensitivity<ExpiryWindow>;
template class RISKMGR_DLL PerNameRiskPropertySensitivity<BoxedInt>;
template class RISKMGR_DLL PerNameRiskPropertySensitivity<ExpiryPair>;
template class RISKMGR_DLL PerNameRiskPropertySensitivity<ExpiryAndStrike>;
template class RISKMGR_DLL PerNameRiskPropertySensitivity<Void>;


template <>
CClassConstSP const GeneralPerNameRiskPropertySensitivity<Void>::TYPE = CClass::registerClassLoadMethod(
    "GeneralPerNameRiskPropertySensitivity<Void>", typeid(GeneralPerNameRiskPropertySensitivity<Void>), load);

template <>
CClassConstSP const GeneralPerNameRiskPropertySensitivity<ExpiryWindow>::TYPE = CClass::registerClassLoadMethod(
    "GeneralPerNameRiskPropertySensitivity<ExpiryWindow>", typeid(GeneralPerNameRiskPropertySensitivity<ExpiryWindow>), load);

template <>
CClassConstSP const GeneralPerNameRiskPropertySensitivity<ExpiryPair>::TYPE = CClass::registerClassLoadMethod(
    "GeneralPerNameRiskPropertySensitivity<ExpiryPair>", typeid(GeneralPerNameRiskPropertySensitivity<ExpiryPair>), load);

template <>
CClassConstSP const GeneralPerNameRiskPropertySensitivity<ExpiryAndStrike>::TYPE = CClass::registerClassLoadMethod(
    "GeneralPerNameRiskPropertySensitivity<ExpiryAndStrike>", typeid(GeneralPerNameRiskPropertySensitivity<ExpiryAndStrike>), load);

template <>
CClassConstSP const GeneralPerNameRiskPropertySensitivity<BoxedInt>::TYPE = CClass::registerClassLoadMethod(
    "GeneralPerNameRiskPropertySensitivity<BoxedInt>", typeid(GeneralPerNameRiskPropertySensitivity<BoxedInt>), load);

template class GeneralPerNameRiskPropertySensitivity<Void>;
template class GeneralPerNameRiskPropertySensitivity<ExpiryWindow>;
template class GeneralPerNameRiskPropertySensitivity<ExpiryPair>;
template class GeneralPerNameRiskPropertySensitivity<ExpiryAndStrike>;
template class GeneralPerNameRiskPropertySensitivity<BoxedInt>;

DRLIB_END_NAMESPACE
