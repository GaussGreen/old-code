/**
 * @file RiskQuantityFactorySensitivity.cpp
 */

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/HypotheticalQuantity.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/LazyRiskQuantityFactory.hpp"
#include "edginc/AbstractPropertyTweakHypothesis.hpp"
#include "edginc/RiskQuantityFactorySensitivity.hpp"
#include "edginc/RiskQuantityEvaluator.hpp"
#include "edginc/Results.hpp"
#include "edginc/RiskMapping.hpp"

DRLIB_BEGIN_NAMESPACE

RiskQuantityFactorySensitivity::RiskQuantityFactorySensitivity(
        CClassConstSP type,
        const string& outputName,
        const string& outputName2):
    Sensitivity(type),
    outputNames(new StringArray(1, outputName))
{
    if (outputName2 != "") {
        StringArraySP::constCast(outputNames)->push_back(outputName2);
    }
}

RiskQuantityFactorySensitivity::~RiskQuantityFactorySensitivity() {}

void RiskQuantityFactorySensitivity::setCurrentHypothesis(
        AtomicHypothesisArraySP history,
        AbstractPropertyTweakHypothesisConstSP current,
        double oldValue) {
    marketDataName = current->getMarketDataName();
}

void RiskQuantityFactorySensitivity::calculate(TweakGroup* tweakGroup, CResults* results) {
    try {
        // Believe it or not this is a VERY thin wrapper round RiskQuantityEvaluator::storeResults()

        RiskQuantityEvaluator().storeResults(
            IRiskQuantityFactoryArray::SP(1, IRiskQuantityFactorySP::attachToRef(this)),
            MultiTweakGroup::SP(
                IInstrumentCollection::singleton(tweakGroup->getInstrument()),
                tweakGroup->getModelSP()),
            CControlSP::attachToRef(control),
            CResultsArray::SP(1, CResultsSP::attachToRef(results)));
    }
    catch (exception& e) {
        throw ModelException(&e, "RiskQuantityFactorySensitivity::calculate");
    }
}

const string& RiskQuantityFactorySensitivity::getSensOutputName() const {
    return (*outputNames)[0];
}

void RiskQuantityFactorySensitivity::scaleResult(Results* results,
                                                 double scaleFactor) const {
    Sensitivity::scaleResult(results, scaleFactor);

    for (int d = 1; d < outputNames->size(); ++d) {
        if (getPacketName() == getSensOutputName()) {
            results->scale((*outputNames)[d], scaleFactor);
        }
        else {
            results->scale(getPacketName(), (*outputNames)[d], scaleFactor);
        }
    }
}

void RiskQuantityFactorySensitivity::addResult(Results* results,
                                               const Results* resultsToAdd,
                                               double scaleFactor) const {
    Sensitivity::addResult(results, resultsToAdd,
                           scaleFactor); // do 1st order derive

    for (int d = 1; d < outputNames->size(); ++d) {
        if (getPacketName() == getSensOutputName()) {
            results->add((*outputNames)[d], resultsToAdd, scaleFactor);
        } else {
            results->add(getPacketName(), (*outputNames)[d], resultsToAdd,
                         scaleFactor);
        }
    }
}

/**/void setOriginatingSensitivity(NamedRiskQuantityArraySP rqs, const ICompatibilitySensitivity* sens) {
/**/    for (int r = 0; r < rqs->size(); ++r) {
/**/        HypotheticalQuantityArrayConstSP hqs = (*rqs)[r]->riskQuantity->parameters();
/**/        for (int h = 0; h < hqs->size(); ++h) {
/**/            IHypothesisConstSP hyp = (*hqs)[h]->hypothesis();
/**/            for (int a = 0; a < hyp->numAtomics(); ++a) {
/**/                const AbstractPropertyTweakHypothesis* t = dynamic_cast<const AbstractPropertyTweakHypothesis*>(hyp->atomic(a).get());
/**/                if (t) ((AbstractPropertyTweakHypothesis *)t)->_originatingSensitivity.reset(const_cast<ICompatibilitySensitivity*>(sens));
/**/            }
/**/        }
/**/    }
/**/}

NamedRiskQuantityArraySP RiskQuantityFactorySensitivity::riskQuantities(
        MultiTweakGroupConstSP world,
        RiskMappingConstSP riskMapping) const {

    try {
        NamedRiskQuantityArraySP rqs(nameRiskQuantities(world, riskMapping));

        if (rqs->empty()) {
            for (int d = 0; d < outputNames->size(); ++d) {
                TRACE_BLOCK("Creating NotApplicable NamedRiskQuantity("
                            << "packet="
                            << (d == 0 ? getPacketName() : (*outputNames)[d])
                            << ")");

                rqs->push_back(NamedRiskQuantity::SP(
                    RiskQuantity::notApplicable(),
                    // FIXME probably wrong in some sense
                    IResultsIdentifier::SP(
                        d == 0 ? getPacketName() : (*outputNames)[d],
                        OutputNameConstSP())));
            }
        }

/**/    setOriginatingSensitivity(rqs, this);

        return rqs;
    }
    catch (exception& e) {
        throw ModelException(
            e, "RiskQuantityFactorySensitivity::riskQuantities()");
    }
}

LazyRiskQuantityFactoryArraySP RiskQuantityFactorySensitivity::lazies(MultiTweakGroupConstSP world) const {
    return LazyRiskQuantityFactoryArraySP();
}

void RiskQuantityFactorySensitivity::load(CClassSP& clazz) {
    REGISTER(RiskQuantityFactorySensitivity, clazz);
    SUPERCLASS(Sensitivity);
    IMPLEMENTS(IRiskQuantityFactory);
    IMPLEMENTS(ICompatibilitySensitivity);
    FIELD(outputNames, "outputNames");
    FIELD_MAKE_TRANSIENT(outputNames);
}

class ArrayRiskQuantityFactorySensitivity:
    public RiskQuantityFactorySensitivity {

public:

    ArrayRiskQuantityFactorySensitivity(
            const string& packetName,
            NamedRiskQuantityArraySP rqs):
        RiskQuantityFactorySensitivity(TYPE, packetName),
        rqs(rqs) 
    {}

    bool discreteShift() const {
        return true;
    }

    static IObject* defaultOne() {
        return new ArrayRiskQuantityFactorySensitivity(
            "", NamedRiskQuantityArraySP());
    }

    static void load(CClassSP& clazz) {
        REGISTER(ArrayRiskQuantityFactorySensitivity, clazz);
        SUPERCLASS(RiskQuantityFactorySensitivity);
        EMPTY_SHELL_METHOD(defaultOne);
        FIELD(rqs, "rqs");
    }

    static CClassConstSP const TYPE;

    NamedRiskQuantityArraySP rqs;

    NamedRiskQuantityArraySP nameRiskQuantities(
            MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping)
                const {
        return rqs;
    }
};

CClassConstSP const ArrayRiskQuantityFactorySensitivity::TYPE = CClass::registerClassLoadMethod(
    "ArrayRiskQuantityFactorySensitivity", typeid(ArrayRiskQuantityFactorySensitivity), load);

RiskQuantityFactorySensitivitySP RiskQuantityFactorySensitivity::singleton(
        const string& packetName,
        NamedRiskQuantityConstSP rq) {
    return RiskQuantityFactorySensitivitySP(
        new ArrayRiskQuantityFactorySensitivity(
            packetName,
            NamedRiskQuantityArray::SP(
                1,
                NamedRiskQuantitySP::dynamicCast(IObjectSP(rq->clone())))));
}

CClassConstSP const RiskQuantityFactorySensitivity::TYPE = CClass::registerClassLoadMethod(
    "RiskQuantityFactorySensitivity", typeid(RiskQuantityFactorySensitivity), load);

DRLIB_END_NAMESPACE
