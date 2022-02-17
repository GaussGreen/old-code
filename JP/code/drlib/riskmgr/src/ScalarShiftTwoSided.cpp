#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/AbstractPropertyTweakHypothesis.hpp"
#include "edginc/HypotheticalQuantity.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/LazyRiskQuantityFactory.hpp"
#include "edginc/ScalarShiftTwoSided.hpp"
#include "edginc/Results.hpp"
#include "edginc/RiskMapping.hpp"

DRLIB_BEGIN_NAMESPACE
ScalarShiftTwoSided::~ScalarShiftTwoSided(){}

/** Overridden to perform two sided tweak */
void ScalarShiftTwoSided::calculate(TweakGroup*  tweakGroup,
                                    Results*     results){
    calculateTwoSidedDeriv(getSecondOrderSensOutputName(), tweakGroup, results);
}

double ScalarShiftTwoSided::divisor() const {
    if (_constantDivisor){
        return 1.0/_sensitivityUnit;
    }
    double shiftSize = getShiftSize();
    if (Maths::isZero(shiftSize)) {
        throw ModelException("GenericScalarTwoSidedShift::divisor",
                             "Shift size is zero");
    }
    // We use this as a hook to apply sensitivityUnit
    return shiftSize / _sensitivityUnit;
} 

/** Overridden to scale 1st and 2nd order derivatives */
void ScalarShiftTwoSided::scaleResult(Results*     results,
                                      double       scaleFactor) const{
    ScalarShift::scaleResult(results, scaleFactor); // do 1st order derive
    const string& secondOrderName = getSecondOrderSensOutputName();
    if (getPacketName() == getSensOutputName()) {
        results->scale(secondOrderName, scaleFactor);
    } else {
        results->scale(getPacketName(), secondOrderName, scaleFactor);
    }
}

/** Overridden to add 1st and 2nd order derivatives */
void ScalarShiftTwoSided::addResult(Results*           results,     // (M)
                                    const Results*     resultsToAdd,
                                    double             scaleFactor) const{
    ScalarShift::addResult(results, resultsToAdd,
                           scaleFactor); // do 1st order derive
    const string& secondOrderName = getSecondOrderSensOutputName();
    if (getPacketName() == getSensOutputName()) {
        results->add(secondOrderName, resultsToAdd, scaleFactor);
    } else {
        results->add(getPacketName(), secondOrderName, resultsToAdd,
                     scaleFactor);
    }
}

/** Returns the shift(s) which have been made for the current pricing
    call */
ScalarShiftArray ScalarShiftTwoSided::getComponentShifts() const{
    // need array of ScalarShiftConstSP to make this work - but this breaks
    // the array template - would need const array template...
    ScalarShiftTwoSided* shift = const_cast<ScalarShiftTwoSided*>(this);
    return ScalarShiftArray(1, ScalarShiftSP::attachToRef(shift));
}

/** Note ScalarShiftTwoSided is abstract. Create a scalar shift of
    type clazz, which uses outputName (eg VEGA_PARALLEL) to
    identify results and with given shiftSize */
ScalarShiftTwoSided::ScalarShiftTwoSided(
        CClassConstSP clazz,
        bool constantDivisor,
        double sensitivityUnit,
        const string& outputName,
        const IScalarRiskProperty* property,
        double coefficient):
    ScalarShift(clazz, outputName, coefficient),
    _constantDivisor(constantDivisor),
    _sensitivityUnit(sensitivityUnit),
    property(property),
    _currentHypothesisConst(
        dynamic_cast<const AbstractPropertyTweakHypothesis*>(
            property->axisFor(OutputNameConstSP(),
                              VoidSP())->hypothesis(1.).get()))
{
    if (!_currentHypothesisConst) {
        IRiskAxisConstSP ax = property->axisFor(OutputNameConstSP(),
                                                VoidSP());
        throw ModelException("ScalarShiftTwoSided::ScalarShiftTwoSided",
            "'property' must actually be a concrete RiskProperty. "
            "It's " + property->getClass()->getName() + " "
            "and it returned a " + ax->getClass()->getName() + " "
            "and that returned a " + ax->hypothesis(1.)->getClass()->getName() + " "
            "Sorry.  --William.");
    }

    AbstractPropertyTweakHypothesisSP::constCast(_currentHypothesisConst)->
        setAxisCoefficient(coefficient);
}

CClassConstSP ScalarShiftTwoSided::shiftInterface() const {
    return currentHypothesis()->shiftableInterface();
}

CClassConstSP ScalarShiftTwoSided::restorableShiftInterface() const {
    return currentHypothesis()->restorableInterface();
}

bool ScalarShiftTwoSided::nameMatches(const OutputName& name,
                                      IObjectConstSP obj) {
    return name.equals(currentHypothesis()->sensName(obj));
}

void ScalarShiftTwoSided::appendName(OutputNameArray& namesList,
                                     IObjectConstSP obj) {
    namesList.push_back(OutputNameSP(new OutputName(
        currentHypothesis()->sensName(obj))));
}

bool ScalarShiftTwoSided::shift(IObjectSP obj) {
    TweakOutcome outcome = currentHypothesis()->sensShift(obj);
    if (outcome.hasOldValue()) setInitialValue(outcome.oldValue());
    return outcome.tweakMembers();
}

void ScalarShiftTwoSided::restore(IObjectSP obj) {
    currentHypothesis()->sensRestore(obj);
}

void ScalarShiftTwoSided::setCurrentHypothesis(AtomicHypothesisArraySP history,
                                               AbstractPropertyTweakHypothesisConstSP current,
                                               double oldValue) {
    _currentHypothesisConst = current;
    _currentHypothesisMutable.reset(0);
    ScalarShift::setMarketDataName(current->getMarketDataName());
    ScalarShift::setShiftSize(current->axisCoefficient());
    if (oldValue != -666e66) setInitialValue(oldValue);
}

AbstractPropertyTweakHypothesisSP ScalarShiftTwoSided::currentHypothesisMutable() {
    if (!_currentHypothesisMutable) {
        _currentHypothesisMutable.reset(_currentHypothesisConst.clone());
        _currentHypothesisConst.reset(0);
    }
    return _currentHypothesisMutable;
}

AbstractPropertyTweakHypothesisConstSP ScalarShiftTwoSided::currentHypothesis() const {
    return !_currentHypothesisMutable ?
               _currentHypothesisConst : _currentHypothesisMutable;
}

void ScalarShiftTwoSided::setShiftSize(double shiftSize) {
    currentHypothesisMutable()->setAxisCoefficient(shiftSize);
    ScalarShift::setShiftSize(shiftSize);
}

void ScalarShiftTwoSided::setMarketDataName(OutputNameConstSP name) {
    currentHypothesisMutable()->setMarketDataName(name);
    ScalarShift::setMarketDataName(name);
}

OutputNameArrayConstSP ScalarShiftTwoSided::names(const IObject* world) const {

    if (hasOverrideNames()) {
        return OutputName::trim(overrideNames());
    }
    else {
        OutputNameArrayConstSP them =
                property->subjectNames(IObjectConstSP::attachToRef(world));
        const OutputNameArray &it = *them;
        for (int i = 0; i < it.size(); ++i)
            ASSERT(it[i].get() != NULL);

        return OutputName::trim(them);
    }

    //    return OutputName::trim(hasOverrideNames() ?
    //        overrideNames() :
    //        property->subjectNames(IObjectConstSP::attachToRef(world)));
}

// FIXME this is blatantly cut'n'pasted from ScalarShift and is a grotesque hack

ScalarShiftTwoSidedConstSP ScalarShiftTwoSided::alteredControl(MultiTweakGroupConstSP tweakGroup) const {
    try {
        ScalarShift*  alteredScalarShift = NULL; 

        SensControlSP alteredModelControl(
            tweakGroup->getModel()->AlterControl(this));

        SensControlSP alteredControl;

        if (tweakGroup->getInstruments()->size() == 1 /* FIXME yikes! */ ) {
            alteredControl.reset((*tweakGroup->getInstruments())[0]->AlterControl(
                                     getModel(), this));
        }

        if (!alteredControl){
            //instrument takes priority over control
            alteredControl = alteredModelControl; 
        }
        if ( !(!alteredControl) ) {
            alteredScalarShift = 
                dynamic_cast<ScalarShift*>(alteredControl.get());
            // set control/alorithm if not set
            if (!alteredScalarShift->control){
                alteredScalarShift->control = control;
            }
            if (!alteredScalarShift->algorithm){
                alteredScalarShift->algorithm = algorithm;
            }
        }

        // FIXME um ...

        return !alteredControl ? ScalarShiftTwoSidedConstSP::attachToRef(this) : ScalarShiftTwoSidedConstSP::dynamicCast(alteredControl);
    }
    catch (ModelException& e) {
        throw ModelException(e, "ScalarShiftTwoSided::alteredControl");
    }
}

static void setOriginatingSensitivity(NamedRiskQuantity& rq, const ICompatibilitySensitivity* sens) {
    HypotheticalQuantityArrayConstSP hqs = rq.riskQuantity->parameters();
    for (int h = 0; h < hqs->size(); ++h) {
        IHypothesisConstSP hyp = (*hqs)[h]->hypothesis();
        for (int a = 0; a < hyp->numAtomics(); ++a) {
            const AbstractPropertyTweakHypothesis* t = dynamic_cast<const AbstractPropertyTweakHypothesis*>(hyp->atomic(a).get());
            if (t) ((AbstractPropertyTweakHypothesis *)t)->_originatingSensitivity.reset(const_cast<ICompatibilitySensitivity*>(sens));
        }
    }
}

NamedRiskQuantityArraySP ScalarShiftTwoSided::riskQuantities(
        MultiTweakGroupConstSP world,
        RiskMappingConstSP riskMapping) const {

    NamedRiskQuantityArraySP rqs(new NamedRiskQuantityArray());

    OutputNameArrayConstSP tweakables = names(world.get());

    if (hasOverrideNames()){
        // manually cope with any explicit names that aren't present.
        // This should be a virtual method on SensControl/Sensitivity - 
        // also see comments in SensMgr::names(SensControl*, Results*)
        OutputNameArrayConstSP allNames(
            riskMapping->subjectNames(property, world));
        OutputNameArraySP extraNames(OutputName::difference(tweakables,
                                                            allNames));
        for (int i = 0; i < extraNames->size(); i++){
            rqs->push_back(NamedRiskQuantity::SP(
                RiskQuantity::notApplicable(),
                IResultsIdentifier::SP(getPacketName(), (*extraNames)[i])));
            rqs->push_back(NamedRiskQuantity::SP(
                RiskQuantity::notApplicable(),
                IResultsIdentifier::SP(getSecondOrderSensOutputName(),
                                       (*extraNames)[i])));
        }
    }

    string stepSizeName = getSensOutputName() + Results::SHIFT_SIZE_POSTFIX;

    for (int t = 0; t < tweakables->size(); ++t) {
        OutputNameConstSP prevmdn = getMarketDataName();
        ScalarShiftTwoSided* mthis = const_cast<ScalarShiftTwoSided*>(this);
        mthis->setMarketDataName((*tweakables)[t]); // FIXME should undo
        ScalarShiftTwoSidedConstSP altered;
        try {
            altered = alteredControl(world);
        }
        catch (...) {
            mthis->setMarketDataName(prevmdn);
            throw;
        }
        mthis->setMarketDataName(prevmdn);

        IRiskAxisConstSP axis = altered->property->axisFor((*tweakables)[t], VoidSP());

        rqs->push_back(NamedRiskQuantity::SP(
            IScalarDerivative::twoSided()->riskQuantity(
               IResultsFunction::price(), axis, altered->getShiftSize()),
            IResultsIdentifier::SP(altered->getSensOutputName(), (*tweakables)[t]),
            _sensitivityUnit));

        setOriginatingSensitivity(*rqs->back(), altered.get());

        rqs->push_back(NamedRiskQuantity::SP(
            IScalarDerivative::second()->riskQuantity(
               IResultsFunction::price(), axis, altered->getShiftSize()),
            IResultsIdentifier::SP(altered->getSecondOrderSensOutputName(),
                                   (*tweakables)[t]),
            _sensitivityUnit * _sensitivityUnit));

        setOriginatingSensitivity(*rqs->back(), altered.get());

        rqs->push_back(NamedRiskQuantity::SP(
            RiskQuantity::constant(altered->getShiftSize()),
            IResultsIdentifier::SP(stepSizeName, (*tweakables)[t])));
    }

    if (rqs->empty()) {
        rqs->push_back(NamedRiskQuantity::SP(
            RiskQuantity::notApplicable(),
            IResultsIdentifier::SP(getPacketName(), OutputNameConstSP())));
        rqs->push_back(NamedRiskQuantity::SP(
            RiskQuantity::notApplicable(),
            IResultsIdentifier::SP(getSecondOrderSensOutputName(), OutputNameConstSP())));
    }

    return rqs;
}

LazyRiskQuantityFactoryArraySP ScalarShiftTwoSided::lazies(MultiTweakGroupConstSP world) const {
    return LazyRiskQuantityFactoryArraySP();
}
void ScalarShiftTwoSided::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ScalarShiftTwoSided, clazz);
    SUPERCLASS(ScalarShift);
    IMPLEMENTS(Additive);
    IMPLEMENTS(ICompatibilitySensitivity);
    IMPLEMENTS(IRiskQuantityFactory);
    FIELD(_constantDivisor, "_constantDivisor");
    FIELD_MAKE_TRANSIENT(_constantDivisor);
    FIELD(_sensitivityUnit, "_sensitivityUnit");
    FIELD_MAKE_TRANSIENT(_sensitivityUnit);
    FIELD(property, "property");
    FIELD_MAKE_TRANSIENT(property);
    FIELD(_currentHypothesisConst, "_currentHypothesisConst");
    FIELD_MAKE_TRANSIENT(_currentHypothesisConst);
    FIELD(_currentHypothesisMutable, "_currentHypothesisMutable");
    FIELD_MAKE_TRANSIENT(_currentHypothesisMutable);
}

CClassConstSP const ScalarShiftTwoSided::TYPE = CClass::registerClassLoadMethod(
    "ScalarShiftTwoSided", typeid(ScalarShiftTwoSided), load);

DRLIB_END_NAMESPACE
