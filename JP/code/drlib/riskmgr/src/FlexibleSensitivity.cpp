/**
 * @file FlexibleSensitivity.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/OutputName.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/RiskMapping.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/LazyRiskQuantityFactory.hpp"
#include "edginc/FieldSensitivityDefinition.hpp"
#include "edginc/RiskQuantityFactorySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * A greek defined with respect to the value of a named field on some class
 *
 * This class is specific to the "flexible scenarios and greeks" facility: it's
 * a client-friendly vehicle for specifying a scalar or pointwise greek with
 * respect to a the value of a field designated through the QLib reflection
 * system.  For instance, a client could use it to construct a greek gave
 * similar results to ParSpreadRhoPointwise without adding anything to QLib
 * (specifying the class "ParSpreadCurve" and the field "parSpreads" by name).
 *
 * If you want to construct greeks on the fly inside QLib, you'll find it much
 * more convenient to use ScalarRiskPropertySensitivity::SP() and
 * PerNameRiskPropertySensitivity<ExpiryWindow>::SP() directly, together with
 * fieldRiskProperty::scalar(), fieldRiskProperty::pointwise() etc.
 *
 *
 * <H3>Background</H3>
 *
 * Historically we have supported new greeks by hard-coding them in C++.  This
 * gives a lot of flexibility, but often (e.g. ParSpreadRhoParallel,
 * LegalBasisAdditivePointwise) the implementation works just by tweaking the
 * values of a concrete field on a particular class of market object.
 * FlexibleSensitivity is a means of reducing time-to-market and cost
 * in these cases.
 *
 * ScalarFlexiblePerturbation and PerEntryFlexiblePerturbation provide a similar
 * facility for generic Scenario's.
 *
 * For an externally-oriented overview of the framework, see "Flexible
 * scenarios and greeks" in the QLib doc db:
 *
 * -     Notes://PPUSMC006/8525703B0051D832/55B689E36191AD7285256DFD00470F6B/272C882E688B88B3852570CB00403228
 *
 *
 * <H3>Architecture</H3>
 *
 * FlexibleSensitivity can be used day-to-day in RiskMgrInterface alongside
 * predefined sensitivities like Delta; similarly ScalarFlexiblePerturbation
 * and PerEntryFlexiblePerturbation can be used in ScenarioInterface alongside
 * predefined scenarios like SpotShift.  However, the only parameter which can
 * be varied is essentially the shift size (or scenario level).  The details of
 * which fields are to be perturbed and which derivatives calculated are
 * defined in two further types IFieldRiskPropertyDefinition and
 * FieldSensitivityDefinition.  Additions and alterations to the "repertore" of
 * available definitions can be made quickly in Pyramid, but require signoff
 * from QR, and in practice will almost always be performed by QR (since
 * knowledge of QLib internals is required).
 *
 * Inside QLib, IFieldRiskPropertyDefinition is a factory for three kinds of
 * generic IRiskProperty:
 *
 *    -  fieldRiskProperty::scalar() represents a scalar-valued field, or
 *       an array- or matrix-valued field all of whose elements are
 *       tweaked at once;
 *
 *    -  fieldRiskProperty::parallel() represents an array- or matrix-valued
 *       field tweaked at a given set of expiries;
 *
 *    -  fieldRiskProperty::pointwise() represents an array- or matrix-valued
 *       field tweaked for each available expiry.
 *
 * The first two of these are used respectively by ScalarFlexiblePerturbation
 * and PerEntryFlexiblePerturbation to obtain IHypothesis's implementing the
 * requested scenario.
 *
 * FieldSensitivityDefinition is a factory for sensitivitys which obtains one
 * of the three IRiskProperty's above from a IFieldRiskPropertyDefinition,
 * combines it with an IScalarDerivative and IResultsFunction, and feeds them
 * into a standard PerNameRiskPropertySensitivity (the same class from which
 * hard-coded greeks such as VegaParallel and ParSpreadRhoPointwise derive).
 *
 *
 * <H3>See also</H3>
 *
 * For background information on the "declarative" sensitivities framework
 * in which the flexible scenarios and greeks are implemented, see
 * IRiskQuantityFactory.
 */

class FlexibleSensitivity: public Sensitivity,
                           public virtual IRiskQuantityFactory,
                           public virtual Additive {

    FlexibleSensitivity(): Sensitivity(TYPE) {}

    static IObject* defaultOne() {
        return new FlexibleSensitivity();
    }

    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FlexibleSensitivity, clazz);
        clazz->setPublic();
        SUPERCLASS(Sensitivity);
        IMPLEMENTS(IRiskQuantityFactory);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultOne);

        FIELD(definition, "definition");
        FIELD(shiftSize, "shiftSize");
        FIELD_MAKE_OPTIONAL(shiftSize);
        FIELD(packetName, "packetName");
        FIELD_MAKE_OPTIONAL(packetName);
        FIELD(packetName2, "packetName2");
        FIELD_MAKE_OPTIONAL(packetName2);

        FIELD(_underlying, "underlying");
        FIELD_MAKE_TRANSIENT(_underlying);
    }

public:

    static CClassConstSP const TYPE;

private:

    /**
     * Details of derivand, derivative, field w.r.t. which it's taken etc.
     */

    FieldSensitivityDefinitionConstSP definition;

    /**
     * Shift size (overriding that in 'definition' if present)
     */

    CDoubleConstSP shiftSize;

    /**
     * Name of packet where we store numbers in Results (overriding that in
     * 'definition' if present)
     */

    CStringConstSP packetName;

    /**
     * Name of packet where we store numbers in Results (overriding that in
     * 'definition' if present) for any second derivative calculated
     * by the sensitivity
     */

    CStringConstSP packetName2;

    // Structurally it's a "delegator", passing all method calls to an
    // underlying implementation obtained from
    // FieldSensitivityDefinition::sensitivity()

    RiskQuantityFactorySensitivitySP _underlying;

public:

    ~FlexibleSensitivity() {}

    void validatePop2Object() {
        try {
            double s = !shiftSize ? 0 : shiftSize->doubleValue();

            _underlying = definition->sensitivity(
                toTweak,
                !shiftSize ? NULL : &s,
                !packetName ? NULL : &packetName->stringValue(),
                !packetName2 ? NULL : &packetName2->stringValue());
        }
        catch (exception& e) {
            throw ModelException(
                e, "FlexibleSensitivity::validatePop2Object()");
        }
    }

    /**
     * @name Sensitivity implementation
     */

    //@{

    virtual const string& getSensOutputName() const {
        return _underlying->getSensOutputName();
    }

    virtual bool discreteShift() const {
        return _underlying->discreteShift();
    }

    virtual const string& getPacketName() const {
        return _underlying->getPacketName();
    }

    virtual void scaleResult(Results*     results,     // (M)
                             double       scaleFactor) const {
        _underlying->scaleResult(results, scaleFactor);
    }

    virtual void addResult(Results*           results,     // (M)
                           const Results*     resultsToAdd,
                           double             scaleFactor) const {
        _underlying->addResult(results, resultsToAdd, scaleFactor);
    }

    virtual double getSensPrice(Results*     results,
                                CInstrument* inst,
                                IModel*      model,
                                Control*     control) {
        return _underlying->getSensPrice(results, inst, model, control);
    }
    
    virtual void calculateSens(IModel*          algorithm,
                               CInstrument*     instrument,
                               Control*         control,
                               Results*         results) {
        this->algorithm = algorithm;
        this->control = control;
        try {
            _underlying->calculateSens(algorithm, instrument, control, results);
        }
        catch (...) {
            this->algorithm = 0;
            this->control = 0;
            throw;
        }
        this->algorithm = 0;
        this->control = 0;
    }

    virtual void getMarket(const IModel* model, const MarketData* market) {
        _underlying->getMarket(model, market);
    }

    virtual bool hasOverrideNames() const {
        return _underlying->hasOverrideNames();
    }

    virtual OutputNameArrayConstSP overrideNames() const {
        return _underlying->overrideNames();
    }

    virtual void storeOverrideNames(OutputNameArraySP names) {
        _underlying->storeOverrideNames(names);
    }

    virtual Sensitivity * spawn(IModel* model) const {
        return _underlying->spawn(model);
    }

    virtual void removeOverrideNames(const string& packetName,
                                     const Results* results) {
        _underlying->removeOverrideNames(packetName, results);
    }

    virtual void setAlgorithm(IModel* a) {
        Sensitivity::setAlgorithm(a);
        _underlying->setAlgorithm(a);
    }

    virtual void setControl(Control* c) {
        Sensitivity::setControl(c);
        _underlying->setControl(c);
    }

    virtual void calculate(TweakGroup*      tweakGroup,
                           Results*         results) {
        _underlying->calculate(tweakGroup, results);
    }

    //@}

    /**
     * @name IRiskQuantityFactory implementation
     */

    //@{

    NamedRiskQuantityArraySP riskQuantities(
            MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping) const {
        return _underlying->riskQuantities(world, riskMapping);
    }

    LazyRiskQuantityFactoryArraySP lazies(
            MultiTweakGroupConstSP world) const {
        return _underlying->lazies(world);
    }

    //@}
};

CClassConstSP const FlexibleSensitivity::TYPE = CClass::registerClassLoadMethod(
    "FlexibleSensitivity", typeid(FlexibleSensitivity), load);

bool FlexibleSensitivityLinkIn() {
    return FlexibleSensitivity::TYPE != NULL;
}

DRLIB_END_NAMESPACE
