/**
 * @file ScalarFlexiblePerturbation.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Perturbation.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/IHypothesis.hpp"
#include "edginc/IFieldRiskPropertyDefinition.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * An IPerturbation which specifies a change to the value of a named field on
 * some class
 *
 * A client-friendly vehicle for specifying a ScenarioShift under which a field
 * designated through the QLib reflection system has its value adjusted.  The
 * field can be of type double, DoubleArray or CDoubleMatrix (in the latter
 * cases all its entries will be adjusted), it doesn't have to be directly
 * attached to the "anchor" class, and there are several operations available
 * for shifting or setting it.  See IFieldRiskPropertyDefinition.
 *
 * If you want to construct scenarios on the fly inside QLib, you'll find it
 * much more convenient to use ScalarFieldRiskProperty directly.
 * 
 * See FlexibleSensitivity for an overview of "flexible scenarios and
 * greeks", and IRiskQuantityFactory for background info on the
 * "declarative sensitivities" framework in which they're implemented.
 */

class ScalarFlexiblePerturbation: public CObject,
                                  public virtual IPerturbation {

    ScalarFlexiblePerturbation():
        CObject(TYPE),
        _coeff(1.)
    {}

    static IObject* defaultOne() {
        return new ScalarFlexiblePerturbation();
    }

    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ScalarFlexiblePerturbation, clazz);
        clazz->setPublic();
        SUPERCLASS(CObject);
        IMPLEMENTS(IPerturbation);
        EMPTY_SHELL_METHOD(defaultOne);

        FIELD(property, "property");
        FIELD(argument, "argument");

        FIELD(_property, "_property");
        FIELD_MAKE_TRANSIENT(_property);
        FIELD(_coeff, "_coeff");
        FIELD_MAKE_TRANSIENT(_coeff);
    }

    void validatePop2Object() {
        try {
            // FIXME Oh dear, a bit of a hackette --- explain ...

            IObjectConstSP argOverCoeff;
            
            if (!!argument && CDouble::TYPE->isInstance(argument)) {
                _coeff = CDoubleConstSP::dynamicCast(argument)->doubleValue();
                argOverCoeff = CDouble::SP(1.);
            }
            else {
                argOverCoeff = argument;
            }

            _property = property->scalarProperty(argOverCoeff);
        }
        catch (exception& e) {
            throw ModelException(
                e, "ScalarFlexiblePerturbation::validatePop2Object()");
        }
    }

public:

    static CClassConstSP const TYPE;

private:

    /**
     * Details of the class, its field and how we tweak it
     *
     * IFieldRiskPropertyDefinition is just a factory for FieldRiskProperty's,
     * specified in a client-friendly way.
     */

    IFieldRiskPropertyDefinitionConstSP property; // $required

    /**
     * Shift size passed to TweakFunction::additive(),
     * TweakFunction::multiplicative() etc., or level passed to
     * TweakFunction::setter()
     */

    IObjectConstSP argument;                     // $optional

    // extracted from 'property', in validatePop2Object

    IScalarRiskPropertyConstSP _property;        // $transient
    double _coeff;                               // $transient

public:

    /**
     * IPerturbation implementation
     *
     * Just delegates to the ScalarFieldRiskProperty returned from
     * IFieldRiskPropertyDefinition::scalarProperty().
     */

    bool findAndShift(IObjectSP objectToShift,
                      OutputNameConstSP name) {
        IHypothesisConstSP hyp = _property->axisFor(name)->hypothesis(_coeff);
        // applyScenario should be a const method
        return IHypothesisSP::constCast(hyp)->applyScenario(objectToShift);
    }

    bool applyBeforeGetMarket() const {
        return property->applyBeforeGetMarket();
    }

    ~ScalarFlexiblePerturbation() {}
};

CClassConstSP const ScalarFlexiblePerturbation::TYPE = CClass::registerClassLoadMethod(
    "ScalarFlexiblePerturbation", typeid(ScalarFlexiblePerturbation), load);

bool ScalarFlexiblePerturbationLinkIn() {
    return ScalarFlexiblePerturbation::TYPE != NULL;
}

DRLIB_END_NAMESPACE
