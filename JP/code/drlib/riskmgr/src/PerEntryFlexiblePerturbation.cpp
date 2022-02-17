/**
 * @file PerEntryFlexiblePerturbation.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/Perturbation.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/IFieldRiskPropertyDefinition.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * An IPerturbation which specifies a change to the value of a named
 * array-field on some class, at indices corresponding to a given set of
 * qualifiers
 *
 * This is client-friendly vehicle for specifying a ScenarioShift under which a
 * field designated through the QLib reflection system, of type DoubleArray or
 * CDoubleMatrix and indexed by some "qualifier" such as Expiry, has its value
 * adjusted at indices corresponding to a given set of qualifiers.  It doesn't
 * have to be directly attached to the "anchor" class, and there are several
 * operations available for shifting or setting it.  See
 * IFieldRiskPropertyDefinition.
 *
 * If you want to construct scenarios on the fly inside QLib, you'll find it
 * much more convenient to use ParallelFieldRiskProperty directly.
 * 
 * See FlexibleSensitivity for an overview of "flexible scenarios and
 * greeks", and IRiskQuantityFactory for background info on the
 * "declarative sensitivities" framework in which they're implemented.
 */

class PerEntryFlexiblePerturbation: public CObject,
                                    public virtual IPerturbation {

    PerEntryFlexiblePerturbation(): CObject(TYPE) {}

    static IObject* emptyShell() {
        return new PerEntryFlexiblePerturbation();
    }

    static void load(CClassSP& clazz) {
        REGISTER(PerEntryFlexiblePerturbation, clazz);
        clazz->setPublic();
        SUPERCLASS(CObject);
        IMPLEMENTS(IPerturbation);
        EMPTY_SHELL_METHOD(emptyShell);

        FIELD(property, "property");
        FIELD(arguments, "arguments");
        FIELD(qualifiers, "qualifiers");

        FIELD(_property, "_property");
        FIELD_MAKE_TRANSIENT(_property);
    }

    void validatePop2Object() {
        try {
            _property = property->parallelProperty(qualifiers, arguments);

            if (arguments->getLength() != qualifiers->getLength()) {
                throw ModelException(
                    "'arguments' are a different length than 'qualifiers'");
            }
        }
        catch (exception& e) {
            throw ModelException(
                e, "PerEntryFlexiblePerturbation::validatePop2Object()");
        }
    }

public:

    static CClassConstSP const TYPE;

private:

    /**
     * Details of the class, its field and how we tweak it
     *
     * IFieldRiskPropertyDefinition is just a factory for FieldRiskProperty's,
     * specified in a client-friendly way
     */

    IFieldRiskPropertyDefinitionConstSP property; // $required

    /**
     * Expiries to be tweaked on the field
     *
     * The association of qualifiers to indices is via the field named
     * in FieldRiskPropertyDefinition::qualifierField.
     */

    IArrayConstSP qualifiers;                     // $required

    /**
     * Shift size passed to TweakFunction::additive(),
     * TweakFunction::multiplicative() etc., or level passed to
     * TweakFunction::setter(), for each expiry
     */

    IArrayConstSP arguments;                     // $required

    // extracted from 'property', in validatePop2Object

    IScalarRiskPropertyConstSP _property;        // $transient

public:

    /**
     * IPerturbation implementation
     *
     * Just delegates to the ParallelFieldRiskProperty returned from
     * IFieldRiskPropertyDefinition::paralleProperty()
     */

    bool findAndShift(IObjectSP objectToShift,
                      OutputNameConstSP name) {
        IHypothesisConstSP hyp = _property->axisFor(name)->hypothesis(1.);
        return IHypothesisSP::constCast(hyp)->applyScenario(objectToShift);
    }

    bool applyBeforeGetMarket() const {
        return property->applyBeforeGetMarket();
    }

    ~PerEntryFlexiblePerturbation() {}
};

CClassConstSP const PerEntryFlexiblePerturbation::TYPE = CClass::registerClassLoadMethod(
    "PerEntryFlexiblePerturbation", typeid(PerEntryFlexiblePerturbation), load);

bool PerEntryFlexiblePerturbationLinkIn() {
    return PerEntryFlexiblePerturbation::TYPE != NULL;
}

DRLIB_END_NAMESPACE
