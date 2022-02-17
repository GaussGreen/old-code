/**
 * @file PropertyRiskAxis.hpp
 */

#ifndef DRLIB_PropertyRiskAxis_H
#define DRLIB_PropertyRiskAxis_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/RiskAxis.hpp"
#include "edginc/IQualifiedRiskAxis.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(OutputName)
FORWARD_DECLARE(IAbstractRiskProperty)

template <class TAG>
class RiskProperty;

/**
 * A property of a particular market object, with respect to which risk can be
 * estimated.
 *
 * By IRiskAxis we mean a 1D manifold along which we can change the state of
 * the world.  A PropertyRiskAxis<TAG> is a risk axis which parameterises
 * changes to a particular property of a particular market name.  For
 * instance PropertyRiskAxis<Spot>("IBM") denotes IBM spot.
 *
 * The key method is hypothesis(), which gives you a
 * PropertyTweakHypothesis<TAG> representing a change to this property
 * of your chosen magnitude.
 *
 * The only subtlety is Qualifier---see RiskProperty::axisFor() for an
 * explanation.
 *
 * For an overview of how these classes fit into the "declarative"
 * sensitivities framework, see the IRiskQuantityFactory class documentation
 * (particularly <I>The new setup: property tags</I>).
 */

template <class TAG>
class PropertyRiskAxis:
    public CObject,
    public virtual IQualifiedRiskAxis<typename TAG::Qualifier> {

public:

    static CClassConstSP const TYPE;

    /**
     * Tag class for the property.
     *
     * Could be used to specify e.g. proportional vs absolute tweaks, or
     * whatever. For most properties it's an empty class --- see e.g. Spot.
     */

    typedef TAG Tag;
    typedef smartConstPtr<Tag> TagConstSP;

    /**
     * Info needed to identify a 1D risk axis along TAG on a given market
     * data name.
     *
     * For term-structured properties (like VolPointwise), this is
     * ExpiryWindow.  For scalar properties (like Spot) it's Void.  See
     * RiskProperty::axisFor();
     */

    typedef typename TAG::Qualifier Qualifier;
    typedef smartConstPtr<Qualifier> QualifierConstSP;

private:

    PropertyRiskAxis():
        CObject(TYPE)
    {}

    static IObject* defaultOne() {
        return new PropertyRiskAxis();
    }

    typedef IQualifiedRiskAxis<typename TAG::Qualifier> Super;

    static void load(CClassSP& clazz) {
        REGISTER(PropertyRiskAxis, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(Super);
        EMPTY_SHELL_METHOD(defaultOne);
        FIELD(_marketDataName, "marketDataName");
        FIELD(_qualifier, "qualifier");
        FIELD(_tag, "tag");
        RiskAxis::registerSimple(
            clazz->getName(), "_marketDataName", "_qualifier");
    }

    RiskAxisConstSP frozen() const {
        return RiskAxis::frozenSimple(this);
    }

    OutputNameConstSP _marketDataName;
    QualifierConstSP _qualifier;
    TagConstSP _tag;

public:

    /**
     * The market name to which we parameterise perturbations.
     */

    OutputNameConstSP marketDataName() const {
        return _marketDataName;
    }

    QualifierConstSP qualifier() const {
        return _qualifier;
    }

    TagConstSP tag() const {
        return _tag;
    }

    IAbstractRiskPropertyConstSP abstractProperty() const {
        return IAbstractRiskPropertyConstSP(
            new RiskProperty<TAG>(_tag));
    }

    /**
     * Constructor.
     *
     * Called by RiskProperty::axisFor(), which see.
     * Now allows construction without a qualifier or parameter
     */

    PropertyRiskAxis(OutputNameConstSP marketDataName,
                     QualifierConstSP qualifier = QualifierConstSP(),
                     TagConstSP tag = TagConstSP()):
        CObject(TYPE),
        _marketDataName(marketDataName),
        _qualifier(qualifier),
        _tag(tag)
    {}

    /**
     * A hypothesis that our property is changed by a given amount for our
     * market name.
     *
     * Returns a PropertyTweakHypothesis<TAG>.
     */

    IHypothesisConstSP hypothesis(double coeff) const {
        return coeff == 0 ?
            IHypothesis::noop() :
            IHypothesisSP(new PropertyTweakHypothesis<TAG>(
                coeff, _marketDataName, _qualifier, _tag));
    }

    string toString() const {
        return getClass()->getName() + "(" +
               (!_marketDataName ? "" : _marketDataName->toString()) +
               (!_qualifier ? "" : " " + _qualifier->toString()) + ")";
    }
private:
    DLL_FIX_FOR_TEMPLATE_TYPE; // work around for VC71 and dlls
};

#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
//// work around problem of trying to export static fields of templates
//// across dlls. For each array template used there needs to be an
//// explicit specialisation of TYPE defined which specifies the name and
//// load method
template <class TAG> CClassConstSP const PropertyRiskAxis<TAG>::TYPE =
CClass::templateRegisterClass(typeid(PropertyRiskAxis<TAG>));
#endif

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
