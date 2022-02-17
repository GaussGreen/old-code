/**
 * @file PropertyTweakHypothesis.hpp
 */

#ifndef QLIB_PropertyTweakHypothesis_H
#define QLIB_PropertyTweakHypothesis_H

#include "edginc/IRiskAxis.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/PropertyTweak.hpp"
#include "edginc/AbstractPropertyTweakHypothesis.hpp"

DRLIB_BEGIN_NAMESPACE

template <class TAG>
class PropertyRiskAxis;

/**
 * A function which maps the "world" to an alternate world in which some "risk
 * property" of a named market object has been modified: for instance, "IBM 3M
 * vol is 3bp bigger".
 *
 * This is the primary implementation of IHypothesis for the "declarative"
 * sensitivities framework (see IRiskQuantityFactory for an overview).  An
 * IHypothesis maps a world (= in practice a TweakGroup instrument/model
 * assembly) to an "alternate world"; a PropertyTweakHypothesis is a hypothesis
 * under which a "risk property"---like spot, vol or some other pricer
 * input---takes a different value for a particular market name.
 *
 * The property to be tweaked is specified by the TAG template parameter: for
 * instance Spot or VolPointwise.  These mostly don't include any actual
 * methods or fields of their own; they just parameterise a family of
 * templates.  See <I>The new setup: property tags</I> in the
 * IRiskQuantityFactory class documentation for an overview.
 *
 * The mechanism for a tweak is:
 *
 *    -  We find all objects of the right type and name, i.e. those
 *       implementing ITweakableWithRespectTo<TAG> and returning
 *       our marketDataName from ITweakableWithRespectTo<TAG>::sensName().
 *
 *    -  We call ITweakableWithRespectTo<TAG>::sensShift() on each
 *       of them, allowing them to change their internal state in a suitable
 *       way.
 *
 *    -  [When we've finished with the alternate world, we check
 *       which objects implement IRestorableWithRespectTo<TAG>
 *       and call IRestorableWithRespectTo<TAG>::sensRestore().]
 *
 * PropertyTweakHypothesis<TAG> is in practice created by
 * PropertyRiskAxis<TAG>, which in turn is created by
 * RiskProperty<TAG>.
 *
 * To minimise the volume of templated code, this class and the others in the
 * family have non-templated abstract base classes: see
 * AbstractPropertyTweakHypothesis, and <I>Abstract interfaces</I> in the
 * IRiskQuantityFactory class documentation.  In particular the key method
 * offered by this class is AbstractPropertyTweakHypothesis::appliedTo().
 */

template <class TAG>
class PropertyTweakHypothesis: public AbstractPropertyTweakHypothesis,
                               public PropertyTweak<TAG> {

public:

    static CClassConstSP const TYPE;

    /**
     * Tag class for the property to be tweaked
     *
     * Could be used to specify e.g. proportional vs absolute tweaks, or
     * whatever.  For most properties it's an empty class --- see e.g. Spot.
     */

    typedef TAG Tag;
    DECLARE(Tag)

    /**
     * Info needed in addition to TAG type and market data name to fully
     * specify the hypothesis.
     *
     * For term-structured properties (like VolPointwise), this is
     * ExpiryWindow.  For scalar properties (like Spot) it's Void.
     */

    typedef typename TAG::Qualifier Qualifier;
    DECLARE(Qualifier)

private:

    PropertyTweakHypothesis(const PropertyTweakHypothesis &);
    PropertyTweakHypothesis& operator =(const PropertyTweakHypothesis &);

    static IObject* defaultTweakHypothesis() {
        return new PropertyTweakHypothesis(0., OutputNameConstSP());
    }

    static void load(CClassSP& clazz) {
        REGISTER(PropertyTweakHypothesis, clazz);
        clazz->setPublic();
        SUPERCLASS(AbstractPropertyTweakHypothesis);
        EMPTY_SHELL_METHOD(defaultTweakHypothesis);
        FIELD(tag, "tag");
        FIELD_MAKE_OPTIONAL(tag);
        FIELD(qualifier, "qualifier");
        FIELD(marketDataName, "marketDataName");
        FIELD(coefficient, "coefficient");
    };

protected:


    /**
     * The interface which designates market objects affected by this hypothesis.
     *
     * It's ITweakableWithRespectTo<TAG>.
     */

    CClassConstSP shiftableInterface() const {
        return ITweakableWithRespectTo<TAG>::TYPE;
    }

    /**
     * Make a change to TAG of a particular object.
     *
     * Calls ITweakableWithRespectTo<TAG>::sensShift() on @a object.
     */

    TweakOutcome sensShift(IObjectSP object) const {
        try {
            return dynamic_cast<ITweakableWithRespectTo<TAG>&>(*object).
                sensShift(*this);
        }
        catch (exception& e) {
            throw ModelException(e, "PropertyTweakHypothesis::sensShift");
        }
    }

    /**
     * The interface which designates market objects which can be
     * mutated in place to apply this hypothesis.
     *
     * It's IRestorableWithRespectTo<TAG>.
     */

    CClassConstSP restorableInterface() const {
        return IRestorableWithRespectTo<TAG>::TYPE;
    }

    /**
     * Undo the effect of the last sensShift().
     *
     * Calls IRestorableWithRespectTo<TAG>::sensRestore() on @a object.
     */

    void sensRestore(IObjectSP object) const {
        try {
            dynamic_cast<IRestorableWithRespectTo<TAG>&>(*object).
                sensRestore(*this);
        }
        catch (exception& e) {
            throw ModelException(e, "PropertyTweakHypothesis::sensRestore");
        }
    }

    /**
     * "Market data name" of an object.
     *
     * Calls ITweakableWithRespectTo<TAG>::sensName() on @a object.  Some
     * objects take control of applying property tweaks to their
     * sub-components, and will return the latter's name in place of their own.
     * An object falls into the "domain" of a PropertyTweakHypothesis if
     * its sensName matches getMarketDataName().
     */

    string sensName(IObjectConstSP object) const {
        try {
            return dynamic_cast<const ITweakableWithRespectTo<TAG>&>(*object).
                sensName(this->tag.get());
        }
        catch (exception& e) {
            throw ModelException(e, "PropertyTweakHypothesis::sensName");
        }
    }

    /**
     * "Market data name" of the object whose TAG changes under this
     * hypothesis.
     *
     * Originates as @a subjectName argument to RiskProperty::axisFor().
     *
     * See sensName().
     */

    OutputNameConstSP getMarketDataName() const {
        return this->marketDataName;
    }

    /**
     * Info needed in addition to getMarketDataName() to fully specify the
     * hypothesis.
     *
     * E.g. for pointwise tweaks (term-specific hypotheses) this is an
     * ExpiryWindow.  It's really of type TAG::Qualifier.  You shouldn't
     * need to know this or what type this is; it's only used in toString().
     *
     * Originates as @a qualifier argument to RiskProperty::axisFor().
     */

    IObjectConstSP getQualifier() const {
        return this->qualifier;
    }

public:

    /**
     * Constructor.
     *
     * (These objects are mostly constructed via the factory method
     * PropertyRiskAxis<TAG>::hypothesis().)
     *
     * @param coefficient    How much the market object's TAG
     *                       gets shifted under the hypothesis.
     *                       Units are not necessarily absolute,
     *                       the exact semantics being specified by
     *                       the individual object's
     *                       ITweakableWithRespectTo<TAG>::sensShift()
     *                       implementation.
     *
     * @param marketDataName Names the market object affected by the
     *                       hypothesis.  May be null to mean "all".
     *
     * @param qualifier      Info needed in addition to TAG type and
     *                       @a marketDataName to fully specify the hypothesis.
     *                       Its type, Qualifier, is typedef'd to the
     *                       TAG::Qualifier specified in the tag class for
     *                       the property.  For pointwise tweaks
     *                       (maturity-specific hypotheses) this is an
     *                       ExpiryWindow --- see e.g. VolPointwise. For
     *                       non-term-structured tweaks it's Void --- see
     *                       e.g. Spot --- and may be omitted.
     *
     * @param tag            Tag object for the property being tweaked:
     *                       could be used to specify e.g. proportional vs
     *                       absolute tweaks, or whatever.  For most properties
     *                       it's an empty class -- see e.g. Spot --- and may
     *                       be omitted.
     */

    PropertyTweakHypothesis(
            double coefficient,
            OutputNameConstSP marketDataName = OutputNameConstSP(),
            QualifierConstSP qualifier = QualifierConstSP(),
            TagConstSP tag = TagConstSP()):
        AbstractPropertyTweakHypothesis(TYPE),
        PropertyTweak<TAG>(coefficient, marketDataName, qualifier, tag)
    {}

    /**
     * Constructor, from a preexisting PropertyTweak
     */

    PropertyTweakHypothesis(const PropertyTweak<TAG>& tweak):
        AbstractPropertyTweakHypothesis(TYPE),
        PropertyTweak<TAG>(tweak.coefficient, tweak.marketDataName,
                           tweak.qualifier, tweak.tag)
    {}

    /**
     * The "risk axis" along which this hypothesis shifts the world.
     *
     * It's a PropertyRiskAxis<TAG>.
     */

    IRiskAxisConstSP axis() const {
        return IRiskAxisConstSP(
            new PropertyRiskAxis<TAG>(this->marketDataName, 
                                      this->qualifier,
                                      this->tag));
    }

    /**
     * The coefficient governing the amount by which this hypothesis
     * shifts the world along axis().
     */

    double axisCoefficient() const {
        return this->coefficient;
    }
    /**
     * Support for Control::getCurrentSensitivity() compatibility
     *
     * These methods are provided so that we can implement a compatibility
     * mechanism for model/instrument logic which relies on interrogating
     * Control::getCurrentSensitivity() --- see <I>Current status</I> in the
     * IRiskQuantityFactory class documentation.
     */

    //@{

    void setMarketDataName(OutputNameConstSP name) {
        this->marketDataName = name;
    }

    void setAxisCoefficient(double c) {
        this->coefficient = c;
    }

    //@}
private:
    DLL_FIX_FOR_TEMPLATE_TYPE; // work around for VC71 and dlls
};

#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
//// work around problem of trying to export static fields of templates
//// across dlls. For each array template used there needs to be an
//// explicit specialisation of TYPE defined which specifies the name and
//// load method
template <class TAG> CClassConstSP const 
PropertyTweakHypothesis<TAG>::TYPE =
CClass::templateRegisterClass(typeid(PropertyTweakHypothesis<TAG>));
#endif

DRLIB_END_NAMESPACE

#include "edginc/PropertyRiskAxis.hpp"

#endif
