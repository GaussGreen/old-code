/**
 * @file RiskProperty.hpp
 */

#ifndef DRLIB_RiskProperty_H
#define DRLIB_RiskProperty_H

#include "edginc/TweakNameListID.hpp"
#include "edginc/IRiskProperty.hpp"
#include "edginc/PropertyRiskAxis.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/TweakNameResolver.hpp"
#include "edginc/RiskMappingMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * A property which a market object can have, and with respect to which risk
 * can be estimated.
 *
 * Instances of this template represent "properties" which can be imputed to
 * market data, and feed into the pricing process.  For instance,
 * RiskProperty<VolPointwise> denotes term-structured volatility.
 *
 * Associating it with a name ("IBM") and a qualifier ("3 month") gives you a
 * PropertyRiskAxis<VolPointwise>.  Our key method is the axisFor() factory
 * method, which returns you that as an abstract IRiskAxis.
 *
 * [Adding a coefficient for how far to shift the world along that axis ("1e-4")
 * gives you a PropertyTweakHypothesis<VolPointwise>---essentially a scenario
 * ("IBM 3 month vol is 1bp higher") against which you can measure price
 * sensitivity.  See PropertyRiskAxis::hypothesis().]
 *
 * Note that the TAG template parameter is just a "tag", for instance Spot or
 * VolPointwise.  These mostly don't include any actual methods or fields of
 * their own; they just parameterise a family of templates.
 *
 * RiskProperty's primary use is as the @a property argument to
 * NameRiskPropertySensitivity::deriv().
 *
 * To reduce the volume of templated code, the abstract interface IRiskProperty
 * is used wherever possible.
 * 
 * For an overview of how these classes fit into the "declarative"
 * sensitivities framework, see the IRiskQuantityFactory class documentation
 * (particularly <I>The new setup: property tags</I>).
 */

template <class TAG>
class RiskProperty:
        public CObject,
        public virtual IRiskProperty<typename TAG::Qualifier> {

public:

    static CClassConstSP const TYPE;

    typedef TAG Tag;
    typedef smartConstPtr<Tag> TagConstSP;

    /**
     * Info needed in addition to TAG type and market data name to make
     * a one-dimensional "risk axis" from the property.
     *
     * See axisFor().  For term-structured properties (like VolPointwise), this
     * is ExpiryWindow.  For scalar properties (like Spot) it's Void.
     */

    typedef typename TAG::Qualifier Qualifier;
    DECLARE(Qualifier)

private:

    typedef ITweakableWithRespectTo<TAG> Tweakable;

    static IObject* emptyShell() { return new RiskProperty(); }

    typedef IRiskProperty<Qualifier> Super;

    static void load(CClassSP& clazz) {
        REGISTER(RiskProperty, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(Super);
        EMPTY_SHELL_METHOD(emptyShell);
        FIELD(_tag, "tag");
    }

    /**
     * Used in subjectNames()
     */

    struct NameListID: public ITweakNameListID {

        TagConstSP tag;

        NameListID(TagConstSP tag): tag(tag) {}

        typedef ITweakableWithRespectTo<TAG> Tweakable;

        CClassConstSP shiftInterface() const {
            return Tweakable::TYPE;
        }

        void appendName(OutputNameArray& namesList, IObjectConstSP obj) {
            string it = dynamic_cast<const Tweakable&>(*obj).sensName(
                            tag.get());
            if (it != "")
                namesList.push_back(OutputNameSP(new OutputName(it)));
        }
    };

    /**
     * Used in subjectQualifiers()
     */

    struct TweakNameResolver: public ITweakNameResolver {
        TagConstSP tag;
        OutputNameConstSP name;

        TweakNameResolver(TagConstSP tag,
                          OutputNameConstSP name):
            tag(tag),
            name(name)
        {}

        OutputNameConstSP getMarketDataName() const { return name; }

        bool nameMatches(const OutputName& name, IObjectConstSP obj) /* const */ {
            return name.equals(dynamic_cast<const Tweakable&>(*obj).
                                   sensName(tag.get()));
        }
    };

    TagConstSP _tag;

public:

    /**
     * Constructor.
     *
     * @param tag      Tag object for the property: could be used to specify
     *                 e.g. proportional vs absolute tweaks, or whatever.  For
     *                 most properties it's an empty class -- see e.g. Spot ---
     *                 and may be omitted.
     */

    RiskProperty(TagConstSP tag = TagConstSP()):
        CObject(TYPE),
        _tag(tag)
    {}

    /**
     * Constructor returning a smartPtr.
     */

    static smartPtr<RiskProperty> SP(
            TagConstSP tag = TagConstSP()) {
        return smartPtr<RiskProperty>(new RiskProperty(tag));
    }

    /**
     * Whether the property is continuous, i.e. tweaks to it can be made
     * arbitrarily small.
     *
     * Used by RiskPropertySensitivity::discreteShift().
     */

    bool discrete() const {
        return TAG::discrete;
    }

    /**
     * The interface designating market objects which can have this property.
     *
     * Having the property means you can tweak it, so this is
     * ITweakableWithRespectTo<TAG>.  For instance, Equity has
     * RiskProperty<Spot> because it implements ITweakableWithRespectTo<Spot>.
     */

    CClassConstSP subjectInterface() const {
        return Tweakable::TYPE;
    }

    /**
     * All the market data names in the world which have this "property".
     */

    OutputNameArrayConstSP subjectNames(IObjectConstSP world) const {
        NameListID namesListID(_tag);
        return OutputName::trim(SensMgrConst(world).allNames(&namesListID));
    }

    /**
     * The "qualifiers" defining the instances of this property on a given
     * name.
     *
     * Some properties, like Spot, are scalar for each name in the market;
     * others exist in multiple instances for each name, and require a
     * qualifier for full specification.  The primary use is for term
     * structured ("vector") properties like VolPointwise, but others can be
     * imagined.  See axisFor().
     *
     * This method finds the object in @a world which implements
     * subjectInterface() and is called @a name, and calls
     * ITweakableWithRespectTo<TAG> on it.
     *
     * So, for term structured properties (e.g. VolPointwise) it returns an
     * ExpiryArrayConstSP giving the expiries at which vol is defined for
     * @a name.
     *
     * For scalar properties (e.g. Spot) it returns Void.
     *
     * @throw
     *   ModelException if there is not exactly one matching object in @a world
     */

    QualifierArrayConstSP subjectQualifiers(IObjectConstSP world,
                                            OutputNameConstSP name) const {
        TweakNameResolver nameres(_tag, name);
        return dynamic_cast<const Tweakable&>(
                   *SensMgrConst(world).theFirst(Tweakable::TYPE, &nameres)
                 ).sensQualifiers(_tag.get());
    }

    /**
     * A "risk axis" aligned along a property of a given market name.
     *
     * By IRiskAxis we mean a 1D manifold along which we can change the state of
     * the world.  Moving along the risk axis returned by this method means
     * tweaking our property on a particular name.
     * 
     * Some properties exist in multiple instances for each name, and need
     * qualification.  Canonically, for term structured properties like
     * VolPointwise, you need to say which expiry you're tweaking, and
     * also how it fits into the overall sequence, so Qualifier is
     * ExpiryWindow (typedef'd from VolPointwise::Qualifier), and
     * @a qualifier must be an ExpiryWindowConstSP.
     *
     * For scalar properties like Spot, Qualifier is Void, and you
     * can just set @a qualifier to be VoidConstSP().
     * name @a subjectName, qualified if necessary by @a qualifier (e.g. expiry
     * for term-structured properties), by <I>x</I>.  I.e. calling
     * ITweakableWithRespectTo<TAG>::sensShift() with <I>x</I> as the
     * PropertyTweak::coefficient and @a qualifier as the
     * PropertyTweak::qualifier.
     *
     * Specifically, then, a move of coefficient <I>x</I> along the returned
     * IRiskAxis involves calling
     * ITweakableWithRespectTo<TAG>::sensShift() on all objects of type
     * subjectInterface() called @a subjectName, with <I>x</I> as the
     * PropertyTweak::coefficient and @a qualifier as the
     * PropertyTweak::qualifier.  See PropertyTweakHypothesis for the
     * implementation.
     */

    IRiskAxisConstSP axisFor(OutputNameConstSP subjectName,
                             QualifierConstSP qualifier) const {
        return IRiskAxisConstSP(new PropertyRiskAxis<TAG>(
            subjectName, qualifier, _tag));
    }

    /**
     * Allow the property to generate a RiskMappingMatrix on the fly.
     *
     * Default implementation returns a null SP
     */

    RiskMappingMatrixConstSP riskMappingMatrix(IObjectConstSP world,
                                               OutputNameConstSP name) const {
        //default implementation returns void
        //=> no dynamic risk mapping matrix for this property
        return RiskMappingMatrixSP();
    }

    /**
     * The property's parameters
     */

    TagConstSP tag() const {
        return _tag;
    }

private:
    DLL_FIX_FOR_TEMPLATE_TYPE; // work around for VC71 and dlls

    /**
     * For error messages
     */

    string toString() const {
        return getClass()->getName();
    }
};

#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
//// work around problem of trying to export static fields of templates
//// across dlls. For each array template used there needs to be an
//// explicit specialisation of TYPE defined which specifies the name and
//// load method
template <class TAG> CClassConstSP const RiskProperty<TAG>::TYPE =
CClass::templateRegisterClass(typeid(RiskProperty<TAG>));
#endif

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
