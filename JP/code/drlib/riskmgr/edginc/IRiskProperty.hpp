/**
 * @file IRiskProperty.hpp
 */

#ifndef DRLIB_IRiskProperty_H
#define DRLIB_IRiskProperty_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Void.hpp"
#include "edginc/IAbstractRiskProperty.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IHypothesis)
FORWARD_DECLARE(IRiskAxis)
FORWARD_DECLARE(Expiry)
FORWARD_DECLARE(MultiTweakGroup)

class ExpiryWindow;
class ExpiryPair;
class ExpiryAndStrike;
class BoxedInt;
class Sensitivity; // for adaptiveCoefficients() below :(

/**
 * A property which a market object can have, and with respect to which risk
 * can be estimated
 *
 * Implementations of this interface represent "properties" which can be
 * imputed to market data, and feed into the pricing process.  For instance,
 * RiskProperty<VolPointwise> denotes term-structured volatility.
 *
 * Most properties are implemented using the RiskProperty<PROPERTY> template,
 * which see for documentation of the methods.  fieldRiskProperty
 * is a "dynamically typed" counterpart.
 *
 * The methods which don't require knowledge of the QUALIFIER type are lifted
 * into a non-templated interface IAbstractRiskProperty so that we don't have
 * to put templated methods everywhere.
 *
 * For an overview of how these classes fit into the "declarative"
 * sensitivities framework, see the IRiskQuantityFactory class documentation.
 */

template <class QUALIFIER>
class RISKMGR_DLL IRiskProperty: public virtual IAbstractRiskProperty {

    static void load(CClassSP& clazz);

public:

    /**
     * See RiskProperty<PROPERTY>::Qualifier.
     */

    typedef QUALIFIER Qualifier;
    DECLARE(Qualifier)

    static CClassConstSP const TYPE;

    IRiskProperty();
    ~IRiskProperty();

    /**
     * See RiskProperty<PROPERTY>::axisFor().
     */

    virtual IRiskAxisConstSP axisFor(
        OutputNameConstSP subjectName,
        QualifierConstSP qualifier = QualifierConstSP()) const = 0;

    /**
     * See RiskProperty<PROPERTY>::subjectQualifiers().
     */

    virtual QualifierArrayConstSP subjectQualifiers(
        IObjectConstSP world,
        OutputNameConstSP name) const = 0;

    virtual BoolArrayConstSP mayHaveEffect(
        MultiTweakGroupConstSP world,
        OutputNameConstSP name,
        QualifierArrayConstSP qualifiers,
        const Sensitivity* sensitivity) const;

    virtual DoubleArrayConstSP adaptiveCoefficients(
        MultiTweakGroupConstSP world,
        OutputNameConstSP name,
        QualifierArrayConstSP qualifiers,
        double targetCoefficient) const;

     /**
     * Tweak this property on all objects in the world which have it.
     *
     * Used by a few basket-like market data objects to tweak their components
     * w.r.t. Spot.
     */

    void tweakAllSubjects(IObjectSP world, QualifierConstSP qualifier,
                          double coefficient) const;

    /**
     * An IRiskProperty considered under a given scenario
     *
     * Takes an IRiskProperty @a property, and returns one representing that
     * same property "after" a given scenario has been applied to the world.
     * The IHypothesis's ultimately arising from property via axisFor() and
     * IRiskAxis::hypothesis() are prefixed by @a hypothesis (implementation
     * uses IRiskAxis::conditioned()).
     *
     * For instance, prefixing a RiskProperty<Spot> with a time-shift
     * hypothesis would give you a "delta next day" property.  Or,
     * the greek LiquiditySpreadRhoPointwise applies this to
     * RiskProperty<LiquiditySpreadPointwise> to obtain the property
     * "liquidity spread with par spread pricing turned off".
     */

    static smartPtr<IRiskProperty> conditioned(
        IHypothesisConstSP hypothesis, smartConstPtr<IRiskProperty> property);
};

#ifndef QLIB_IRISKPROPERTY_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL IRiskProperty<Void>);
EXTERN_TEMPLATE(class RISKMGR_DLL IRiskProperty<ExpiryWindow>);
EXTERN_TEMPLATE(class RISKMGR_DLL IRiskProperty<ExpiryPair>);
EXTERN_TEMPLATE(class RISKMGR_DLL IRiskProperty<ExpiryAndStrike>);
EXTERN_TEMPLATE(class RISKMGR_DLL IRiskProperty<BoxedInt>);
#endif

/**
 * A non-term-structured property which a market object can have.
 */

typedef IRiskProperty<Void> IScalarRiskProperty;
DECLARE(IScalarRiskProperty)

/**
 * A term-structured property which a market object can have.
 */

typedef IRiskProperty<ExpiryWindow> IExpiryRiskProperty;
DECLARE(IExpiryRiskProperty)

/**
 * A 2d-term-structured property which a market object (i.e CDS vol objects) can have.
 */

typedef IRiskProperty<ExpiryPair> IExpiryPairRiskProperty;
DECLARE(IExpiryPairRiskProperty)

/**
 * A 2d-term-structured property which a market object (i.e vol surfs) can have.
 */

typedef IRiskProperty<ExpiryAndStrike> IExpiryAndStrikeRiskProperty;
DECLARE(IExpiryAndStrikeRiskProperty)

/**
 * An array-valued property which a market object can have.
 */

typedef IRiskProperty<BoxedInt> IIntRiskProperty;
DECLARE(IIntRiskProperty)

DRLIB_END_NAMESPACE

#endif
