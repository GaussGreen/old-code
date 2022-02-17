/**
 * @file PerNameRiskPropertySensitivity.hpp
 */

#ifndef QLIB_PerNameRiskPropertySensitivity_H
#define QLIB_PerNameRiskPropertySensitivity_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IRiskProperty.hpp"
#include "edginc/RiskPropertySensitivity.hpp"
#include "edginc/IPerNameSensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IScalarDerivative)
FORWARD_DECLARE(IResultsFunction)

/**
 * A greek calculated by tweaking a "risk property" (like spot or vol) on each
 * of the "names" in the market in turn
 *
 * This is part of the "declarative" sensitivities framework: see
 * IRiskQuantityFactory for an overview.  It's the base class for greeks
 * estimated by tweaking a specified RiskProperty (like Spot or VolPointwise)
 * on all or some of the "names" which it finds in the world.  The key
 * method is nameRiskQuantities().
 *
 * Why's this class a template?  Because the abstract IRiskProperty interface
 * under which we store our 'property' field is parameterised by the type of
 * "qualifier" it requires (ExpiryWindow for term-structured properties like
 * VolPointwise, Void for scalars like Spot).
 *
 *
 * <H3>Implementations</H3>
 *
 * Implementations at the time of writing include
 *
 *    -  ScalarRiskPropertySensitivity, for "scalar" greeks
 *
 *       -  ParSpreadRhoParallel.cpp
 *       -  CRMeanReversionParallel.cpp
 *       -  ...
 *
 *    -  PerNameRiskPropertySensitivity<ExpiryWindow>, for term-structured
 *       greeks
 *
 *       -  VegaPointwise.cpp
 *       -  ParSpreadRhoPointwise.cpp
 *       -  ...
 *
 *    -  PerNameRiskPropertySensitivity<BoxedInt>, for array-valued greeks
 *       indexed by plain int rather than Expiry
 *
 *    -  PerNameRiskPropertySensitivity<ExpiryPair>, for greeks indexed by
 *       two expiries
 *
 *       - CRVegaPointwise.cpp
 *
 * PerNameRiskPropertySensitivity corresponds roughly to SensControlPerName in
 * the existing Sensitivity framework, but also subsumes VectorShift,
 * i.e. greeks indexed by ExpiryWindow, as well as array-valued greeks indexed
 * by plain integers, and greeks indexed by ExpiryPair.
 * ScalarRiskPropertySensitivity corresponds to ScalarShift (it's basically an
 * instantiation of this template class but adds some stuff for compatibility).
 *
 * As a special case, for the moment, Delta and BasketDelta are <I>not</I>
 * descended from PerNameRiskPropertySensitivity, but rather from ScalarShift via
 * ScalarShiftTwoSided.  This is because there are some fiddly mechanisms (like
 * MonteCarlo "quick greeks") which are written in terms of ScalarShift and I
 * don't want to have to fix them up before committing this stuff.
 */

template <class QUALIFIER>
class PerNameRiskPropertySensitivity: public RiskPropertySensitivity<QUALIFIER>,
                                      public virtual IPerNameSensitivity {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

    // Duplicate from parent class to avoid "implicit typename" warnings
    
    typedef QUALIFIER Qualifier;
    DECLARE(Qualifier)

protected:

    /**
     * Constructor
     *
     * This is really the first half of the constructor, the second half being
     * RiskPropertySensitivity::deriv() --- that's called after the object's
     * reflection FIELDs have been filled in so that it can refer to their
     * values.  The fields of the Deriv object you return from deriv() are
     * listed here as if they were constructor arguments, so that you can see
     * everything in one place.
     *
     * @param type             Reflection type of final subclass
     *
     * @param shiftSize        How big to make the @a property tweaks used to
     *                         estimate @a derivative
     *
     * @param outputName       Packet name under which quantities computed for
     *                         this greek will be stored in the Results
     *                         (e.g. "PAR_SPREAD_RHO_PARALLEL" for
     *                         ParSpreadRhoParallel.cpp)
     *
     * @param outputName2      Packet name under which quantities computed for
     *                         the optional second derivative will be stored
     *                         in the results --- see the alternative
     *                         constructor on Deriv
     *
     * @param derivand         The quantity whose derivative we want to
     *                         estimate: generally fair value (i.e.
     *                         IResultsFunction::price()) but it could be
     *                         NAKED_BOND_PRICE or anything
     * 
     * @param property         The IRiskProperty with respect to which the
     *                         @a derivative is to be estimated (e.g.
     *                         RiskProperty<ParSpreadParallel> for
     *                         ParSpreadRhoParallel.cpp)
     *
     * @param derivative       The type of derivative to estimate, e.g.
     *                         IScalarDerivative::oneSided(),
     *                         IScalarDerivative::twoSided(),
     *                         IScalarDerivative::second()
     *
     * @param sensitivityUnit  Scaling factor applied to computed derivatives
     *                         before they're placed in the Results.  Normally
     *                         it should just be 1.0, but people want to see
     *                         some sensitivities reported in terms of basis
     *                         points (0.0001).
     *
     * @param derivative2      An optional other derivative to compute
     *                         (used by two-sided greeks which provide also
     *                         gamma)
     *
     * @param sensitivityUnit2 Goes with derivative2
     */

    PerNameRiskPropertySensitivity(CClassConstSP type,
                                   double shiftSize,
                                   const string& outputName,
                                   const string& outputName2 = "");

    /**
     * The risk quantities which this PerNameRiskPropertySensitivity says
     * should be computed for a given base-state world (one per name).
     *
     * Called by RiskQuantityFactorySensitivity::riskQuantities() to implement
     * IRiskQuantityFactory::riskQuantities().
     */

    virtual NamedRiskQuantityArraySP nameRiskQuantities(
        MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping) const;

public:

    ~PerNameRiskPropertySensitivity();

    /**
     * IPerNameSensitivity implementation
     *
     * These methods allow PerNameRiskPropertySensitivity to be a drop-in
     * replacement for SensControlPerName in model/instrument logic which (for
     * the moment) relies on Control::getCurrentSensitivity().  See
     * IPerNameSensitivity for more info.
     */

    //@{

    virtual OutputNameConstSP getMarketDataName() const;

    virtual OutputNameArrayConstSP allNames(const IObject* object) const;

    virtual CClassConstSP shiftInterface() const;

    //@}

    /**
     * A dynamically constructed PerNameRiskPropertySensitivity
     *
     * Analogously with those in the SensControl framework, sensitivities in
     * the new "declarative" framework are defined as subclasses of (mostly)
     * this class or ScalarRiskPropertySensitivity, using the "protected"
     * constructor above.  However, you can actually make a fully working
     * PerNameRiskPropertySensitivity without further inheritance (the
     * constructor arguments are enough).
     */

    static smartPtr<PerNameRiskPropertySensitivity> SP(
        const string& outputName,
        IResultsFunctionConstSP derivand,
        smartConstPtr<IRiskProperty<QUALIFIER> > property,
        IScalarDerivativeConstSP derivative,
        double sensitivityUnit,
        double coefficient);

    /**
     * A dynamically constructed PerNameRiskPropertySensitivity which computes
     * two derivatives for each name/qualifier
     */

    static smartPtr<PerNameRiskPropertySensitivity> SP(
        const string& outputName1,
        const string& outputName2,
        IResultsFunctionConstSP derivand,
        smartConstPtr<IRiskProperty<QUALIFIER> > property,
        IScalarDerivativeConstSP derivative1,
        IScalarDerivativeConstSP derivative2,
        double sensitivityUnit1,
        double sensitivityUnit2,
        double coefficient);

    /**
     * Equivalent of VectorShift::setExpiryToTweak().  This is only used in
     * Duration and is considered temporary.
     */

    smartPtr<PerNameRiskPropertySensitivity> withQualifier(
            QualifierConstSP qualifier) const;
};

/**
 * Base class for term-structured greeks
 */

#ifndef QLIB_PERNAMERISKPROPERTYSENSITIVITY_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL PerNameRiskPropertySensitivity<ExpiryWindow>);
EXTERN_TEMPLATE(class RISKMGR_DLL PerNameRiskPropertySensitivity<BoxedInt>);
EXTERN_TEMPLATE(class RISKMGR_DLL PerNameRiskPropertySensitivity<ExpiryPair>);
EXTERN_TEMPLATE(class RISKMGR_DLL PerNameRiskPropertySensitivity<ExpiryAndStrike>);
EXTERN_TEMPLATE(class RISKMGR_DLL PerNameRiskPropertySensitivity<Void>);
#endif

typedef PerNameRiskPropertySensitivity<ExpiryWindow> VectorRiskPropertySensitivity;
DECLARE(VectorRiskPropertySensitivity)

/**
 * Base class for array-valued greeks
 */

typedef PerNameRiskPropertySensitivity<BoxedInt> ElementwiseRiskPropertySensitivity;
DECLARE(ElementwiseRiskPropertySensitivity)

/**
 * Base class for greeks indexed by two expiries
 */

typedef PerNameRiskPropertySensitivity<ExpiryPair> ExpiryPairRiskPropertySensitivity;
DECLARE(ExpiryPairRiskPropertySensitivity)

/**
 * Base class for greeks indexed by ExpiryAndStrike
 */

typedef PerNameRiskPropertySensitivity<ExpiryAndStrike> ExpiryAndStrikeRiskPropertySensitivity;
DECLARE(ExpiryAndStrikeRiskPropertySensitivity)

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
