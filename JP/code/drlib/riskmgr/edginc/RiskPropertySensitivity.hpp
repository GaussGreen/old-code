/**
 * @file RiskPropertySensitivity.hpp
 */

#ifndef QLIB_RiskPropertySensitivity_H
#define QLIB_RiskPropertySensitivity_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IRiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/RiskQuantityFactorySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IScalarDerivative);
FORWARD_DECLARE(IResultsFunction);

/**
 * A greek calculated by tweaking a "risk property" (like spot or vol) on
 * "names" in the market.
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
 * Subclasses at the time of writing include
 *
 *    -  PerNameRiskPropertySensitivity, for greeks reported with respect
 *       to tweaks to each name individually
 *
 *       -  ScalarRiskPropertySensitivity, for "scalar" greeks
 *
 *          -  ParSpreadRhoParallel.cpp
 *          -  CRMeanReversionParallel.cpp
 *          -  ...
 *
 *       -  PerNameRiskPropertySensitivity<ExpiryWindow>, for term-structured
 *          greeks
 *
 *          -  VegaPointwise.cpp
 *          -  ParSpreadRhoPointwise.cpp
 *          -  ...
 *
 *    -  AllNamesRiskPropertySensitivity, for greeks reported with
 *       respect to a tweak to all names at the same time
 *
 * PerNameRiskPropertySensitivity corresponds roughly to SensControlPerName in
 * the existing Sensitivity framework, while ScalarRiskPropertySensitivity and
 * PerNameRiskPropertySensitivity<ExpiryWindow> correspond to ScalarShift and
 * VectorShift.
 *
 * As a special case, for the moment, Delta and BasketDelta are <I>not</I>
 * descended from RiskPropertySensitivity, but rather from ScalarShift via
 * ScalarShiftTwoSided.  This is because there are some fiddly mechanisms (like
 * MonteCarlo "quick greeks") which are written in terms of ScalarShift and I
 * don't want to have to fix them up before committing this stuff.
 */

template <class QUALIFIER>
class RiskPropertySensitivity: public RiskQuantityFactorySensitivity {

    RiskPropertySensitivity(const RiskPropertySensitivity& rhs);
    RiskPropertySensitivity& operator=(const RiskPropertySensitivity& rhs);
    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

    typedef QUALIFIER Qualifier;
    DECLARE(Qualifier)

protected:

    void ensureDeriv() const;
    mutable IResultsFunctionConstSP _derivand;
    mutable smartConstPtr<IRiskProperty<QUALIFIER> > _property;
    mutable IScalarDerivativeArrayConstSP _derivatives;
    mutable DoubleArrayConstSP _sensitivityUnits;

    /* const - except for ImpliedScalarShift */ double shiftSize;
    QualifierArrayConstSP predefinedQualifiers; // $unregistered
    bool hasPredefinedQualifiers; // $unregistered

    /**
     * Constructor (part 1)
     *
     * deriv() is conceptually the "other half" of this constructor.  We have
     * to split it up because the Deriv you specify may depend on reflection
     * FIELDs which haven't been filled in at construct time.
     *
     * @param type             Reflection type of final subclass
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
     * @param shiftSize        How big to make the @a property tweaks used to
     *                         estimate @a derivative
     */

    RiskPropertySensitivity(CClassConstSP type,
                            double shiftSize,
                            const string& outputName,
                            const string& outputName2 = "");

public: // VC6 needs this

    /**
     * Specification for the derivative to be calculated by a
     * RiskPropertySensitivity
     *
     * This what subclasses must return from deriv(), which is conceptually the
     * "second half" of the constructor RiskPropertySensitivity().
     */

    class Deriv {
    public:
        // these are explicit to avoid the need for the compiler to generate
        // them everytime it sees this class
        Deriv(const Deriv& rhs);
        Deriv& operator=(const Deriv& rhs);

        IResultsFunctionConstSP derivand;
        smartConstPtr<IRiskProperty<QUALIFIER> > property;
        IScalarDerivativeArrayConstSP derivatives;
        DoubleArrayConstSP sensitivityUnits;

        /**
         * Constructor
         *
         * These parameters are really the "second half" of RiskPropertySensitivity().
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

        Deriv(IResultsFunctionConstSP derivand,
              smartConstPtr<IRiskProperty<QUALIFIER> > property,
              IScalarDerivativeConstSP derivative,
              double sensitivityUnit,
              IScalarDerivativeConstSP derivative2 = IScalarDerivativeConstSP(),
              double sensitivityUnit2 = 0.);

        ~Deriv();

    };

protected:

    /**
     * Constructor (part 2): specification for the derivative to be calculated
     * by this RiskPropertySensitivity
     *
     * The Deriv::Deriv() constructor arguments are essentially the second half
     * of the RiskPropertySensitivity() constructor.  We have to split it up
     * because the Deriv you specify may depend on reflection FIELDs which
     * haven't been filled in at construct time.
     */

    virtual Deriv deriv() const = 0;

public:

    IResultsFunctionConstSP derivand() const;
    smartConstPtr<IRiskProperty<QUALIFIER> > property() const;
    IScalarDerivativeArrayConstSP derivatives() const;
    DoubleArrayConstSP sensitivityUnits() const;

    ~RiskPropertySensitivity();

    /**
     * Sensitivity implementation
     */

    //@{

    /**
     * Whether the property against which the sensitivity is measured is
     * continuous, i.e. tweaks to it can be made arbitrarily small.
     *
     * Returns RiskProperty::discrete() from the @a property argument
     * to RiskPropertySensitivity() constructor.
     */

    bool discreteShift() const;

    //@}
};

#ifndef QLIB_RISKPROPERTYSENSITIVITY_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL RiskPropertySensitivity<Void>);
EXTERN_TEMPLATE(class RISKMGR_DLL RiskPropertySensitivity<ExpiryWindow>);
EXTERN_TEMPLATE(class RISKMGR_DLL RiskPropertySensitivity<ExpiryPair>);
EXTERN_TEMPLATE(class RISKMGR_DLL RiskPropertySensitivity<ExpiryAndStrike>);
EXTERN_TEMPLATE(class RISKMGR_DLL RiskPropertySensitivity<BoxedInt>);
EXTERN_TEMPLATE(class RISKMGR_DLL RiskPropertySensitivity<Void>::Deriv);
EXTERN_TEMPLATE(class RISKMGR_DLL 
                RiskPropertySensitivity<ExpiryWindow>::Deriv);
EXTERN_TEMPLATE(class RISKMGR_DLL RiskPropertySensitivity<BoxedInt>::Deriv);
EXTERN_TEMPLATE(class RISKMGR_DLL RiskPropertySensitivity<ExpiryPair>::Deriv);
EXTERN_TEMPLATE(class RISKMGR_DLL RiskPropertySensitivity<ExpiryAndStrike>::Deriv);
#endif

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
