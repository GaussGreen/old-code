/**
 * @file CrossRiskPropertySensitivity.hpp
 */

#ifndef QLIB_CrossRiskPropertySensitivity_H
#define QLIB_CrossRiskPropertySensitivity_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IRiskProperty.hpp"
#include "edginc/RiskQuantityFactorySensitivity.hpp"
#include "edginc/IPerNameSensitivity.hpp"
#include "edginc/ICrossDerivative.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IScalarDerivative)
FORWARD_DECLARE(IResultsFunction)

/**
 * A cross-gamma calculated by tweaking pairs of elements of some "risk
 * property", on each of the "names" in the market in turn
 *
 * This is part of the "declarative" sensitivities framework: see
 * IRiskQuantityFactory for an overview.  It's the base class for cross gammas.
 * The key method is nameRiskQuantities().
 *
 * Why's this class a template?  Because the abstract IRiskProperty interface
 * under which we store our 'property0' and 'property1' fields is parameterised
 * by the type of "qualifier" it requires (ExpiryWindow for term-structured
 * properties like VolPointwise, Void for scalars like Spot).
 *
 *
 * <H3>Implementations</H3>
 *
 * Implementations at the time of writing include
 *
 *    -  VSCurveCrossGamma.cpp
 *
 * CrossRiskPropertySensitivity corresponds roughly to CrossGamma in the
 * existing Sensitivity framework, but also subsumes greeks indexed by
 * ExpiryWindow.  We could implement array-valued greeks indexed by plain
 * integers, and even greeks indexed by ExpiryPair, very easily.
 *
 * You'll notice that we don't implement ITwoSidedDeriv --- so we don't
 * play right with "quick greeks".  Ideally we would simplify the quick
 * greeks framework by having it look at the IHypothesis's w.r.t. which
 * the world has changed, but for now ...
 */

template <class QUALIFIER0, class QUALIFIER1>
class CrossRiskPropertySensitivity: public RiskQuantityFactorySensitivity,
                                    public virtual IPerNameSensitivity {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

    typedef QUALIFIER0 Qualifier0;
    DECLARE(Qualifier0)
    typedef QUALIFIER1 Qualifier1;
    DECLARE(Qualifier1)

protected:

    void ensureDeriv() const;
    mutable IResultsFunctionConstSP _derivand;                     // $transient
    mutable smartConstPtr<IRiskProperty<QUALIFIER0> > _property0;  // $transient
    const double shiftSize;                                        // $required
    mutable smartConstPtr<IRiskProperty<QUALIFIER1> > _property1;  // $transient
    const double shiftSize1;                                       // $optional
    mutable ICrossDerivativeConstSP _derivative;                   // $transient
    mutable double _sensitivityUnit;                               // $transient

    IResultsFunctionConstSP derivand() const;
    smartConstPtr<IRiskProperty<QUALIFIER0> > property0() const;
    smartConstPtr<IRiskProperty<QUALIFIER1> > property1() const;
    ICrossDerivativeConstSP derivative() const;
    double sensitivityUnit() const;

    /**
     * Constructor (part 1)
     *
     * CrossRiskPropertySensitivity::deriv() is conceptually the "other half" of
     * this constructor.  We have to split it up because the Deriv you specify
     * may depend on reflection FIELDs which haven't been filled in at
     * construct time.
     *
     * @param type             Reflection type of final subclass
     *
     * @param outputName       Packet name under which quantities computed for
     *                         this greek will be stored in the Results
     *                         (e.g. "VSCURVE_CROSS_GAMMA" for
     *                         VSCurveCrossGamma.cpp)
     *
     * @param shiftSize0       How big to make the @a property0 tweaks
     *                         used to estimate @a derivative
     *
     * @param shiftSize1       How big to make the @a property1 tweaks
     *                         used to estimate @a derivative
     */

    CrossRiskPropertySensitivity(CClassConstSP type,
                                 double shiftSize0,
                                 double shiftSize1,
                                 const string& outputName);

    /**
     * Specification for the derivative to be calculated by a
     * RiskPropertySensitivity
     *
     * This what subclasses must return from deriv(), which is conceptually the
     * "second half" of the constructor RiskPropertySensitivity().
     */

    struct Deriv {
        IResultsFunctionConstSP derivand;
        smartConstPtr<IRiskProperty<QUALIFIER0> > property0;
        smartConstPtr<IRiskProperty<QUALIFIER1> > property1;
        ICrossDerivativeConstSP derivative;
        double sensitivityUnit;

        /**
         * Constructor
         *
         * These parameters are really the "second half" of
         * CrossRiskPropertySensitivity().
         *
         * @param derivand         The quantity whose cross-gamma we want
         *                         to estimate: generally fair value (i.e.
         *                         IResultsFunction::price()) but it could be
         *                         NAKED_BOND_PRICE or anything
         * 
         * @param property0        The first IRiskProperty with respect to
         *                         which the @a derivative is to be estimated
         *
         * @param property1        The second IRiskProperty with respect to
         *                         which the @a derivative is to be estimated
         *
         * @param derivative       The type of derivative to estimate, e.g.
         *                         ICrossDerivative::cross()
         *
         * @param sensitivityUnit  Scaling factor applied to computed derivatives
         *                         before they're placed in the Results.  Think
         *                         about multiplying your scaling factor
         *                         for property0 by that for property1.
         */
        /* inlined (ie definition in header file) as haven't worked
         * out how to get MS compiler to export this when building
         * separate dlls. */
        Deriv(IResultsFunctionConstSP derivand,
              smartConstPtr<IRiskProperty<QUALIFIER0> > property0,
              smartConstPtr<IRiskProperty<QUALIFIER1> > property1,
              ICrossDerivativeConstSP derivative,
              double sensitivityUnit);
    };

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

    /**
     * The risk quantities which this CrossRiskPropertySensitivity says
     * should be computed for a given base-state world
     *
     * Called by RiskQuantityFactorySensitivity::riskQuantities() to implement
     * IRiskQuantityFactory::riskQuantities().
     */

    virtual NamedRiskQuantityArraySP nameRiskQuantities(
        MultiTweakGroupConstSP world, RiskMappingConstSP riskMapping) const;

public:

    ~CrossRiskPropertySensitivity();

    bool discreteShift() const;

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
     * Constructor returning SP
     */

    static smartPtr<CrossRiskPropertySensitivity> SP(
        const string& outputName,
        IResultsFunctionConstSP derivand,
        smartConstPtr<IRiskProperty<QUALIFIER0> > property0,
        smartConstPtr<IRiskProperty<QUALIFIER1> > property1,
        ICrossDerivativeConstSP derivative,
        double sensitivityUnit,
        double shiftSize0, double shiftSize1);
};

#ifndef QLIB_CROSSRISKPROPERTYSENSITIVITY_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL CrossRiskPropertySensitivity<Void _COMMA_ Void>);
EXTERN_TEMPLATE(class RISKMGR_DLL 
                CrossRiskPropertySensitivity<ExpiryWindow _COMMA_ ExpiryWindow>);
EXTERN_TEMPLATE(struct RISKMGR_DLL CrossRiskPropertySensitivity<Void _COMMA_ Void>::Deriv);

EXTERN_TEMPLATE(struct RISKMGR_DLL 
                CrossRiskPropertySensitivity<ExpiryWindow _COMMA_ ExpiryWindow>::Deriv);
#endif

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
