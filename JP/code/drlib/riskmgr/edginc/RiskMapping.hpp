/**
 * @file RiskMapping.hpp
 */

#ifndef QLIB_RiskMapping_H
#define QLIB_RiskMapping_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IRiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ExpiryWindow)
FORWARD_DECLARE(ExpiryPair)
FORWARD_DECLARE(ExpiryAndStrike)
FORWARD_DECLARE(IModel)
FORWARD_DECLARE(MarketData)
FORWARD_DECLARE(IInstrumentCollection)
FORWARD_DECLARE(IHypothesis)
FORWARD_DECLARE(RiskMappingMatrix)
FORWARD_DECLARE(RiskMapping)

/**
 * A linear mapping between IRiskProperty's --- transforms tweaks of
 * Black-Scholes delta/vega/... to tweaks of dynamics parameters objects.
 *
 * This class essentially says how to "fake" tweaks to each of a set of
 * Black-Scholes quantities ("IBM 3M vol", "IBM 6M vol", "IBM spot") using
 * tweaks to a parametric path model (e.g. SRMEQ::Vol "IBM 3M smileA1" or
 * VolSVJ "IBM meanReversRate").
 *
 *
 * <H3>Overview</H3>
 *
 * [See also "Risk mapping" in the QLib doc db.]
 *
 * Prices computed using models which are parameterised in terms of objects
 * like like SRMEQ::Vol or VolSVJ look as if they are insensitive to good old
 * fashioned Spot/VolParallel/... --- e.g. Delta and VegaParallel come out
 * "not applicable".  Or non-zero, but misleading, which is even worse.
 *
 * What we [will soon] do about this is outlined below:
 *
 *    -  Dynamics parameters like VolSVJ or SRMEQ::Vol are "marked" maybe once
 *       a day by a batch Calibrator process and written back to Pyramid
 *
 *    -  After calibrating each parameters object, Pyramid calls
 *       CalibratedRiskMappingMatrix, obtaining a RiskMappingMatrix to
 *       go with it:
 *          -  EAS gives us a list of desired Black-Scholes greeks (Delta,
 *             VegaParallel, ...)
 *          -  By pricing the same grid of vanillas using both closed-form and
 *             calibrated parametric models, we figure out how to tweak
 *             the parameter's fields in order to make their prices change in
 *             the same way as they do under each desired Black-Scholes tweak
 *          -  We return a RiskMappingMatrix summarising the relationship
 *             and it gets written down to Pyramid
 *
 *    -  When pricing exotics using parametric models, EAS feeds us back
 *       the calibrated RiskMappingMatrix for each parameter object
 *       involved, as part of the market data
 *
 *    -  They're assembled into a single RiskMapping object (see fromMarket())
 *
 *    -  Greeks written using the new-style "declarative" framework are
 *       passed through this transform:
 *          -  During their evaluation
 *             (by RiskQuantityEvaluator::storeResults()),
 *             every IHypothesis they generate is intercepted
 *             (in HypothesisTree::evaluate())
 *             and passed through mapped()
 *          -  So instead of e.g. performing a (BS) vol tweak, to which
 *             VolSVJ is not sensitive, we make a combination of
 *             tweaks to VolSVJ parameters which achieve the same effect
 *
 *
 * <H3>How to enable risk mapping</H3>
 *
 *    -  Ensure that the greeks you want mapped are written in the
 *       "declarative" sensitivities framework, i.e. are implementations of
 *       IRiskQuantityFactory (which see for an overview).  Delta,
 *       VegaPointwise, VegaParallel and some others already are.
 *
 *    -  Ensure that the dynamics parameter defining the model for
 *       which you want risk mapping implements the IDynamicsParameter
 *       interface.  VolSV, VolSVJ, VolSVCJ, SRMEQ::Vol and SRMFX::Vol
 *       already do.
 *
 *    -  Get Analytics to schedule a batch process which obtains
 *       CalibratedRiskMappingMatrix for your desired greeks whenever
 *       a parameter for your model is calibrated.
 *
 *    -  Get clients to pass that RiskMappingMatrix as market data
 *       along with the dynamics parameter.
 *
 *
 * <H3>Summary of classes involved</H3>
 *
 * <DL>
 * <DT>RiskMappingMatrix        
 *    <DD>How to fake BS tweaks in terms of tweaks to the fields of a particular
 *        parameter object
 *
 * <DT>CalibratedRiskMappingMatrix
 *    <DD>"Addin" which generates RiskMappingMatrix's for a given parameter
 *        object
 *
 * <DT>RiskMapping
 *    <DD>Union of RiskMappingMatrix's: the key method is mapped()
 *
 * <DT>RiskQuantityEvaluator
 *    <DD>Performs calculation of all requested greeks; intercepts tweaks and
 *        passes them through RiskMapping
 * </DL>
 */

class RISKMGR_DLL RiskMapping: public CObject {

    RiskMapping(const RiskMapping& rhs);
    RiskMapping& operator=(const RiskMapping& rhs);
    static IObject* defaultOne();
    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

private:

    mutable RiskMappingMatrixArrayConstSP matrices; // $required

    // a mechanism by which the matrices field may be
    // dynamically extended
    void extendMatrices(RiskMappingMatrixConstSP mtx) const;

public:

    /**
     * A RiskMapping which is effectively the union of a set of
     * RiskMappingMatrix's
     *
     * (Called via fromMarket() by RiskQuantityEvaluator::getMarket().)
     */

    RiskMapping(RiskMappingMatrixArrayConstSP matrices);

    /**
     * An identity RiskMapping, which just maps properties and hypotheses to
     * themselves
     */

    RiskMapping();

    /**
     * An identity RiskMapping, which just maps properties and hypotheses to
     * themselves
     */

    static RiskMappingSP SP();

    /**
     * A RiskMapping assembled from all relevant RiskMappingMatrix's provided
     * in the market
     */

    static RiskMappingSP fromMarket(IModelConstSP model,
                                    IInstrumentCollectionConstSP instruments,
                                    MarketDataConstSP market);

    ~RiskMapping();

    /**
     * The market names in the world which carry the given property, or for
     * which we can map the (Black-Scholes) property to a parametric
     * pseudo-property
     *
     * Looks at the IRiskAxis labelling each column of all the matrices; if its
     * IRiskAxis::abstractProperty() matches, its IRiskAxis::marketDataName()
     * is included.
     *
     * If no matching columns are found, just returns @a property's
     * IAbstractRiskProperty::subjectNames().
     */

    OutputNameArrayConstSP subjectNames(
        IAbstractRiskPropertyConstSP property,
        IObjectConstSP world) const;

    /**
     * The qualifiers available for tweaking a given property in the world, or
     * for which we can map the (Black-Scholes) property to a parametric
     * pseudo-property
     *
     * Looks at the IRiskAxis labelling each column of all the matrices; if
     * its IRiskAxis::abstractProperty() and IRiskAxis::marketDataName()
     * match, and it's an IQualifiedRiskAxis<QUALIFIER>, then its
     * expiry is included.
     *
     * If no matching columns are found, just returns @a property's
     * IRiskProperty<QUALIFIER>::subjectQualifiers().
     */

    template <class Q>
    smartConstPtr<array<smartPtr<Q>, Q> > subjectQualifiers(
        smartConstPtr<IRiskProperty<Q> > property,
        OutputNameConstSP name,
        IObjectConstSP world) const;

    /**
     * A hypothesis, constructed in terms of model parameter tweaks, which has
     * approximately the same effect as a given "Black-Scholes" hypothesis,
     * while leaving other Black-Scholes quantities approximately equal.
     */

    IHypothesisConstSP mapped(IHypothesisConstSP hypothesis) const;

    /**
     * For testing
     */

    void storeRequestResults(CControlSP control, ResultsSP results) const;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
