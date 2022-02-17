/**
 * @file RiskMappingMatrix.hpp
 */

#ifndef QLIB_RiskMappingMatrix_H
#define QLIB_RiskMappingMatrix_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/IRiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CDoubleMatrix)
FORWARD_DECLARE(IHypothesis)
FORWARD_DECLARE(AbstractPropertyTweakHypothesis)
FORWARD_DECLARE(AtomicHypothesis)
FORWARD_DECLARE(RiskAxis)
FORWARD_DECLARE(IRiskAxis)
FORWARD_DECLARE(RiskMappingMatrix)

/**
 * A linear mapping between IHypothesis's --- transforms tweaks of
 * Black-Scholes delta/vega/... to tweaks of a dynamics parameter object.
 *
 * This class essentially says how to "fake" tweaks to each of a set of
 * Black-Scholes quantities ("IBM 3M vol", "IBM 6M vol", "IBM spot") using
 * tweaks to a parametric path model (e.g. SRMEQ::Vol "IBM 3M smileA1" or
 * VolSVJ "IBM meanReversRate").
 *
 * It's returned from CalibratedRiskMappingMatrix during a batch "risk mapping
 * calibration", stored in Pyramid, and passed back to us in the market data
 * along with the model parameters.
 *
 * Each RiskMappingMatrix captures the BS-to-param tweak mapping for a single
 * model parameter object---i.e. for a single name.  To handle multi-factor
 * models like SRM3, we combine all the relevant RiskMappingMatrix's in the
 * market data into an overall RiskMapping.
 *
 * Then when we're asked to compute Black-Scholes greeks, we can use the
 * RiskMapping to create scenarios under which (say) IBM 3M vol changes "other
 * things being equal", in spite of the fact that the model is not
 * independently parameterised in terms of that quantity.  This allows us to
 * compute a vega directly comparable with those obtained using non-parametric
 * pricers.
 *
 * The key method is mapped().
 *
 * See RiskMapping for more information.
 */

// $public
class RISKMGR_DLL RiskMappingMatrix: public MarketObject {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

private:

    static IObject* defaultRiskMappingMatrix();
    RiskMappingMatrix();
    void validate();
    void validatePop2Object();

    /**
     * Fields as stored in Pyramid: our "frozen" representation.  See RiskAxis.
     */

    // $required(Unique (but otherwise not meaningful\) market data name)
    string name;
    RiskAxisArrayConstSP inputBasis;                   // $required
    CDoubleArraySP inputDistanceScaling;               // $required
    RiskAxisArrayConstSP outputBasis;                  // $required
    CDoubleMatrixSP coefficients;                      // $required
	CDoubleMatrixSP coefficients2;                     // $optional
    CDoubleArraySP errorBars;                          // $optional

    /**
     * Unfrozen representation of basis
     */

    IRiskAxisArrayConstSP _inputBasis;                 // $transient
    IRiskAxisArrayConstSP _outputBasis;                // $transient

    friend class RiskMappingMatrixCache; // see RiskMapping.cpp

    /**
     * Cache IRiskProperty's along which bases are aligned
     */

    IAbstractRiskPropertyArraySP _inputProperties;     // $transient

public:

    /**
     * Constructor.
     *
     * @param name          The market name for the RiskMappingMatrix
     *
     * @param inputBasis    The (Black-Scholes) risk axes to which the
     *                      columns of the matrix correspond (typically
     *                      PropertyRiskAxis)
     *
     * @param inputDistanceScaling  For each of the inputBasis risk axes:
     *                      the ratio between the nominal shift coefficient
     *                      (shift size) requested, and the distance by which
     *                      the resulting hypothesis says it changes the world.
     *                      E.g.  Equity::sensShift(const PropertyTweak<Spot>&)
     *                      treats the the shift size as specifying a
     *                      proportional tweak, but reports the
     *                      TweakOutcome::distance() in absolute terms, so that
     *                      Delta ends up being expressed in absolute units.
     *                      We need to know this so that we can report
     *                      (something like) the same distance when we fake a
     *                      B-S spot shift using dynamics parameter tweaks.
     * 
     * @param outputBasis   The risk axes to which the rows of the matrix
     *                      correspond (typically ScalarFieldRiskAxis
     *                      or ExpiryFieldRiskAxis)
     *
     * @param coefficients  The matrix: inputBasis->size() columns and
     *                      outputBasis->size() rows
     *
     * @param errorBars     errorBars[g] gives a rough idea of the error in
     *                      a first derivative along inputBasis[g] estimated by
     *                      risk-mapping using this matrix.  The error can
     *                      arise from two sources: inability of the model to
     *                      express all necessary scenarios (a stationary VolSV
     *                      can't be made to mimic a vol tweak to a particular
     *                      expiry), and breakdown of the linearity assumptions
     *                      on which the risk mapping methodology rests.  As
     *                      provided by CalibratedRiskMappingMatrix, the
     *                      error bars capture the former problem but not the
     *                      latter.
     */

    RiskMappingMatrix(const string& name,
                      IRiskAxisArrayConstSP inputBasis,
                      CDoubleArrayConstSP inputDistanceScaling,
                      IRiskAxisArrayConstSP outputBasis,
                      const CDoubleMatrix& coefficients,
                      CDoubleArrayConstSP errorBars);

	RiskMappingMatrix(const string& name,
                      IRiskAxisArrayConstSP inputBasis,
                      CDoubleArrayConstSP inputDistanceScaling,
                      IRiskAxisArrayConstSP outputBasis,
                      const CDoubleMatrix& coefficients,
					  const CDoubleMatrix& coefficients2,
                      CDoubleArrayConstSP errorBars);

    /**
     * Constructor returning SP.
     */

    static RiskMappingMatrixSP SP(const string& name,
                                  IRiskAxisArrayConstSP inputBasis,
                                  CDoubleArrayConstSP inputDistanceScaling,
                                  IRiskAxisArrayConstSP outputBasis,
                                  const CDoubleMatrix& coefficients,
                                  CDoubleArrayConstSP errorBars);

	static RiskMappingMatrixSP SP(const string& name,
                                  IRiskAxisArrayConstSP inputBasis,
                                  CDoubleArrayConstSP inputDistanceScaling,
                                  IRiskAxisArrayConstSP outputBasis,
                                  const CDoubleMatrix& coefficients,
								  const CDoubleMatrix& coefficients2,
                                  CDoubleArrayConstSP errorBars);

    ~RiskMappingMatrix();

    /**
     * Append to @a names the market data names on which @a property can be
     * faked using this matrix
     */

    void appendMappedNames(
        IAbstractRiskPropertyConstSP property,
        OutputNameArraySP names) const;

    /**
     * Append to @a qualifiers the qualifiers for which @a property can be
     * faked on the market data name @a name using this matrix
     */

    template <class Q>
    bool appendMappedQualifiers(
        smartConstPtr<IRiskProperty<Q> > property,
        OutputNameConstSP name,
        smartPtr<array<smartPtr<Q>, Q> > qualifiers) const;


    /**
     * A hypothesis, constructed in terms of model parameter tweaks, which has
     * approximately the same effect as a given "Black-Scholes" hypothesis,
     * while leaving other Black-Scholes quantities approximately equal
     */

    IHypothesisConstSP mapped(IHypothesisConstSP hypothesis) const;

    /**
     * The market data name for the matrix, as specified via
     * CalibratedRiskMappingMatrix::name
     */

    string getName() const;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
