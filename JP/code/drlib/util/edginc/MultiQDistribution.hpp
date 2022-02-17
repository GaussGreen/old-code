//----------------------------------------------------------------------------
//
//   Group       : Credit QRD
//
//   Filename    : MultiQDistribution.hpp
//
//   Description : Helper class for calibrating Multi-Q Distributions
//                 And calculating option prices and greeks
//
//   Author      : Charles Morcom
//
//   Date        : August 5, 2005
//
//
//----------------------------------------------------------------------------

#ifndef MULTIQDISTRIBUTION_HPP
#define MULTIQDISTRIBUTION_HPP

#include "edginc/AtomicArray.hpp"
#include "edginc/Object.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IDistribution1D.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(MultiQDistribution);
FORWARD_DECLARE(DiscreteDistribution);

/*===========================================================================*/
/**MultiQDistribution
 * A class to represent a multi-Q distribution. This is a mapping of a normal
 * random variable, X ~ N(muQ, sigQ) to a new random variable,
 * Y = F(K(CH(X)-1)+1),
 * where H(X) piece-wise equal to (exp(qX)-1)/q. 
 * C = E[H(X)], so that F is the forward value of Y. K is determined by
 * the option prices.
 * Vanilla option prices are easily calculable in piecewise closed-form and
 * smile calibration is (relatively) straightforward.
 */
class UTIL_DLL MultiQDistribution : public IDistribution1D,
                                    public CObject {
public:
    static CClassConstSP const TYPE; // for reflection API
    static void load(CClassSP& clazz);

    virtual ~MultiQDistribution();
    // default assignment and copy constructors are OK
    
    /*=======================================================================*/
    /**BSQIntegral
     * Evaluates the integral of exp(qx).phi(x, mu, sigma) on [a,b], where
     * phi is the density of a normal RV with mean mu, standard deviation
     * sigma.
     * BSQIntegral(a,b,mu,sigma,q) = exp(q.mu + (q.sigma)^2/2)[N(B-q.sigma)-N(A-q.sigma)]
     * where B=(b-mu)/sigma and A=(a-mu)/sigma
     */
    static double BSQIntegral(
        double a, 
        double b, 
        double mu, 
        double sigma, 
        double q);

    /*=======================================================================*/
    /**MultiQDistribution
     * Create a distribution directly from a specified ATM option, and
     * a set of q values and deltas.
     * The number of deltas must be one bigger than the number of qs. The
     * middle delta is used as the start point for the left qs and the right qs.
     * The deltas are the cumulative gaussian space probability strikes: 0.5
     * corresponds to the forward.
     */
    MultiQDistribution(
        double              optionTimeToExpiry,
        double              forward,
        double              atmBSVolatility,
        const DoubleArray&  inputQs,
        const DoubleArray&  inputDeltas /**<Must be ordered.*/
        );

    MultiQDistribution(double              optionTimeToExpiry,
                       double              forward,
                       double              atmBSVolatility,
                       const DoubleArray&  inputQs,
                       const DoubleArray&  inputDeltas, /**<Must be ordered.*/
                       const string&       inputType
                       );

  //Create a log normal multiQ distribution
  static MultiQDistributionSP createLogNormalMq(double              optionTimeToExpiry,
						double              forward,
						double              atmBSVolatility);
  
    
    //-----------------------------------------------------------------------
    //  IDistribution1D methods
    //-----------------------------------------------------------------------
    /**Probability density function of the mapped value: P(y<=Y<=y+dy) = mq.pdf(y) dy.
       Fixed - (MC 3/17/2006)*/
    double pdf(double y) const;
    
    /**Cumulative density function of the mapped value. P(Y<=y) = mq.cdf(y) */
    double cdf(double y) const;

    double expectation() const;

    double expectation(const Function1DDouble& payoff) const;

    double variance() const;

    DiscreteDistributionConstSP discretise() const;

    double mgf(double z, int n) const;

    //-----------------------------------------------------------------------
    //  MultiQDistribution methods
    //-----------------------------------------------------------------------
    /**Returns the expected value of the distribution - probably the
       forward price, or spread, depending on what you are using this
       for. */
    double forward() const;
    
    /**Returns the value that a gaussian(muQ, sigQ) 
       x-value is mapped to by the Q-distribution. This is H(x).*/
    double qMap(double xVal) const;
	/**Returns the value sigQ, the stdev of the gaussian variable*/
	double getNormalVol() const;
    /**Computes the price of a vanilla call or put with the M-Q distribution*/
    double vanillaOptionPrice (
        bool isCall, /**<true if call, else put*/
        double         strike  /**<Option strike price */
    ) const;

    /**Computes the price of a binary call or put with the M-Q distribution*/
    double binaryOptionPrice (
        bool isCall, /**<true if call, else put*/
        double         strike  /**<Option strike price */
    ) const;

protected:
    // INTERNALLY CALIBRATED MEMBERS NOT VISIBLE TO REFLECTION INTERFACE
    mutable double f;       /**<Foward of distribution (i.e. expectation)            $unregistered */
    mutable double muQ;     /**<Mean of normal driver RV                             $unregistered */
    mutable double sigQ;    /**<St. Dev. of normal driver RV (and ATM sigma sqrt(T)) $unregistered */
    mutable double atmPx;   /**<Forward price of ATM B-S call or put                 $unregistered */
    mutable double atmVol;  /**<B-S ATM volatility                                   $unregistered */
    mutable double atmT;    /**<Time to expiry of option                             $unregistered */
    mutable double C;       /**<Forward calibration affine coefficient               $unregistered */
    mutable double K;       /**<Option calibration affine coefficient                $unregistered */
    mutable DoubleArray qL; /**<Left Q-coefficients                                  $unregistered */
    mutable DoubleArray qR; /**<Right Q-coefficients                                 $unregistered */
    mutable DoubleArray dL; /**<Left slope D-coefficients                            $unregistered */
    mutable DoubleArray dR; /**<Right slope D-coefficients                           $unregistered */
    mutable DoubleArray kL; /**<Left continuity k-coefficients                       $unregistered */
    mutable DoubleArray kR; /**<Right continuity k-coefficients                      $unregistered */
    mutable DoubleArray xL; /**<Left x-range coefficients                            $unregistered */
    mutable DoubleArray xR; /**<Right x-range coefficients                           $unregistered */

    // INPUT ARGUMENTS FOR REFLECTION INTERFACE
    double inputTimeToExpiry;
    double inputForward;
    double inputATMVolatility;
    string inputCalibrationType;
    DoubleArray inputVolParameters;
    DoubleArray inputStrikeBoundaryParameters;
    
private:
    /**Default constructor for reflection use only */
    MultiQDistribution();
    /**TRUE if calibrateFromInputs has been called*/
    mutable bool calibrated; // $unregistered
    /**TRUE if have tried to calibrate, but failed - signals not to try again*/
    mutable bool calibrationFailed; // $unregistered
    /**TRUE if publi functions should not call calibration: set during
        initialization to prevent infinite loops */
    mutable bool dontCalibrate; // $unregistered
    /**Called to set internal parameters given input parameters */
    void calibrateFromInputs() const;
    /**Calibrate from explicit qs and strike deltas. */
    void calibrateExplicitQAndDeltaInputs() const;
    /**Calibrate distribution from BS vols and strike ratios. */
    void calibrateBSImpliedVolInputs() const;
    /**Set K and C given all other parameters */
    void calibrateAffineAdjustments() const;
    /**Calculate the forward - used for calibration iteration */
    double calculateExpectedH() const;
    /**Ensure all qs are bounded away from zero */
    void ensureNonZeroQs() const;
    static IObject* defaultMultiQDistribution();
    /**Find a point's location in gaussian x-space.*/
    void locateInGaussianSpace(
        double val, /**<Value to find*/
        bool isX, /**<Is val a forward space value, or an x-space value?*/
        bool* isLeft, /**<Returns true if left-side, else false*/
        int* nQ, /**<Number of qs on side*/
        int* idx, /**<(*x)[idx] is just closer to 0 than val */
        double* hVal /**<Value of H(x).*/
        ) const;

    friend class MultiQDistributionAddin;
    friend class MultiQDistributionAddin2;
    friend class MapQsToPriceDiffs;
};

DRLIB_END_NAMESPACE

#endif
