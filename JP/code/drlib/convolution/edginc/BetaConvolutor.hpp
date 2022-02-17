//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : BetaConvolutor.hpp
//
//   Description : Convolution Algorithm
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_BETACONVOLUTOR_HPP
#define QLIB_BETACONVOLUTOR_HPP

#include "edginc/IConvolutor.hpp"

DRLIB_BEGIN_NAMESPACE

// Class responsible to perform the beta distribution approximation of the convolution.
// Would probably internally represent the "DiscreteDistribution"s in a different way
// (discretising "values" according to loss unit).
class CONVOLUTION_DLL BetaConvolutor:
    public CObject,
    public virtual IConvolutor
{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** empty Constructor */
    BetaConvolutor();

    /** Constructor */
    BetaConvolutor(  
        double varCutoff,           // maximum variance value in % to be considered to be 0
        int    maxNbIter,           // maximum nb of iterations in the beta function algorithm
        double precision,           // precision in the beta function algorithm
        double flooring);           // flooring value of the internal variables close to 0 

    /* Destructor */
    virtual ~BetaConvolutor();

    /**
    * "Pure" convolution algorithm that returns the whole convoluted distribution
    * [Implements IConvolutor]
    * */
    // Result could be cached if this method is called more than once.
    // not implemented
    virtual IDistribution1DConstSP convolute(IDistribution1DArrayConstSP distributions) const;

    /** [Implements IConvolutor] */
    virtual double convoluteAndIntegrate(
        IDistribution1DArrayConstSP distributions,
        ICreditLossConfigConstSP lossConfig,
        const DateTime& timepoint) const;

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    /**
    * Compute the expected loss of a senior tranche
    * Assuming that L follows a beta distribution
    * with a given expected loss and variance.
    *
    * we use an "incomplete beta" integ close formula 
    */
    double betaTrancheLoss(
               double K,            /* (I) lower strike   */
               double min,          /* (I) maximum value  */
               double max,          /* (I) minimum value  */
               double mean,         /* (I) expected loss  */
               double var) const;   /* (I) loss variance  */

    /* incomplete beta function continued fraction expansion */
    double betacf(double a, double b, double x) const;
    
    /* incomplete beta function */
    double betai(double a, double b, double x) const;

    // beta function algorithm parameters :
    double varCutoff;           // $optional(maximum variance value in % to be considered to be 0)
    int    maxNbIter;           // $optional(maximum nb of iterations in the beta function algorithm)
    double precision;           // $optional(precision in the beta function algorithm)
    double flooring;            // $optional(flooring value of the internal variables close to 0 
                                //           in the beta function algorithm)

    // default values
    static const double VAR_CUTOFF;
    static const int MAX_NB_ITER;
    static const double PRECISION;
    static const double FLOORING;
};

DRLIB_END_NAMESPACE

#endif

