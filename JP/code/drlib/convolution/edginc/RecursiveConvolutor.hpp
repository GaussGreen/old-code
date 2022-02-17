//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : RecursiveConvolutor.hpp
//
//   Description : Convolution Algorithm
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_RECURSIVECONVOLUTOR_HPP
#define QLIB_RECURSIVECONVOLUTOR_HPP

#include "edginc/IConvolutor.hpp"

DRLIB_BEGIN_NAMESPACE

// Class responsible to perform the recursive convolution.
// Would probably internally represent the "DiscreteDistribution"s in a different way
// (discretising "values" according to loss unit).
class CONVOLUTION_DLL RecursiveConvolutor:
    public CObject,
    public virtual IConvolutor
{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Constructor */
    RecursiveConvolutor(  
                            int                         maxNbSlice,         // max nb of slices for loss distribution
                            int                         nbSubdivision,      // nb by which to divide the GCD of notionals
                            double                      lossUnitOverride,   // loss unit override
                            int                         discretisationType); // discretisation type

    /* Destructor */
    virtual ~RecursiveConvolutor();

    /**
     * "Pure" convolution algorithm that returns the whole convoluted distribution
     * [Implements IConvolutor]
     * */
    // Result could be cached if this method is called more than once.
    virtual IDistribution1DConstSP convolute(IDistribution1DArrayConstSP distributions) const;

    /** [Implements IConvolutor] */
    virtual double convoluteAndIntegrate(
        IDistribution1DArrayConstSP distributions,
        ICreditLossConfigConstSP lossConfig,
        const DateTime& timepoint) const;
    
    // Computes the loss unit
    double lossUnit(IDistribution1DArrayConstSP distributions) const;

private:
    /**
     * "Pure" convolution algorithm that returns the whole convoluted distribution
     * accurately up to a cutoff after which all losses are aggregated into one loss amount
     * */
    virtual IDistribution1DConstSP convolute(IDistribution1DArrayConstSP distributions, double *cutoff) const;

    /** empty Constructor */
    RecursiveConvolutor();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
        
    // Loss units parameters :
    int maxNbSlice;             // $optional(max nb of slices for loss distribution)
    int nbSubdivision;          // $optional(nb by which to divide the GCD of notionals)
    double lossUnitOverride;    // $optional(loss unit override)
    int    discretisationType;  // $optional(discretisation flag)
    
    // default values
    static const int MAX_NB_SLICE;
    static const int NB_SUBDIVISION;
    static const double LOSS_UNIT_OVERRIDE;
    static const int DISCRETISATION_TYPE;
};

DRLIB_END_NAMESPACE

#endif

