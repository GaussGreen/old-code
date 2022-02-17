//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 11-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_RFLMIXTUREDEFAULTSMODEL_HPP
#define QLIB_RFLMIXTUREDEFAULTSMODEL_HPP

#include "edginc/IConditionalDefaultsModel.hpp"
#include ext_hash_map // used to cache calibrated thresholds
#include "edginc/RflMixtureParameters.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * 1D "Credit Metrics" conditional defaults model:
 * 
 * For a name i, define:
 * X_i = b_i.Z + sqrt(1-b_i.b_i).e_i
 * 
 * where
 *   Z   = Market factor (any 1D shape)
 *   e_i = iid N(0,1) variables independent of Z
 *   b_i = "beta" (model parameter)
 * 
 * Defaults happens when X_i < c_i(T) (threshold).
 * NB: X_i doesn't depend on T.
 * 
 * Method conditionalSurvivalProbas computes P(X_i(T) < c_i(T) | Z=z).
 * 
 * Method calibrateThreshold computes c_i such that
 * integral(p(X_i(T) < c_i(T) | Z=z)dz) = defaultProba.
 * 
 * */
class CONVOLUTION_DLL RFLMixtureDefaultsModel:
    public CObject,
    public virtual IConditionalDefaultsModel
{
public:
	/** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Override clone() to copy the local cache */
    virtual IObject* clone() const;

    /** [Implements IConditionalDefaultsModel] */
    virtual ICondLossDistributionsGenKeySP initialise(
        double defaultProba,
        double expectedLoss,
        double notional,
        CreditEngineParametersConstSP modelParameters,
        const DateTime& startDate,
        const DateTime& endDate) const;

    /** [Implements IConditionalDefaultsModel] */
    virtual double integrateCondFunction(
        const MFunctionND* convolutedCond,
        ICondLossDistributionsGenKeyArrayConstSP condKeys,
        const DateTime& time) const; 

    /** [Implements IConditionalDefaultsModel] */
    int marketFactorDimension() const; 

    /** [Implements IConditionalDefaultsModel] */
    CClassConstSP engineParamsType() const;

    /** Public constructor */
    RFLMixtureDefaultsModel(double c,
                     double w,
                     double factorShift = 0.0,
                     double lowerBound = -7.5,
                     double upperBound = 7.5,
                     unsigned int numNodes = 101);

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    // Set value of w
    const double& SetW(double w) { return m_w = w; }

private:
// foundation
    RFLMixtureDefaultsModel();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    
// types
    /** Key for the thresholds cache */
    struct CONVOLUTION_DLL Key
    {
        Key(double defaultProba,
            IObjectConstSP rflParams);
            
        bool operator==(const Key& key) const;
        
        double defaultProba;
        IObjectConstSP rflParams;
    };
    
    /** hash class for use in hashmap's that use Key */
    struct CONVOLUTION_DLL KeyHash {
        size_t operator()(const Key& key) const;
    };
    
    typedef hash_map<Key, double, KeyHash> ThresholdsCache;
    typedef pair<Key, double> MapEntry; // (key, value) pair 
    
// data
    // Mixture parameters
    double m_c, m_w;
    // Factor shift
    double m_factorShift;

    // Quadrature parameters (to be used with 'GridIntegrator')
    double m_lowerBound, m_upperBound;
    int m_numNodes;

    /**
     * Cache calibrated thresholds so we only recompute them when necessary.
     * Cache improves performance:
     * - when calibrating RFL parameters with more than one
     *   tranche on the same portfolio 
     * - when computing greeks (eg: delta pointwise)
     * 
     * It does not alter performance for simple pricing.
     * */
    mutable ThresholdsCache barriersCache; // $unregistered
};

DECLARE(RFLMixtureDefaultsModel);

DRLIB_END_NAMESPACE

#endif /*RFLMIXTUREDEFAULTSMODEL*/
