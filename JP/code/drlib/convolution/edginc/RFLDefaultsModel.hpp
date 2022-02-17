//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 11-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_RFLDEFAULTSMODEL_HPP
#define QLIB_RFLDEFAULTSMODEL_HPP

#include "edginc/IConditionalDefaultsModel.hpp"
#include ext_hash_map // used to cache calibrated thresholds
#include "edginc/CmRflParameters.hpp"

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
class CONVOLUTION_DLL RFLDefaultsModel:
    public CObject,
    public virtual IConditionalDefaultsModel
{
public:
	/** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Virtual destructor */
    virtual ~RFLDefaultsModel();
    
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
        const MFunctionND* condFunction,
        ICondLossDistributionsGenKeyArrayConstSP condKeys,
        const DateTime& time) const; 

    /** [Implements IConditionalDefaultsModel] */
    virtual bool isCondFunctionIntegrationTimeDependent() const;

    /** [Implements IConditionalDefaultsModel] */
    virtual int marketFactorDimension() const;

    /** [Implements IConditionalDefaultsModel] */
    virtual CClassConstSP engineParamsType() const;

    /** Public constructor */
    RFLDefaultsModel(
        IIntegratorSP integrator,
        IMarketFactorModelSP marketFactorModel);

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

private:
    RFLDefaultsModel();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    
    /** Market factor integrator - used for internal thresholds calibration */
    IIntegratorSP integrator;

    /** Market factor model - defines distribution of the market factor */
    IMarketFactorModelSP marketFactorModel;
    
    /** Key for the thresholds cache */
    struct CONVOLUTION_DLL Key
    {
        Key(double defaultProba,
            double betaHist,
            IObjectConstSP rflParams);
            
        bool operator==(const Key& key) const;
        
        double defaultProba;
        double betaHist;
        IObjectConstSP rflParams;
    };
    
    /** hash class for use in hashmap's that use Key */
    struct CONVOLUTION_DLL KeyHash {
        size_t operator()(const Key& key) const;
    };
    
    typedef hash_map<Key, double, KeyHash> ThresholdsCache;
    typedef pair<Key, double> MapEntry; // (key, value) pair 
    
    /**
     * Cache calibrated thresholds so we only recompute them when necessary.
     * Cache improves performance:
     * - when calibrating RFL parameters with more than one
     *   tranche on the same portfolio 
     * - when computing greeks (eg: delta pointwise)
     * 
     * It does not alter performance for simple pricing.
     * */
    mutable ThresholdsCache thresholdsCache; // $unregistered
};

DECLARE(RFLDefaultsModel);

DRLIB_END_NAMESPACE

#endif /*RFLDEFAULTSMODEL*/
