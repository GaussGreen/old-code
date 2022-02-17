//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 11-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CREDITMETRICSDEFAULTSMODEL_HPP
#define QLIB_CREDITMETRICSDEFAULTSMODEL_HPP

#include "edginc/IConditionalDefaultsModel.hpp"

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
class CONVOLUTION_DLL CreditMetricsDefaultsModel:
    public CObject,
    public virtual IConditionalDefaultsModel
{
public:
	/** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Virtual destructor */
    virtual ~CreditMetricsDefaultsModel();

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
    CreditMetricsDefaultsModel(
        IIntegratorSP integrator,
        IMarketFactorModelSP marketFactorModel);

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

private:
    CreditMetricsDefaultsModel();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    
    /** Market factor integrator - used for internal thresholds calibration */
    IIntegratorSP integrator;

    /** Market factor model - defines distribution of the market factor */
    IMarketFactorModelSP marketFactorModel;
};

DECLARE(CreditMetricsDefaultsModel);

DRLIB_END_NAMESPACE

#endif /*CREDITMETRICSDEFAULTSMODEL*/
