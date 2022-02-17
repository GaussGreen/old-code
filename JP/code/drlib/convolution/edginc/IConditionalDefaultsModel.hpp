//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 11-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_ICONDITIONALDEFAULTSMODEL_HPP
#define QLIB_ICONDITIONALDEFAULTSMODEL_HPP

#include "edginc/IMarketFactorValue.hpp"
#include "edginc/IIntegrator.hpp"
#include "edginc/IMarketFactorModel.hpp"
#include "edginc/CreditEngineParameters.hpp"
#include "edginc/ICondLossDistributionsGen.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IConditionalDefaultsModel);

/**
 * Interface for "conditional defaults" models:
 * 
 * For a name i, define:
 * X_i(T) = f(Z(T), e_i, p(T))
 * 
 * where
 *   Z(T) = Market factor (any shape)
 *   e_i  = Random variable independent of Z(T)
 *   p(T) = Model parameter
 *   T    = Time
 * 
 * Defaults happens when X_i(T) < c_i(T) (threshold).
 * 
 * Method conditionalSurvivalProbas computes P(X_i(T) < c_i(T) | Z=z).
 * 
 * Method calibrateThreshold computes c_i such that
 * integral(p(X_i(T) < c_i(T) | Z=z)dz) = defaultProba.
 * 
 * Note that in the general case, calibrateThreshold method may need
 * an IMarketFactorModel and an IIntegrator.
 * 
 * Possible implementations would be:
 * - Credit Metrics
 * - CCM
 * - RFL
 * - Stochastic correlation
 * - Generic copula products
 * - ...
 * 
 * */
class CONVOLUTION_DLL IConditionalDefaultsModel: public virtual IObject {
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /**
     * Give a chance to do some market factor independent
     * initialisation (eg: Calibrates a threshold c_i(T) such that
     * integral(p(X_i(T) < c_i(T) | Z=z)dz) = defaultProba).
     * Creates an IKey that will be used to compute loss distributions conditional
     * on some market factor value.
     * */
    virtual ICondLossDistributionsGenKeySP initialise(
        double defaultProba,
        double expectedLoss,
        double notional,
        CreditEngineParametersConstSP modelParameters,
        const DateTime& startDate,
        const DateTime& endDate) const = 0;

    /**
     * Performs the integration of whatever conditional function is passed in, 
     * over the market factor, using relevant integrator and market factor "model".
     * Input "ICondLossDistributionsGenKeyArrayConstSP condELKeys" can be used to
     * optimise integration (see for example implementation in CompositeCopulaDefaultsModel).
     * Typical the conditional functions are the  expected loss (function "convolutedCondEL")
     * or the survival probability.
     * 
     * NB: Market factor model and associated integrator form part of the "IConditionalDefaultsModel"
     * */
    virtual double integrateCondFunction(
        const MFunctionND* condFunction,
        ICondLossDistributionsGenKeyArrayConstSP condKeys,
        const DateTime& time) const = 0; 

    /** Returns whether the "integrateCondFunction" or the underlying
        marketFactorModel are time dependant */
    virtual bool isCondFunctionIntegrationTimeDependent() const { return false; }

    /**
     * Returns the market factor dimension. Eg:
     * - 1 for 1D Gaussian Credit Metrics
     * - 2 for Composite Copula Model using 1D Credit Metrics
     * */
    virtual int marketFactorDimension() const = 0;
        
    /**
     * Returns the type of CreditEngineParameters expected
     * by this IConditionalDefaultsModel
     * */
    virtual CClassConstSP engineParamsType() const = 0;        

    /**
      Returns measure of "non-desirability" of model parameters.
      This function will be called by the objective function class used in calibration by optimization.
      Note that the function name may not be appropriate in all contexts; it is aimed at the situation
      where the calibration problem is underdetermined and some non-linear smoothing requirement is added.
    */
    virtual double ParameterRoughness() const  { return 0.0; }     
};

DRLIB_END_NAMESPACE

#endif /*QLIB_ICONDITIONALDEFAULTSMODEL_HPP*/
