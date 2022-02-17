//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CCMParameters.hpp
//
//   Description : Per name additional parameters needed for CCM model.
//
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//      This class is now deprecated: Use CcmOnlyParameters instead.
//      All non-essential methods in this class have been removed or throw
//      exceptions (note that instances of this class are automatically
//      converted to the new-style parameters in 
//      CreditEngineParameters::convertToNewParamStyle().
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//
//   Author      : Antoine Gregoire
//
//   Date        : March 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_CCM_PARAMETERS_HPP
#define QLIB_CCM_PARAMETERS_HPP

#include "edginc/CreditEngineParameters.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

/** Per name additional parameters needed for CCM model */
class MARKET_DLL CCMParameters: public CreditEngineParameters,
                                public virtual Calibrator::IAdjustable
{
public:
    static CClassConstSP const TYPE;
    
    /** Destructor */
    ~CCMParameters();

    /** Check the model parameters */
    virtual void validatePop2Object();

    /** Pull out the component assets & correlations from the market data */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns the name of this object. This is the name with which
        it is stored in the market data cache and is the name with
        which results (eg tweaks) should be reported against */
    virtual string getName() const;
    
    /** Returns whether the name is modelling stochastic recovery
        - true : stochastic recovery
        - false : fixed recovery.
        The default implementation in this class returns false */
    bool hasStochasticRecovery() const;

    /** Computes dependence survival proba between startDate and endDate */ 
    double dependenceSurvivalProba(
        const DateTime& startDate, const DateTime& endDate) const;

    /** Returns independence factor */
    double getIndependenceFactor() const;

    /**
     * Returns catastrophic recovery factor f such that
     * catastrophic loss = notional * (1 - f * marketRecovery)
     * */
    double getCatastrophicRecoveryFactor() const;
    
    // [TO BE RETIRED]
    double getBetaTweak() const;
    
    // [TO BE RETIRED]
    double getQM() const;

    // [TO BE RETIRED]
    double getRecDispersion() const;
    
    // [TO BE RETIRED]
    double getBetaR() const;
        
    /** Builds a default CCMParameters object with all fields defaulted to 0
     * (corresponds to Credit Metrics model) */
    static CCMParameters* buildDefault();

    friend class Bootstrapper;

private:
    friend class StopCompilerWarning;

    CCMParameters();
    CCMParameters(const CCMParameters& rhs); // don't use
    CCMParameters& operator=(const CCMParameters& rhs); // don't use

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    string   name;         //// typically expect this to be the CDS Name

// Yes, this must be private but it is used in CCMCopulaModel - should be refactored one day
// and use dependenceSurvivalProba(...) method instead of accessing the field directly.
public:
    ICDSParSpreadsWrapper    seniorCurve;  ///<! spread assigned to the dependence copula

private:
    double   qm;           ///<! skewed gaussian copula skew
    double   indepFactor;  ///<! non cata spread assigned to independence
    double   recDispersion;///<! recovery dispersion (impact variance)
    double   recCataFactor;///<! catastrophic recovery factor
    double   betaR;        ///<! recovery correlation
    double   betaTweak;    ///<! relative beta tweak from historical beta
};

// Support for smart pointers
typedef smartPtr<CCMParameters> CCMParametersSP;
typedef smartConstPtr<CCMParameters> CCMParametersConstSP;
#ifndef QLIB_CCMPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CCMParameters>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CCMParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CCMParameters>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CCMParameters>);
#endif

// Support for wrapper
typedef MarketWrapper<CCMParameters> CCMParametersWrapper;
#ifndef QLIB_CCMPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CCMParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CCMParameters>);
#endif

DRLIB_END_NAMESPACE

#endif
