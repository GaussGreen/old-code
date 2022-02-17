//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CcmOnlyParameters.hpp
//
//   Description : Per name additional parameters needed for CCM model.
//                 Based on (rather, copy of) the now deprecated CCMParameters.
//                 The difference is that this class derives from 
//                 RationalisedCreditEngineParameters.
//
//   Author      : Antoine Gregoire
//
//   Date        : March 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_CCMONLYPARAMETERS_HPP
#define QLIB_CCMONLYPARAMETERS_HPP

#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/CCMParameters.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

/** Per name additional parameters needed for CcmOnly model */
class MARKET_DLL CcmOnlyParameters: public RationalisedCreditEngineParameters,
                                    public virtual Calibrator::IAdjustable
{
public:
    static CClassConstSP const TYPE;
    
    /** Destructor */
    ~CcmOnlyParameters();

    /** Check the model parameters */
    virtual void validatePop2Object();

    /** Build a CcmOnlyParameters object based on an old-style
        CCMParameters object - which must not be null */
    CcmOnlyParameters(CCMParametersConstSP ccmParam);

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
        
    /** Builds a default CcmOnlyParameters object with all fields defaulted to 0
     * (corresponds to Credit Metrics model) */
    static CcmOnlyParameters* buildDefault();

    friend class Bootstrapper;

private:
    friend class StopCompilerWarning;

    CcmOnlyParameters();
    CcmOnlyParameters(const CcmOnlyParameters& rhs); // don't use
    CcmOnlyParameters& operator=(const CcmOnlyParameters& rhs); // don't use

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    string   name;         //// typically expect this to be the CDS Name

// Yes, this must be private but it is used in CcmOnlyCopulaModel - should be refactored one day
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
typedef smartPtr<CcmOnlyParameters> CcmOnlyParametersSP;
typedef smartConstPtr<CcmOnlyParameters> CcmOnlyParametersConstSP;
#ifndef QLIB_CCMONLYPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CcmOnlyParameters>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CcmOnlyParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CcmOnlyParameters>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CcmOnlyParameters>);
#endif

// Support for wrapper
typedef MarketWrapper<CcmOnlyParameters> CcmOnlyParametersWrapper;
#ifndef QLIB_CCMONLYPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CcmOnlyParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CcmOnlyParameters>);
#endif

DRLIB_END_NAMESPACE

#endif
