//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Description : Container class for Credit Metrics + CMM + RFL 
//                 per name model params 
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CM_CMM_RFL_PARAMETERS_HPP
#define QLIB_CM_CMM_RFL_PARAMETERS_HPP

#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/CmOnlyParameters.hpp"
#include "edginc/CcmOnlyParameters.hpp"
#include "edginc/RflOnlyParameters.hpp"

DRLIB_BEGIN_NAMESPACE

/** Contains parameters for Credit Metrics + Base Correlation.
 * Represents 'model oriented' parameters for a given name. */
class MARKET_DLL CmCcmRflParameters : public RationalisedCreditEngineParameters {

public:
    /** TYPE (for reflection) */        
    static CClassConstSP const TYPE;

    /** Destructor */
    ~CmCcmRflParameters();

    /** Build a CmCcmRflParameters object out of the 
        individual components */
    CmCcmRflParameters(CmOnlyParametersConstSP  cmParam,
                       CcmOnlyParametersConstSP ccmParam,
                       RflOnlyParametersConstSP rflParam);

    /** Populates the object with the market data that this object needs */
    virtual void getMarket(const IModel* model, const MarketData* market);
    
    /** Returns the name of this object */
    virtual string getName() const;
    
    /** Returns the engine parameters of the specified type, if any,
        in this CreditEngineParameters object. The parameter specifies what
        type of engine parameters are required. Null is returned if this 
        object is not of the required type */
    virtual CreditEngineParametersConstSP getEngineParams(
        CClassConstSP engineParamsType) const;
    
private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    CmCcmRflParameters();
    
    /** Default constructor */
    static IObject* defaultConstructor();

    // ----------------
    // FIELDS
    // ----------------   
    /** Name of this market object */
    string name;
    
    /** Parameters for Credit Metrics model */
    CmOnlyParametersWrapper cmParameters;

    /** Parameters for CCM model */
    CcmOnlyParametersWrapper ccmParameters;

    /** Parameters for RFL model */
    RflOnlyParametersWrapper rflParameters;
};

typedef smartPtr<CmCcmRflParameters> CmCcmRflParametersSP;
typedef smartConstPtr<CmCcmRflParameters> CmCcmRflParametersConstSP;

// Support for wrapper
typedef MarketWrapper<CmCcmRflParameters> 
    CmCcmRflParametersWrapper;

#ifndef QLIB_CMCCMRFLPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CmCcmRflParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CmCcmRflParameters>);
#endif

DRLIB_END_NAMESPACE

#endif
