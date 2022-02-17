//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Description : Container class for Credit Metrics + CMM + Base Correlation 
//                 per name model params 
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CM_CMM_BASE_CORRELATION_PARAMETERS_HPP
#define QLIB_CM_CMM_BASE_CORRELATION_PARAMETERS_HPP

#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/CmOnlyParameters.hpp"
#include "edginc/CcmOnlyParameters.hpp"
#include "edginc/BaseCorrelationOnlyParameters.hpp"

DRLIB_BEGIN_NAMESPACE

/** Contains parameters for Credit Metrics + Base Correlation.
 * Represents 'model oriented' parameters for a given name. */
class MARKET_DLL CmCcmBaseCorrelationParameters : 
    public RationalisedCreditEngineParameters 
{
public:
    /** TYPE (for reflection) */        
    static CClassConstSP const TYPE;

    /** Destructor */
    ~CmCcmBaseCorrelationParameters();

    /** Build a CmCcmBaseCorrelationParameters object out of the 
        individual components */
    CmCcmBaseCorrelationParameters(CmOnlyParametersConstSP cmParam,
                                   CcmOnlyParametersConstSP ccmParam,
                                   BaseCorrelationOnlyParametersConstSP bcParam);

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

    CmCcmBaseCorrelationParameters();
    
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

    /** Parameters for BC model */
    BaseCorrelationOnlyParametersWrapper bcParameters;
};

typedef smartPtr<CmCcmBaseCorrelationParameters> CmCcmBaseCorrelationParametersSP;
typedef smartConstPtr<CmCcmBaseCorrelationParameters> CmCcmBaseCorrelationParametersConstSP;

// Support for wrapper
typedef MarketWrapper<CmCcmBaseCorrelationParameters> 
    CmCcmBaseCorrelationParametersWrapper;

#ifndef QLIB_CMCCMBASECORRELATIONPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CmCcmBaseCorrelationParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CmCcmBaseCorrelationParameters>);
#endif

DRLIB_END_NAMESPACE

#endif
