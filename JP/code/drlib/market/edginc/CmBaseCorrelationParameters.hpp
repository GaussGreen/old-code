//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Description : Container class for Credit Metrics + Base Correlation 
//                 per name model params 
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CM_BASE_CORRELATION_PARAMETERS_HPP
#define QLIB_CM_BASE_CORRELATION_PARAMETERS_HPP

#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/CmOnlyParameters.hpp"
#include "edginc/BaseCorrelationOnlyParameters.hpp"

DRLIB_BEGIN_NAMESPACE

/** Contains parameters for Cm  Base Correlation
 * Represents 'model oriented' parameters for a given name. */
class MARKET_DLL CmBaseCorrelationParameters : 
    public RationalisedCreditEngineParameters 
{
public:
    /** TYPE (for reflection) */        
    static CClassConstSP const TYPE;

    /** Destructor */
    ~CmBaseCorrelationParameters();

    /** Build a CmBaseCorrelationParameters object based on an old-style
        BaseCorrelationParameters object - which must not be null */
    CmBaseCorrelationParameters(CmOnlyParametersConstSP ccmParam,
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

    CmBaseCorrelationParameters();
    
    /** Default constructor */
    static IObject* defaultConstructor();

    // ----------------
    // FIELDS
    // ----------------   
    /** Name of this market object */
    string name;
    
    /** Parameters for Cm model */
    CmOnlyParametersWrapper cmParameters;

    /** Parameters for BC model */
    BaseCorrelationOnlyParametersWrapper bcParameters;
};

typedef smartPtr<CmBaseCorrelationParameters> CmBaseCorrelationParametersSP;
typedef smartConstPtr<CmBaseCorrelationParameters> CmBaseCorrelationParametersConstSP;

// Support for wrapper
typedef MarketWrapper<CmBaseCorrelationParameters> 
    CmBaseCorrelationParametersWrapper;

#ifndef QLIB_CMBASECORRELATIONPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CmBaseCorrelationParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CmBaseCorrelationParameters>);
#endif

DRLIB_END_NAMESPACE

#endif
