//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Description : Container class for Credit Metrics + CCM 
//                 per name model params 
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CM_CCM_PARAMETERS_HPP
#define QLIB_CM_CCM_PARAMETERS_HPP

#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/CmOnlyParameters.hpp"
#include "edginc/CcmOnlyParameters.hpp"

DRLIB_BEGIN_NAMESPACE

/** Contains parameters for Credit Metrics + CCM 
 * Represents 'model oriented' parameters for a given name. */
class MARKET_DLL CmCcmParameters : public RationalisedCreditEngineParameters {

public:
    /** TYPE (for reflection) */        
    static CClassConstSP const TYPE;

    /** Destructor */
    ~CmCcmParameters();

    /** Build a CmCcmParameters object based on an old-style
        CcmParameters object - which must not be null */
    CmCcmParameters(CmOnlyParametersConstSP  cmParam,
                    CcmOnlyParametersConstSP ccmParam);

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

    CmCcmParameters();
    
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
};

typedef smartPtr<CmCcmParameters> CmCcmParametersSP;
typedef smartConstPtr<CmCcmParameters> CmCcmParametersConstSP;

// Support for wrapper
typedef MarketWrapper<CmCcmParameters> 
    CmCcmParametersWrapper;

#ifndef QLIB_CMCCMPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CmCcmParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CmCcmParameters>);
#endif

DRLIB_END_NAMESPACE

#endif
