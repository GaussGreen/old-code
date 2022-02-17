//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Description : Container class for Credit Metrics + RFL 
//                 per name model params 
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CM_RFL_PARAMETERS_HPP
#define QLIB_CM_RFL_PARAMETERS_HPP

#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/CmOnlyParameters.hpp"
#include "edginc/RflOnlyParameters.hpp"

DRLIB_BEGIN_NAMESPACE

/** Contains parameters for Cm  Base Correlation
 * Represents 'model oriented' parameters for a given name. */
class MARKET_DLL CmRflParameters : public RationalisedCreditEngineParameters {

public:
    /** TYPE (for reflection) */        
    static CClassConstSP const TYPE;

    /** Destructor */
    ~CmRflParameters();

    /** Build a CmRflParameters object out of the individual components */
    CmRflParameters(CmOnlyParametersConstSP  ccmParam,
                    RflOnlyParametersConstSP rflParam);

    /** Populates the object with the market data that this object needs */
    virtual void getMarket(const IModel* model, const MarketData* market);
    
    /** Returns the name of this object */
    virtual string getName() const;
    
    /** Returns the cmParameters of this object */
    virtual CmOnlyParametersConstSP getCmParameters() const;   

    /** Returns the cmParameters of this object */
    virtual RflOnlyParametersConstSP getRflParameters() const;   

    /** Returns the engine parameters of the specified type, if any,
        in this CreditEngineParameters object. The parameter specifies what
        type of engine parameters are required. Null is returned if this 
        object is not of the required type */
    virtual CreditEngineParametersConstSP getEngineParams(
        CClassConstSP engineParamsType) const;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    CmRflParameters();
    
    /** Default constructor */
    static IObject* defaultConstructor();

    // ----------------
    // FIELDS
    // ----------------   
    /** Name of this market object */
    string name;
    
    /** Parameters for CM model */
    CmOnlyParametersWrapper cmParameters;

    /** Parameters for RFL model */
    RflOnlyParametersWrapper rflParameters;
};

typedef smartPtr<CmRflParameters> CmRflParametersSP;
typedef smartConstPtr<CmRflParameters> CmRflParametersConstSP;

// Support for wrapper
typedef MarketWrapper<CmRflParameters> CmRflParametersWrapper;

#ifndef QLIB_CMRFLPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CmRflParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CmRflParameters>);
#endif

DRLIB_END_NAMESPACE

#endif
