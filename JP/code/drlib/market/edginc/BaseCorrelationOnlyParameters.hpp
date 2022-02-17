//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Description : Container class for Base Correlation specific per name
//                 model params 
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_BASE_CORRELATION_ONLY_PARAMETERS_HPP
#define QLIB_BASE_CORRELATION_ONLY_PARAMETERS_HPP

#include "edginc/IndexWeights.hpp"
#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE(BaseCorrelationParameters);


/** Contains parameters for Base Correlation.
 * Represents 'model oriented' parameters for a given name. */
class MARKET_DLL BaseCorrelationOnlyParameters : 
    public RationalisedCreditEngineParameters 
{
public:
    /** TYPE (for reflection) */        
    static CClassConstSP const TYPE;

    /** Destructor */
    ~BaseCorrelationOnlyParameters();

    /** public constructor */
    BaseCorrelationOnlyParameters(string name, IndexWeightsSP indexWeights);

    /** Build a BaseCorrelationOnlyParameters object based on an old-style
        BaseCorrelationParameters object - which must not be null */
    BaseCorrelationOnlyParameters(BaseCorrelationParametersConstSP bcParam);

    /** Populates the object with the market data that this object needs */
    virtual void getMarket(const IModel* model, const MarketData* market);
    
    /** Returns the name of this object */
    virtual string getName() const;
    
    /** Returns the index weights associated to this name */
    IndexWeightsConstSP getIndexWeights() const;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    BaseCorrelationOnlyParameters();
    
    /** Default constructor */
    static IObject* defaultConstructor();

    // ----------------
    // FIELDS
    // ----------------
    /** Name of this market object */
    string name;
    
    /** The index weights associated to this name */
    IndexWeightsSP indexWeights;
};

DECLARE(BaseCorrelationOnlyParameters);

// Support for wrapper
typedef MarketWrapper<BaseCorrelationOnlyParameters> 
    BaseCorrelationOnlyParametersWrapper;

#ifndef QLIB_BASECORRELATIONONLYPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<BaseCorrelationOnlyParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<BaseCorrelationOnlyParameters>);
#endif

DRLIB_END_NAMESPACE

#endif
