//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Description : BaseCorrelationParameters contains 
//                 parameters for Base Correlation and optional CCMParameters. 
//                 Represents 'model oriented' parameters for a given name.
//
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//      This class is now deprecated: Use Cm(Ccm)BaseCorrelationParameters 
//      instead.
//      All non-essential methods in this class have been removed or throw
//      exceptions (note that instances of this class are automatically
//      converted to the new-style parameters in 
//      CreditEngineParameters::convertToNewParamStyle().
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//
//   Author      : Antoine Gregoire
//
//   Date        : May 2005
//
//----------------------------------------------------------------------------

#ifndef EDR_BASE_CORRELATION_PARAMETERS_HPP
#define EDR_BASE_CORRELATION_PARAMETERS_HPP

#include "edginc/CompositeCreditEngineParameters.hpp"
#include "edginc/CCMParameters.hpp"
#include "edginc/IndexWeights.hpp"

DRLIB_BEGIN_NAMESPACE

/** BaseCorrelationParameters contains parameters for Base Correlation and
 * optional CCMParamters.
 * Represents 'model oriented' parameters for a given name. */
class MARKET_DLL BaseCorrelationParameters : public CompositeCreditEngineParameters {

public:
    /** TYPE (for reflection) */        
    static CClassConstSP const TYPE;

    /** Destructor */
    ~BaseCorrelationParameters();

    /** Populates the object with the market data that this object needs */
    virtual void getMarket(const IModel* model, const MarketData* market);
    
    /** Returns the name of this object */
    virtual string getName() const;

    /** Returns the index weights of this object */
    IndexWeightsConstSP getIndexWeights() const;
        
    /** Returns the CCMParameters of this object */
    CCMParametersConstSP getCCMEngineParams() const;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Only build instances of that class using reflection */
    BaseCorrelationParameters();
    
    /** Default constructor */
    static IObject* defaultConstructor();

    // ---------------
    // OPTIONAL FIELDS
    // ---------------
    
    /** Parameters for CCM model */
    CCMParametersWrapper ccmParameters;

    // ----------------
    // MANDATORY FIELDS
    // ----------------
    
    /** Name of this market object */
    string name;
    
    /** The index weights associated to this name */
    IndexWeightsSP indexWeights;
};

typedef smartPtr<BaseCorrelationParameters> BaseCorrelationParametersSP;
typedef smartConstPtr<BaseCorrelationParameters> BaseCorrelationParametersConstSP;

DRLIB_END_NAMESPACE

#endif
