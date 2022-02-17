//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Description : Abstract base class for per name engine specific 
//                 credit parameters.
//
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//      This class and its subtypes are now deprecated - Use types 
//      CmCcmBaseCorrelationParameters or CmCcmRflParameters instead.
//      All non-essential methods in this class have been removed or throw
//      exceptions (note that instances of this class are automatically
//      converted to the new-style parameters in 
//      CreditEngineParameters::convertToNewParamStyle().
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//
//----------------------------------------------------------------------------

#ifndef EDR_COMPOSITECREDITENGINEPARAMETERS_HPP
#define EDR_COMPOSITECREDITENGINEPARAMETERS_HPP

#include "edginc/CreditEngineParameters.hpp"

DRLIB_BEGIN_NAMESPACE

class CompositeCreditEngineParameters;
typedef smartPtr<CompositeCreditEngineParameters> CompositeCreditEngineParametersSP;
typedef smartConstPtr<CompositeCreditEngineParameters> CompositeCreditEngineParametersConstSP;

/** Abstract base class for per name engine specific credit parameters*/
class MARKET_DLL CompositeCreditEngineParameters: public CreditEngineParameters
{
public:
    static CClassConstSP const TYPE;

    virtual ~CompositeCreditEngineParameters();
    
    /** Returns the engine parameters of the specified type, if any,
        in this CompositeCreditEngineParameters object. The parameter specify what
        type of engine parameters are required. Null is returned if the
        relevant type is not found (or there are no embedded parameters in
        this object) */
    virtual CreditEngineParametersConstSP getEngineParams(
        CClassConstSP engineParamsType) const;

protected:
    CompositeCreditEngineParameters(const CClassConstSP& clazz);

private:
    CompositeCreditEngineParameters();

    CompositeCreditEngineParameters(const CompositeCreditEngineParameters& rhs); // don't use
    CompositeCreditEngineParameters& operator=(
        const CompositeCreditEngineParameters& rhs); // don't use
    static void load(CClassSP& clazz);
};


DRLIB_END_NAMESPACE
#endif
