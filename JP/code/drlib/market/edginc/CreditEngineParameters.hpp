//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Description : Abstract base class for per name engine specific 
//                 credit parameters.
//                 The main method in this class is 
//                 getEngineParams(engineParamsType), which returns the 
//                 parameters of the required type if available within this
//                 object. Note there is currently no support for models to
//                 fetch "pre-processed" parameters similar to equity's
//                 volatilities, since this is not anticipated as a requirement
//                 for credit.
//
//   Date        : 28 Oct 2004
//
//----------------------------------------------------------------------------

#ifndef EDR_CREDITENGINEPARAMETERS_HPP
#define EDR_CREDITENGINEPARAMETERS_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE


class CreditEngineParameters;
typedef smartPtr<CreditEngineParameters> CreditEngineParametersSP;
typedef smartConstPtr<CreditEngineParameters> CreditEngineParametersConstSP;

// support for wrapper class
typedef MarketWrapper<CreditEngineParameters> CreditEngineParametersWrapper;
typedef smartPtr<CreditEngineParametersWrapper> CreditEngineParametersWrapperSP;


/** Abstract base class for per name engine specific credit parameters*/
class MARKET_DLL CreditEngineParameters: public MarketObject {

public:
    static CClassConstSP const TYPE;

    virtual ~CreditEngineParameters();
    
    /** Returns the engine parameters of the specified type, if any,
        in this CreditEngineParameters object. The parameter specify what
        type of engine parameters are required. An exception is thrown if 
        the relevant type is not found  */
    virtual CreditEngineParametersConstSP getEngineParams(
        CClassConstSP engineParamsType) const;

    CreditEngineParametersConstSP getEngineParams() const 
          { return CreditEngineParametersConstSP::attachToRef(this); }

    /** If required, converts the engineParams from an old-style parameter
        into a new-style parameter. If the engineParams are old-style 
        parameters it will assume that beta and decretionBeta are NOT null;
        they can be null if using new-style parameters
        NOTE: Beta sensitivities used to be reported under the name of the 
        portfolio (which is the asset's name, really) .
        Therefore, for backwards compatibility reasons, in case of old-style
        model parameters need to name the CmOnlyParameters after the
        PortfolioName, so that those greeks appear under the right name */
    static CreditEngineParametersWrapper convertToNewParamStyle(
        CreditEngineParametersWrapper engineParams,
        const string nameOfPortfolio,
        CDoubleSP beta, 
        CDoubleSP decretionBeta);

    /* Return badness of parameters, ie, some non-negative measure of the 
       undesirability of the parameters, eg, a measure of "smoothness". 
       This function is intended for use in calibration. */
    virtual double ParameterBadness() const { return 0.0; }

protected:
    CreditEngineParameters(const CClassConstSP& clazz);

private:
    CreditEngineParameters(const CreditEngineParameters& rhs); // don't use
    CreditEngineParameters& operator=(
        const CreditEngineParameters& rhs); // don't use
    static void load(CClassSP& clazz);
};

// support for arrays
typedef array<CreditEngineParametersSP, 
              CreditEngineParameters> CreditEngineParametersArray;

DRLIB_END_NAMESPACE

#endif
