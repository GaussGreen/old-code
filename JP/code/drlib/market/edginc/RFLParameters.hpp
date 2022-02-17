//----------------------------------------------------------------------------
//
//   Group       : QR Credit Hybrids
//
//   Filename    : RFLParameters.hpp
//
//   Description : RFLParameters contains 
//                 parameters for RFL and optional CCMParameters.
//
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//      This class is now deprecated: Use Cm(Ccm)RflParameters instead
//      All non-essential methods in this class have been removed or throw
//      exceptions (note that instances of this class are automatically
//      converted to the new-style parameters in 
//      CreditEngineParameters::convertToNewParamStyle().
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//
//----------------------------------------------------------------------------

#ifndef QLIB_RFL_PARAMETERS_HPP
#define QLIB_RFL_PARAMETERS_HPP

#include "edginc/CCMParameters.hpp"
#include "edginc/CompositeCreditEngineParameters.hpp"
#include "edginc/MappingFunction.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL RFLParameters : public CompositeCreditEngineParameters,
                                 public virtual Calibrator::IAdjustable
{
public:
    /** TYPE (for reflection) */        
    static CClassConstSP const TYPE;

    /** Destructor */
    ~RFLParameters();

    /** Populates the object with the market data that this object needs */
    virtual void getMarket(const IModel* model, const MarketData* market);
    
    /** Returns the name of this object */
    virtual string getName() const;

    /** Returns the CCMParameters of this object */
    CCMParametersConstSP getCCMEngineParams() const;
    
    // IDEALLY build to get getImappingFunction()
    /** beta = f(market factor) */
    IMappingFunctionSP betaCurve;

    /** boolean flag to correct or not the expected loss */
    bool correctMean;

    /** boolean flag to correct or not the variance */
    bool correctVariance;

	/** boolean to say wheteher the beta curev is a tweak from
     *  the original betaHist or ab absolute value
	 */
	bool isTweak;

    /** boolean to say wheteher the CF computations
     *  should be used in the piecewise flat case
	 */
	bool useClosedForm;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Only build instances of that class using reflection */
    RFLParameters();
    
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
};

typedef smartPtr<RFLParameters> RFLParametersSP;
typedef smartConstPtr<RFLParameters> RFLParametersConstSP;

DRLIB_END_NAMESPACE

#endif //QLIB_RFL_PARAMETERS_HPP
