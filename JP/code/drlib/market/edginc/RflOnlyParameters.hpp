//----------------------------------------------------------------------------
//
//   Group       : QR Credit Hybrids
//
//   Description : Container class for RFL-specific per name model params 
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_RFL_ONLY_PARAMETERS_HPP
#define QLIB_RFL_ONLY_PARAMETERS_HPP

#include "edginc/MappingFunction.hpp"
#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(RFLParameters);

class MARKET_DLL RflOnlyParameters : public RationalisedCreditEngineParameters,
                                     public virtual Calibrator::IAdjustable
{

public:
    /** TYPE (for reflection) */        
    static CClassConstSP const TYPE;

    /** Destructor */
    ~RflOnlyParameters();
    
    /** Build an RflOnlyParameters object based on an old-style
        RFLParameters object - which must not be null */
    RflOnlyParameters(RFLParametersConstSP rflParam);

    /** Returns the name of this object */
    virtual string getName() const;
    
   /** Compute the conditional survival probability for a name */
    virtual void condSurvivalProba(double pind,
                                   double tgauss,
                                   double betaHist,
                                   double M,
                                   double& pm) const;

    /** compute the threshold associated to the survival probability 
        for a name */
    virtual void threshold(double pgauss, 
                           double &tgauss, 
                           double betaHist = 0.0) const;

    /** compute the cumulative proba of X=b*M+... */
    virtual double XCumProba(double x, double betaHist = 0.0) const;

    /** compute the density function of X = beta(M)*M+ ... */
    virtual double XDensProba(double x, double betaHist = 0.0) const;

    /** compute the mean proba of b(M)*M */
    virtual double mean(double betaHist = 0.0) const;

    /** compute the var of b(M)*M */
    virtual double variance(double betaHist = 0.0) const;

    // ----------------
    // FIELDS
    // ----------------  
    // IDEALLY build to get getImappingFunction()
    /** beta = f(market factor) */
    IMappingFunctionSP betaCurve;

    /** boolean flag to correct or not the expected loss */
    bool correctMean;

    /** boolean flag to correct or not the variance */
    bool correctVariance;

	/** boolean to say wheteher the beta curve is a tweak from
     *  the original betaHist or ab absolute value */
	bool isTweak;

    /** boolean to say wheteher the CF computations
     *  should be used in the piecewise flat case */
	bool useClosedForm;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    RflOnlyParameters();
    
    /** Default constructor */
    static IObject* defaultConstructor();

    /** Returns true if underlying "betaCurve" is flat or piecewise flat:
     * in that case, there are closed form formulae for mean, variance and
     * integration of the conditional single-name survival probabilities
     * over the market factor.
     * */
    bool isFlatBetaCurve() const;

    // ----------------
    // FIELDS
    // ----------------  
    /** Name of this market object */
    string name;
};

typedef smartPtr<RflOnlyParameters> RflOnlyParametersSP;
typedef smartConstPtr<RflOnlyParameters> RflOnlyParametersConstSP;


// Support for wrapper
typedef MarketWrapper<RflOnlyParameters> RflOnlyParametersWrapper;

#ifndef QLIB_RFLONLYPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<RflOnlyParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<RflOnlyParameters>);
#endif

DRLIB_END_NAMESPACE

#endif //QLIB_RFL_ONLY_PARAMETERS_HPP
