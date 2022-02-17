//----------------------------------------------------------------------------
//
//   Group       : QR Credit Hybrids
//
//   Description : Container class for RFL-specific per name model params 
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_RFL_MIXTURE_PARAMETERS_HPP
#define QLIB_RFL_MIXTURE_PARAMETERS_HPP

#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include <set>

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE(RflMixtureParameters);

// This implementation supports only piecewise constant factor loading functions
// [Actually we are implementing a bit more. Instead of representing the factor loading function
//  b(M) we look at the whole systemic term: f(M) = b(M)M which is always (even in the multi-factor case
//  a real-valued function. With just one Gaussian factor there is an efficient closed-form for the 
//  unconditional default probability of f is piecewise linear. Note that this is more general than 
//  having piecewise constant b (which corresponds to all line elements in f having intercept zero).]
class MARKET_DLL RflMixtureParameters : public RationalisedCreditEngineParameters,
                                     public virtual Calibrator::IAdjustable
{
public:
    /** TYPE (for reflection) */        
    static CClassConstSP const TYPE;

    /** Build an RflMixtureParameters object based on an old-style
        RFLParameters object - which must not be null 
    RflMixtureParameters(RFLParametersConstSP rflParam);*/

    /** Returns the name of this object */
    virtual string getName() const;
    
    /** Compute the conditional survival probability for a name */
    void condSurvivalProba(double barrier,
                           double factorValue,
                           double& pm) const;

    // Return systemic term for given factor value
    double f(double factorValue) const;  

    // Return unconditional survival probability for Gaussian-mixture distributed factor given
    //  the mixture parameters and the default barrier level. 
    // Note that the mixture has two components one of which has weight w, mean c and variance 1; the other 
    //  has weight 1-w, mean 0 and variance 1.
    // We supply this functionality because this important special case allows a closed form computation
    double survivalProb(double barrier,   // Value of default barrier
                        double c = 0.0,   // displacement of mixture components (usually < 0)
                        double w = 0.0) const;  // weight of displaced component (usually w << 1)
                        
    /** compute the default barrier associated to the survival probability for a name */
    void defaultBarrier(double survProb, 
                        double &barrier, 
                        double c = 0.0,
                        double w = 0.0) const;

    // Write thresholds to provided container
    void Thresholds(set<double>& thresholds) const;

    void fieldsUpdated(const CFieldArray& fields);

    double ParameterBadness() const;

    void validatePop2Object();

    // Set global beta value, ie, set all beta entries to this value
    void SetGlobalBeta(double b) { for (size_t i = 0; i < m_betas.size(); ++i) m_betas[i] = b; }

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    RflMixtureParameters();
    
    /** Default constructor */
    static IObject* defaultConstructor();

// data
    // These three data members represent a piecewise linear f which may be discontinuous
    //  at the thresholds. Between thresholds f is given by a slope and an intercept.
    // Thresholds (note that only m_thresholds[0] is actually a threshold; the rest are additive increments)
    CDoubleArray m_thresholds;
    // Slopes (Note that m_betas.size() == m_thresholds.size()+1)
    CDoubleArray m_betas;
    // Intercepts
    CDoubleArray m_alphas;

    // Badness weight 
    // This really should be in the objective function, but that does not seem practical...
    double m_badnessWeight;

    /** Name of this market object */
    string name;
};

// Support for wrapper
typedef MarketWrapper<RflMixtureParameters> RflMixtureParametersWrapper;

DRLIB_END_NAMESPACE

#endif //QLIB_RFL_MIXTURE_PARAMETERS_HPP
