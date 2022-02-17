//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ICDOQuotesWeights.hpp
//
//   Description : ICDOQuotesWeights is an interface for the definition of weights
//                  of each quote in the calibration objective function
//                  (see trancheIndexLeastSquareFit.[c|h]pp)
//
//   Author      : Sebastien Gay    
//
//   Date        : October 2006
//
//----------------------------------------------------------------------------

#ifndef I_CDO_QUOTES_WEIGHTS_HPP
#define I_CDO_QUOTES_WEIGHTS_HPP

#include "edginc/CDOQuotes.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * QuoteWeights is an interface that is used to compute weights
 * of each tranche quotes in the objective fct of the index calibration
 * */
class PRODUCTS_DLL ICDOQuotesWeights : public virtual IObject {
public:
	
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Destructor */
    virtual ~ICDOQuotesWeights();

    /** Main method : returns a weight for the tranche with a give maturity and given strikes */
    virtual double getWeight(DateTime maturity, double K1, double K2) const = 0;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

DECLARE(ICDOQuotesWeights);

typedef MarketWrapper<ICDOQuotesWeights> ICDOQuotesWeightsWrapper;

DRLIB_END_NAMESPACE

#endif /* I_CDO_QUOTES_WEIGHTS_HPP*/
