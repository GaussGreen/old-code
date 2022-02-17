//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOQuotesMaturityFilter.hpp
//
//   Description : CDOQuotesMaturityFilter is a quotes filter for the calibration.
//                  It assigns a weight equal to 1 for all quotes with the
//                  maturity provided and 0 otherwise.
//
//   Author      : Sebastien Gay    
//
//   Date        : October 2006
//
//----------------------------------------------------------------------------

#ifndef CDO_QUOTES_MATURITY_FILTER_HPP
#define CDO_QUOTES_MATURITY_FILTER_HPP

#include "edginc/ICDOQuotesWeights.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * CDOQuotesMaturityFilter is an implementation of ICDOQuotesWeights
 * that puts a zero weight to all quotes with a maturity different from the
 * input maturity date and a weight equal to one for the quotes with
 * that maturity.
 * */
class PRODUCTS_DLL CDOQuotesMaturityFilter :
            public MarketObject,
            public virtual ICDOQuotesWeights  {
public:
	
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Destructor */
    virtual ~CDOQuotesMaturityFilter();

    /** Public constructor  */
    CDOQuotesMaturityFilter(DateTime maturity);

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** Main method : returns a weight for the tranche with a give maturity and given strikes */
    virtual double getWeight(DateTime maturity, double K1, double K2) const;

    /** Get the market data needed */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Get the name of the market object */
    virtual string getName() const;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Only build instances of that class using reflection */
    CDOQuotesMaturityFilter();

    /** Default constructor */
    static IObject* defaultConstructor();

    /// FIELDS
    /** string name of the market object */
    string name;

    /** maturity filtered */
    DateTime maturity;
};

DECLARE(CDOQuotesMaturityFilter);

DRLIB_END_NAMESPACE

#endif /* CDO_QUOTES_MATURITY_FILTER_HPP*/

