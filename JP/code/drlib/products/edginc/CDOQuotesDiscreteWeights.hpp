//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOQuotesDiscreteWeights.hpp
//
//   Description : CDOQuotesDiscreteWeights is an implementation of
//                  ICDOQuotesWeights that specifies a default weight
//                  and a list of weights for a list of points defined
//                  as (maturity, K1, K2)
//
//   Author      : Sebastien Gay    
//
//   Date        : October 2006
//
//----------------------------------------------------------------------------

#ifndef CDO_QUOTES_DISCRETE_WEIGHTS_HPP
#define CDO_QUOTES_DISCRETE_WEIGHTS_HPP

#include "edginc/ICDOQuotesWeights.hpp"
#include "edginc/TimePoint2D.hpp"
#include <map>

DRLIB_BEGIN_NAMESPACE

/**
 * DiscreteCDOQuoteWeights is an implementation of ICDOQuotesWeights
 * that allow to specify weights on a point per point basis 
 * with a general weight for the points for which the weight
 * is not provided
 * */
class PRODUCTS_DLL CDOQuotesDiscreteWeights :
            public MarketObject,
            public virtual ICDOQuotesWeights  {
public:
	
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Destructor */
    virtual ~CDOQuotesDiscreteWeights();

    /** Public constructor  */
    CDOQuotesDiscreteWeights(
            double defaultWeight,
            DoubleArraySP      lowStrikes,
            DoubleArraySP      highStrikes,
            ExpiryArraySP      expiries,
            DoubleArraySP      weights);

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** Main method : returns a weight for the tranche with a give maturity and given strikes */
    virtual double getWeight(DateTime maturity, double K1, double K2) const;

    /** Access to value date */
    const DateTime& getValueDate() const;

    /** Get the market data needed */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Get the name of the market object */
    virtual string getName() const;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Only build instances of that class using reflection */
    CDOQuotesDiscreteWeights();

    /** Default constructor */
    static IObject* defaultConstructor();

    /// FIELDS

    /** string name of the market object */
    string name;

    /** array containing the list of lower stikes */
    DoubleArraySP lowStrikes;

    /** array containing the list of upper stikes */
    DoubleArraySP highStrikes;

    /** array containing the list of expiries */
    ExpiryArraySP   expiries;

    /** array containing the list of weight */
    DoubleArraySP   weights;

    /** default weight value for points not provided */
    double defaultWeight;

    /** Value date [Transient] */
    DateTime valueDate;

    /// TRANSIENT FIELDS built in the constructor
    map<const TimePoint2D, double, TimePoint2D::CompareDateC1C2> pointWeightMap;
};

DECLARE(CDOQuotesDiscreteWeights);

DRLIB_END_NAMESPACE

#endif /* CDO_QUOTES_DISCRETE_WEIGHTS_HPP*/

