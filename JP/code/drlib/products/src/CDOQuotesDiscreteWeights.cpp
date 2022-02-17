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

#include "edginc/config.hpp"
#include "edginc/CDOQuotesDiscreteWeights.hpp"

DRLIB_BEGIN_NAMESPACE

/** Destructor */
CDOQuotesDiscreteWeights::~CDOQuotesDiscreteWeights() {}

/** Called immediately after object constructed */
void CDOQuotesDiscreteWeights::validatePop2Object() {
    try {
    	if ( lowStrikes->size()  != weights->size() ||
             highStrikes->size() != weights->size() ||
             expiries->size()    != weights->size() )
        {
            throw ModelException("CDOQuotesDiscreteWeights::validatePop2Object",
                        "input arrays have different size");
        }

        // builds the map of (points, weights)
        int n = (*weights).size();
        for (int i=0; i<n; ++i)
        {
            pointWeightMap[ TimePoint2D( ((*expiries)[i])->toDate(valueDate),
                                        (*lowStrikes)[i],
                                        (*highStrikes)[i] ) ]
                            = (*weights)[i] ;
        }
    } catch (exception& e){
        throw ModelException(e, "CDOQuotesDiscreteWeights::validatePop2Object");
    }
}

/** Get the market data needed */
void CDOQuotesDiscreteWeights::getMarket(const IModel* model, const MarketData* market)
{
    try
    {
        market->GetReferenceDate(valueDate);
    }
    catch(exception& e)
    {
        throw ModelException(e, "CDOQuotesDiscreteWeights::getMarket");
    }
}

/** Get the name of the market object */
string CDOQuotesDiscreteWeights::getName() const
{
    return name;
}

/** Access to value date */
const DateTime& CDOQuotesDiscreteWeights::getValueDate() const {
    return valueDate;
}

/** Invoked when Class is 'loaded' */
void CDOQuotesDiscreteWeights::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CDOQuotesDiscreteWeights, clazz);
    SUPERCLASS(MarketObject);
    IMPLEMENTS(ICDOQuotesWeights);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(defaultWeight, "default weight value for point not provided explicitly");
    FIELD(lowStrikes, "array of lower strikes");
    FIELD(highStrikes, "array of higher strikes");
    FIELD(expiries, "array of expiries");
    FIELD(weights, "array of weights corresponding to each point");
}

/** Public constructor  */
CDOQuotesDiscreteWeights::CDOQuotesDiscreteWeights(
            double defaultWeight,
            DoubleArraySP      lowStrikes,
            DoubleArraySP      highStrikes,
            ExpiryArraySP      expiries,
            DoubleArraySP      weights) :
    MarketObject(TYPE),
    defaultWeight(defaultWeight),
    lowStrikes(lowStrikes),
    highStrikes(highStrikes),
    expiries(expiries),
    weights(weights)
{}

/** Main method : returns a weight for the tranche with a give maturity and given strikes */
double CDOQuotesDiscreteWeights::getWeight(DateTime maturity, double K1, double K2) const
{
    const TimePoint2D pt = TimePoint2D(maturity, K1, K2);
    map<const TimePoint2D, double, TimePoint2D::CompareDateC1C2>::const_iterator it = pointWeightMap.find(pt);
    if (it != pointWeightMap.end())
    {
        return (*it).second;
    }
    else
    {
        return defaultWeight;
    }
}


/** Private constructor (only build instances of that class using reflection) */
CDOQuotesDiscreteWeights::CDOQuotesDiscreteWeights() :
    MarketObject(TYPE),
    defaultWeight(0),
    lowStrikes(0),
    highStrikes(0),
    expiries(0),
    weights(0)
{}

/** Default constructor */
IObject* CDOQuotesDiscreteWeights::defaultConstructor() {
    return new CDOQuotesDiscreteWeights();
}

/** TYPE for CDOQuotesDiscreteWeights */
CClassConstSP const CDOQuotesDiscreteWeights::TYPE = CClass::registerClassLoadMethod(
    "CDOQuotesDiscreteWeights", typeid(CDOQuotesDiscreteWeights), CDOQuotesDiscreteWeights::load);
    
/* external symbol to allow class to be forced to be linked in */
bool CDOQuotesDiscreteWeightsLoad(){
    return (CDOQuotesDiscreteWeights::TYPE != 0);
}

DECLARE(CDOQuotesDiscreteWeights);

DRLIB_END_NAMESPACE
