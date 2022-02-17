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

#include "edginc/config.hpp"
#include "edginc/CDOQuotesMaturityFilter.hpp"
#include "edginc/TimePoint2D.hpp"

DRLIB_BEGIN_NAMESPACE

/** Destructor */
CDOQuotesMaturityFilter::~CDOQuotesMaturityFilter() {}

/** Called immediately after object constructed */
void CDOQuotesMaturityFilter::validatePop2Object() {
    try {
    // TODO
    } catch (exception& e){
        throw ModelException(e, "DiscreteCDOQuotesWeights::validatePop2Object");
    }
}

/** Get the market data needed */
void CDOQuotesMaturityFilter::getMarket(const IModel* model, const MarketData* market)
{
    try
    {
        //Nothing to do
    }
    catch(exception& e)
    {
        throw ModelException(e, "CDOQuotesMaturityFilter::getMarket");
    }
}

/** Get the name of the market object */
string CDOQuotesMaturityFilter::getName() const
{
    return name;
}

/** Invoked when Class is 'loaded' */
void CDOQuotesMaturityFilter::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CDOQuotesMaturityFilter, clazz);
    SUPERCLASS(MarketObject);
    IMPLEMENTS(ICDOQuotesWeights);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(maturity, "Maturity Date for which quotes are assigned a non 0 weight");
}

/** Public constructor  */
CDOQuotesMaturityFilter::CDOQuotesMaturityFilter(
            DateTime maturity) :
    MarketObject(TYPE),
    maturity(maturity)
{}

/** Main method : returns a weight for the tranche with a give maturity and given strikes */
double CDOQuotesMaturityFilter::getWeight(DateTime maturity, double K1, double K2) const
{
    const TimePoint2D pt = TimePoint2D(maturity, K1, K2);
    if (maturity == this->maturity)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}


/** Private constructor (only build instances of that class using reflection) */
CDOQuotesMaturityFilter::CDOQuotesMaturityFilter() :
    MarketObject(TYPE),
    maturity()
{}

/** Default constructor */
IObject* CDOQuotesMaturityFilter::defaultConstructor() {
    return new CDOQuotesMaturityFilter();
}

/** TYPE for CDOQuotesMaturityFilter */
CClassConstSP const CDOQuotesMaturityFilter::TYPE = CClass::registerClassLoadMethod(
    "CDOQuotesMaturityFilter", typeid(CDOQuotesMaturityFilter), CDOQuotesMaturityFilter::load);
    
/* external symbol to allow class to be forced to be linked in */
bool CDOQuotesMaturityFilterLoad(){
    return (CDOQuotesMaturityFilter::TYPE != 0);
}

DECLARE(CDOQuotesMaturityFilter);

DRLIB_END_NAMESPACE
