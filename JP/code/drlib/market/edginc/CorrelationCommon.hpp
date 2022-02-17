
#ifndef EDR_CORRELATION_COMMON_HPP
#define EDR_CORRELATION_COMMON_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/CorrelationBase.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class CorrelationCommon;
DECLARE(CorrelationCommon);

class MARKET_DLL CorrelationCommon: public CorrelationBase
{
public:
    static CClassConstSP const TYPE;

    static const string BENCHMARK_EXPIRY;
    static const int BENCHMARK_TIME;

    /** Validation */
    void validatePop2Object() = 0;

	//virtual void getMarket(const IModel* model, const MarketData* market);

    /** returns the correlation's name */
    virtual string getName() const;

    /** Initialises this piece of market data - records the pair of names
        idenitfying the correlation */
    virtual void initialise(MarketData* market);

    /** Searches the supplied array of correlations looking for a
        correlation that is between the names of the two assets
        supplied */
    static CorrelationCommonSP find(const CorrelationCommonArray& correlations,
                       const string&           asset1,
                       const string&           asset2);

    /** Like find(), but user can specify how to handle "not found" situation (true => throw exception; false => return 0.0 */
    static CorrelationCommonSP lookupCorrelation(const CorrelationCommonArray& corrs,
                                    const string& name1,
                                    const string& name2,
                                    bool strictCorr);

    /** Searches the supplied array of correlations looking for a
    correlation that is between the names of the two assets
    supplied.  Returns correlations.end() if the index is not found.
    TO DO:  The following does a linear search.  This isn't very 
    efficient for Sampras where we have correlations between 10,000 assets.*/
    static CorrelationCommonArray::const_iterator findIndex(
        const CorrelationCommonArray& correlations,
        const string&           asset1,
        const string&           asset2);

    static int findCorrArrayIndex(const CorrelationCommonArray& correlations,
                                  const string&                 asset1,
                                  const string&                 asset2);

    /** constructor */
    CorrelationCommon(CClassConstSP clazz,
                      const string& name,
                      const string& nameAsset1,
                      const string& nameAsset2);

    const string& getAsset1Name() const { return asset1; }
    const string& getAsset2Name() const { return asset2; }

protected:
    string          name;   // optional
    string          asset1;
    string          asset2;

protected:
    CorrelationCommon(CClassConstSP clazz);

private:
    static void load(CClassSP& clazz);
    CorrelationCommon();
};


DRLIB_END_NAMESPACE
#endif // CORRELATION_COMMON_HPP
