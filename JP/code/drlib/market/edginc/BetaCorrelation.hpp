
#ifndef QLIB_BETA_CORRELATION_HPP
#define QLIB_BETA_CORRELATION_HPP

#include "edginc/MarketObject.hpp"

DRLIB_BEGIN_NAMESPACE
class BetaCorrelation;
typedef smartPtr<BetaCorrelation> BetaCorrelationSP;
typedef smartConstPtr<BetaCorrelation> BetaCorrelationConstSP;
typedef array<BetaCorrelationSP, BetaCorrelation> BetaCorrelationArray;
#ifndef QLIB_BETA_CORRELATION_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<BetaCorrelation>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<BetaCorrelation>);
EXTERN_TEMPLATE(class MARKET_DLL array<BetaCorrelationSP _COMMA_ BetaCorrelation>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<BetaCorrelation>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<BetaCorrelation>);
INSTANTIATE_TEMPLATE(class MARKET_DLL array<BetaCorrelationSP _COMMA_ BetaCorrelation>);
#endif

/** This class implements a simple beta correlation for multi-asset diffusion
    the meaning of beta is asymmetric -- it is always a statement that
    Tier-2 (otherwise uncorrelated) asset is beta-correlated to one or more
    of more of Tier-1 (strongly correlated) assets. */
class MARKET_DLL BetaCorrelation: public MarketObject {
public:
    static CClassConstSP const TYPE;

    /** Validation */
    void validatePop2Object();

	virtual void getMarket(const IModel* model, const MarketData* market);

    /** returns the correlation as a simple double */
    double getBetaCorrelation() const { return betaCorr; }
    string getTier1Name() const { return assetTier1; }
    string getTier2Name() const { return assetTier2; }

    /** returns the correlation's name */
    virtual string getName() const { return name; }

    /** Initialises this piece of market data - records the pair of names
        idenitfying the correlation */
    virtual void initialise(MarketData* market);

    /** Searches the supplied array of correlations looking for a
        correlation that is between the names of the two assets
        supplied */
    static double findBeta(
        const BetaCorrelationArray& correlations,
        const string&           assetTier1,
        const string&           assetTier2);


    /** Searches the supplied array of correlations looking for a
    correlation that is between the names of the two assets
    supplied.  Returns correlations.end() if the index is not found.
    TO DO:  The following does a linear search.  This isn't very
    efficient for Sampras where we have correlations between 10,000 assets.*/
    static BetaCorrelationArray::const_iterator findIndex(
        const BetaCorrelationArray& correlations,
        const string&           assetTier1,
        const string&           assetTier2);

    /** constructor */
    BetaCorrelation(const string& name,        // optional
                const string& nameAssetTier1,
                const string& nameAssetTier2,
                double        betaCorr);

private:
    static void load(CClassSP& clazz);
    static IObject* defaultBetaCorrelation();
    BetaCorrelation();

    string          name;   // optional
    string          assetTier1;
    string          assetTier2;
    double          betaCorr;
};

DRLIB_END_NAMESPACE
#endif // BETA_CORRELATION_HPP
