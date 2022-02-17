//
//   Author      : Mark A Robson
//
//   Date        : 7th August 2006
//
//----------------------------------------------------------------------------

#ifndef MDF_UTIL_HPP
#define MDF_UTIL_HPP


DRLIB_BEGIN_NAMESPACE
class MarketDataFetcher;

/** As the name suggests a slightly weak class that has some static methods
    for altering MarketDataFetchers. Ideally the class shouldn't exist and 
    the methods should live on MarketDataFetcher itself. However, that class
    is in riskmgr directory and the methods need to manipulate market data
    types. Ideally we move MarketDataFetcher to the market directory ... 
    but that is part of a bigger refactorisation ... */
class MARKET_DLL MDFUtil {
public:
    /** Sets whether correlation skew shall be used or not. */
    static void setCorrSkewMode(MarketDataFetcher& mdf,
                                bool               setUseCorrSkew);

    /** Retrieves whether correlation skew is being used */
    static bool useCorrSkew(const MarketDataFetcher& mdf);
    
    /** Sets whether correlation term structure shall be used or not. */
    static void setCorrTermMode(MarketDataFetcher& mdf,
                                bool               setUseCorrTemStructure);
    
    /** Retrieves whether correlation term structure is being used */
    static bool useCorrTerm(const MarketDataFetcher& mdf);

    /** Specific models can request the support of currency basis, if
        it is supplied */
    static void setUseCurrencyBasis(MarketDataFetcher& mdf,
                                    bool               flag);

    /** Sets whether VarSwapBasis objects are retrieved or not */
    static void setVarSwapBasis(MarketDataFetcher& mdf,
                                bool               flag);

    /** Configure MDF so that only IR Swaption vols are retrieved ie 
        redirects IRVolBase to IRVol and does not retrieve any IRCalib
        type parameters */
    static void setUseSimpleIRVol(MarketDataFetcher& mdf);
private:
    MDFUtil();
    ~MDFUtil();
    MDFUtil(const MDFUtil& rhs);
    MDFUtil& operator=(const MDFUtil& rhs);
};

DRLIB_END_NAMESPACE
#endif
