//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketDataFetcherCDS.hpp
//
//   Description : Helper class for models that use CDS Par Spreads.
//                 It's in the mcarlo directory since we might need this for
//                 a credit mc
//
//   Author      : Mark A Robson
//
//   Date        : 29 Nov 2004
//
//----------------------------------------------------------------------------

#ifndef EDR_MARKETDATAFETCHERCDS_HPP
#define EDR_MARKETDATAFETCHERCDS_HPP

#include "edginc/MarketDataFetcher.hpp"
#include "edginc/CreditIndexMap.hpp"


DRLIB_BEGIN_NAMESPACE

/** Helper class for models to get data out of market cache.
    This class doesn't need to exist any more since it just configures the
    MarketDataFetcher instance and doesn't override any methods.
    Just need to replace constructors with 'static' constructors */        
class MCARLO_DLL MarketDataFetcherCDS: public MarketDataFetcher {
public:
    ~MarketDataFetcherCDS();

    /** Fill fetch stochastic yield curves plus type of smile and model
        specified. Will probably need some sort of vol choice once we know
        what we're doing */
    MarketDataFetcherCDS(bool getSwaptionVols);

	/** Fill fetch stochastic yield curves plus type of smile and model
        specified. Will probably need some sort of vol choice once we know
        what we're doing */
    MarketDataFetcherCDS(bool getSwaptionVols,
			 bool getCreditVols,
			 CClassConstSP creditVolType = NULL);

protected:
    MarketDataFetcherCDS();
private:
    MarketDataFetcherCDS(const MarketDataFetcherCDS& rhs);
    MarketDataFetcherCDS& operator=(const MarketDataFetcherCDS& rhs);
    void initialise(bool getSwaptionVols);
};

typedef smartConstPtr<MarketDataFetcherCDS> MarketDataFetcherCDSConstSP;
typedef smartPtr<MarketDataFetcherCDS> MarketDataFetcherCDSSP;

DRLIB_END_NAMESPACE
#endif
