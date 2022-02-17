//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : MarketDataFetcherCIS.hpp
//
//   Description : Helper class for models that use CDS Basket.
//
//   Author      : Jose Hilera
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#ifndef EDR_MARKETDATAFETCHERCIS_HPP
#define EDR_MARKETDATAFETCHERCIS_HPP

#include "edginc/MarketDataFetcherCDS.hpp"


DRLIB_BEGIN_NAMESPACE

/** Helper class for models to get data out of market cache */        
class MCARLO_DLL MarketDataFetcherCIS: public MarketDataFetcherCDS {
public:
    ~MarketDataFetcherCIS();

    /** Will fetch stochastic yield curves plus type of smile and model
        specified. Will probably need some sort of vol choice once we know
        what we're doing */
    MarketDataFetcherCIS(bool getSwaptionVols);

    /** Will fetch stochastic yield curves plus type of smile and model
        specified. Also accept an extra parameter indicating whether to 
        fetch the index basis or not. */
    MarketDataFetcherCIS(bool getSwaptionVols, 
                         bool calculateIndexBasis,
                         bool useIndexBasis);
 
    /** As above, but allows the index map to be copied in */
    MarketDataFetcherCIS(bool                   getSwaptionVols,
                         ICreditIndexMapConstSP indexMap,
                         const string&          indexMapTreatment,
                         bool                   calculateIndexBasis);
 
    /** Pulls out index basis objects from the cache if required */
    virtual MarketObjectSP fetch(const MarketData*    market,
                                 const string&        name,
                                 const CClassConstSP& type,
                                 const IModel*        model) const;

    /** ensures correct data is retrieved quanto curves. Models should redirect
        calls to getComponentMarketData to this method */
    virtual void getComponentMarketData(const IModel*        model,
                                        const MarketData*    market,
                                        MarketObjectSP       mo) const;

private:
    MarketDataFetcherCIS();
    MarketDataFetcherCIS(const MarketDataFetcherCIS& rhs);
    MarketDataFetcherCIS& operator=(const MarketDataFetcherCIS& rhs);

    CClassConstSP getCDSPSType(CClassConstSP requestedType) const;
    void initialise();

    // Locally apply index basis if applicable
    MarketObjectSP applyIndexBasis(const MarketObjectSP mo) const;

    // Factor out the actual application of index basis
    MarketObjectSP adjustCurve(const MarketObjectSP         mo,
                               const CreditIndexBaseConstSP index) const;

    // Fields
    bool useIndexBasis;                  /* If false, produce a "dummy" index 
                                          * basis. Default: true */
    bool         calculateIndexBasis;    /* If false, fetch the index basis. 
                                          * Default: True*/
    // Copies of information from a CreditMetricsModel (or variety thereof)
    // so that the MDF can make index basis adjustments
    ICreditIndexMapConstSP indexMap;     /* allows single names to become 
                                          * index basis adjusted via its
                                          * corresponding index definition */
    string indexMapTreatment;            /* controls how the index map may 
                                          * be used */
};

typedef smartConstPtr<MarketDataFetcherCIS> MarketDataFetcherCISConstSP;
typedef smartPtr<MarketDataFetcherCIS> MarketDataFetcherCISSP;

DRLIB_END_NAMESPACE
#endif
