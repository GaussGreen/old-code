//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketDataFetcherSRM.hpp
//
//   Description : Helper class for SRM type models to get
//                 data out of market cache
//
//   Author      : Mark A Robson
//
//   Date        : 2 June 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_MARKETDATAFETCHERSRM_HPP
#define EDR_MARKETDATAFETCHERSRM_HPP

#include "edginc/MarketDataFetcher.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

/** Helper class for models to get data out of market cache */        
class MCARLO_DLL MarketDataFetcherSRM : public MarketDataFetcher {
public:
    ~MarketDataFetcherSRM();

    /** default fetch method */
    virtual MarketObjectSP fetch(const MarketData*    market,
                                 const string&        name,
                                 const CClassConstSP& type,
                                 const IModel*        model) const;

    /** Fill fetch stochastic yield curves plus type of smile and model
        specified. */
    MarketDataFetcherSRM(const string&      irCalibSmileType, 
                         const string&      irCalibModelType,
                         bool               getSwaptionVols,
                         bool               useIRVolPair, // FIXME: there is a better approach, check with Guoping
                         const string&      volType,
                         const StringArray& fxVolType,
                         const string&      crCalibSmileType,
                         const string&      cdsVolType);

    /** change the irCalibSmileType */
    void setIrCalibSmileType(const string& irCalibSmileType);

    /** change the irCalibModelType */
    void setIrCalibModelType(const string& irCalibModelType);

    /** Setting this to true, means that a request for a IRVolBase is turned
        into a request for an IRVol otherwise it is turned into one for 
        an IRCalib */
    void setSwaptionVolFlag(bool getSwaptionVols);

private:
    MarketDataFetcherSRM();
    MarketDataFetcherSRM(const MarketDataFetcherSRM& rhs);
    MarketDataFetcherSRM& operator=(const MarketDataFetcherSRM& rhs);
    bool useIRVolPair;
    StringArray fxVolType;
};

typedef smartConstPtr<MarketDataFetcherSRM> MarketDataFetcherSRMConstSP;
typedef smartPtr<MarketDataFetcherSRM> MarketDataFetcherSRMSP;
// derive from MarketDataFetcherSRM to allow easy access to the market
// data retrieved
class MCPathConfigSRM;
	
DRLIB_END_NAMESPACE
#endif
