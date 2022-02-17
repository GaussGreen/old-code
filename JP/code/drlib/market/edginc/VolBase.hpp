//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolBase.hpp
//
//   Description : Abstract vol interface
//
//   Author      : Mark A Robson
//
//   Date        : 13 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef VOLBASE_HPP
#define VOLBASE_HPP

#include "edginc/DateTime.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/GetMarket.hpp"

DRLIB_BEGIN_NAMESPACE
class IVolProcessed;
// for backwardsc compatibility
typedef IVolProcessed  CVolProcessed;
class CVolRequest;
class CAsset;
class FXAsset;
class Correlation;

/** Base class for volalities.
    Apart from being able to
    identify itself through its
    {@link  #getName()} method, it can return an 
    {@link  CVolProcessed}. An 
    {@link  CVolProcessed} represents the
    combination of the instrument and the market data.
    <BR><BR>
    Note that this interface is generally hidden behind the 
    {@link CAsset} and {@link  CVolProcessed}
    interfaces. At a product level, this would be rarely used.
    <BR><BR>
    But what is the {@link CAsset} doing in the interface?
    <BR>
    The idea is to allow the {@link  CVolBase} access to any asset
    related information that it needs. Currently, this is only the spot price
    (for forward starting interpolation) but it is not hard to imagine a 
    situation where forward prices, for example, would be needed.
    The alternatives appear to be hard coding in the interface the precise
    asset data needed (ie spot price) or storing the asset information
    in the {@link VolInterp} object. The latter
    approach has the advantage of simplifying this interface but would mean
    that the asset related information would always be calculated regardless of
    whether the volatility object actually needed it. It would also cause 
    major problems should the data needed from the asset ever depend on the
    volatility.
*/
class MARKET_DLL CVolBase: public MarketObject,
                virtual public IGetMarket {
public:
    static CClassConstSP const TYPE;
    friend class CVolBaseHelper;

    /** Combines market and instrument data together to give a
        Processed Vol. The 'implementation' here routes through 
        VolProcessedDispatch which will fail if no suitable method has been
        registered */
    virtual IVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const CAsset*      asset) const = 0;

    /** Combines market and instrument data together to give a
        Processed Vol. Here the processed volatility is a processed
        struck volatility ie it reflects the combination of this
        CVolBase together with the supplied FX asset and the
        correlation between this CVolBase and the vol of the
        FX. */
    virtual IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      eqAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const = 0;

    
    /** populate from market cache - default implementation provided */
    virtual void getMarket(const IModel* model, const MarketData* market);

    virtual ~CVolBase();

    /** Override the one in MarketObject to allow for automatic conversion
        of vols to fx vols */
    virtual void convert(IObjectSP&    object,
                         CClassConstSP requiredType) const;
protected:
    CVolBase(const CClassConstSP& clazz);
private:
    CVolBase(const CVolBase &rhs);
    CVolBase& operator=(const CVolBase& rhs);
};

typedef smartConstPtr<CVolBase> CVolBaseConstSP;
typedef smartPtr<CVolBase> CVolBaseSP;
#ifndef QLIB_VOLBASE_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CVolBase>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CVolBase>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CVolBase>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CVolBase>);
#endif

// support for wrapper class
typedef MarketWrapper<CVolBase> CVolBaseWrapper;
#ifndef QLIB_VOLBASE_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CVolBase>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CVolBase>);
#endif

DRLIB_END_NAMESPACE

#endif
