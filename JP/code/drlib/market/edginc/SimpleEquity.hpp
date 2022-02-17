//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SimpleEquity.hpp
//
//   Description : Single factor, no currency treatment, equity asset
//
//   Author      : Mark A Robson
//
//   Date        : 15 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef SIMPLE_EQUITY_HPP
#define SIMPLE_EQUITY_HPP
#include "edginc/EquityBase.hpp"
#include "edginc/CanBeRisky.hpp"
#include "edginc/CreditCurve.hpp"
#include "edginc/PhysicalDelivery.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implementation of CAsset for a single stock with no currency treatment */
class MARKET_DLL SimpleEquity: public EquityBase,
                    public virtual ICanBeRisky,
                    public virtual Asset::IStruckable,
                    public virtual ICanPhysicallySettle {
public:
    static CClassConstSP const TYPE;
    friend class SimpleEquityHelper;
    friend class PseudoSimpleEquity;
    friend class IrConverter;

    virtual ~SimpleEquity();

    /** returns the asset name */
    virtual string getName() const;

    SimpleEquity(const Equity*   equity,
                 const CVolBase* vol,
                 CClassConstSP const type = TYPE);

    /** Allows asset to be currency struck.
        Combines market and instrument data together to give a
        Processed Vol. Here the processed volatility is a processed
        struck volatility ie it reflects the combination of this
        asset together with the supplied FX asset and the
        correlation between this CVolBase and the vol of the
        FX. Note that the struckAsset is indeed the struckAsset cf
        'this' which is the non struck asset */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      struckAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;

    /** Populates the object with the market data that this object
        needs.  In particular this will just call the equivalen function
        in EquityBase and then retrieve the AssetHistory(s)*/
    void getMarket(const IModel* model, const MarketData* market);

    /** adds the credit spread to the asset's growth curve */
    virtual void makeRisky(ICreditCurveSP creditSpreads,
        const  DateTime *maturityDate=0);

    // for given trade date, return settle date and name of the asset
    virtual void delivery(const DateTime&               tradeDate, 
                          double                        quantity, 
                          double                        price,
                          PhysicalDeliveryByAssetArray* byAsset) const;
private:
    SimpleEquity();
    SimpleEquity(const SimpleEquity& rhs);
    SimpleEquity& operator=(const SimpleEquity& rhs);

    string name;
};

typedef smartPtr<SimpleEquity> SimpleEquitySP;
typedef smartConstPtr<SimpleEquity> SimpleEquityConstSP;

// support for wrapper class
typedef MarketWrapper<SimpleEquity> SimpleEquityWrapper;

DRLIB_END_NAMESPACE
#endif
