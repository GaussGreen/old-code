//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Asset.hpp
//
//   Description : Asset interface
//
//   Author      : Mark A Robson
//
//   Date        : 15 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_ASSET_HPP
#define EDR_ASSET_HPP
#include "edginc/PriceAsset.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/OutputName.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/AsMultiFactors.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/MarketObservable.hpp"
#include "edginc/Correlation.hpp"

DRLIB_BEGIN_NAMESPACE
class Results;
class AssetNameCollector;

/** Base class for all assets. 
    An asset represents some sort of tradeable item. It has a name which
    can be accessed via
    {@link #GetName()}, a spot price obtainable from
    {@link #GetSpot()},
    an ability to calculate a 
    {@link #SettleDate()}
    given a trade date,
    and an ability to calculate a
    {@link #FwdValue()}, and some sort
    of 
    {@link #CVolProcessed}, accessible via
    {@link #GetProcessedVol()},  which gives access to the asset's volatility.
    <BR><BR>
    Note that the asset does not promise the ability to return a
    {@link CVolBase} directly 
    That would be too strict a requirement for certain
    types of assets such as cross currency baskets. What it does promise is
    that instrument specific data (captured in the 
    {@link CVolRequest}) can be combined with the volatility 
    information within the asset to create an
    {@link CVolProcessed}. The 
    {@link CVolProcessed} interface promises very little
    and which implementation of the {@link CVolProcessed}
    is returned by the {@link GetProcessedVol()}
    method
    depends upon the asset and the {@link CVolRequest}.

*/
class MARKET_DLL CAsset: public MarketObject, 
              public virtual IMarketObservable,
              public virtual IPriceAsset,
              public virtual IAsMultiFactors {
public:
    static CClassConstSP const TYPE;
    virtual ~CAsset();
   
    static const string CCY_TREATMENT_NONE;      // N
    static const string CCY_TREATMENT_VANILLA;   // V (same as N)
    static const string CCY_TREATMENT_STRUCK;    // S
    static const string CCY_TREATMENT_PROTECTED; // P

    /** What an asset needs to do in order for it to be currency struck 
        (it should be derived from Asset but Asset is not an interface and
        C++ wouldn't let us define it here) */
    class MARKET_DLL IStruckable{
    public:
        static CClassConstSP const TYPE;
        virtual ~IStruckable();

        /** Combines market and instrument data together to give a
            Processed Vol. Here the processed volatility is a processed
            struck volatility ie it reflects the combination of this
            asset together with the supplied FX asset and the
            correlation between this CVolBase and the vol of the
            FX. Note that the struckAsset is indeed the struckAsset cf
            'this' which is the non struck asset */
        virtual IVolProcessed* getProcessedVol(
            const CVolRequest* volRequest,
            const CAsset*      struckAsset,
            const FXAsset*     fxAsset,
            const Correlation* eqFXCorr) const = 0;
    };
    /** What additional methods an asset offers if it is currency struck.
        Again, as for IStruckable, it should really be derived from an Asset
        interface */
    class MARKET_DLL IStruck{
    public:
        static CClassConstSP const TYPE;
        virtual ~IStruck();
        /** Returns the fx spot */
        virtual double getFXSpot() const = 0;
        /** Returns the fx forward value */
        virtual double fxFwdValue(const DateTime& date) const = 0;
    };
	
    class MARKET_DLL IQuanto {
    public:
        static CClassConstSP const TYPE;
        virtual ~IQuanto();

        /* Calculate the expected spot price of underlying at given date
         * with no currency protection adjustment */
        virtual double unadjustedFwdValue(const DateTime& date) const = 0;

        /* Calculate the expected spot price of underlying at given dates
         * ie no currency protection adjustment is made */
        virtual void unadjustedFwdValue(const DateTimeArray& dates,
                                        CDoubleArray&        result) const;

        /* Return the Asset/FX correlation */
        virtual const Correlation* getCorrelation() const = 0;

        /* Returns a processed vol - which combines the FX vol market data with
         * the instrument data in the volRequest */
        virtual IVolProcessed* getProcessedFXVol( const CVolRequest* volRequest ) const = 0;
    };

    /** simple class to indicate DDE capability
     **/
    class MARKET_DLL ICanHaveDDE
    {
    public:
        static CClassConstSP const TYPE;
        virtual bool isDDE() const = 0; // return true if want to turn on DDE
        virtual ~ICanHaveDDE();
    };
    
    /////////////////////////////////////////////////////////////////////////
    ////////// Note: Many method now live in IPriceAsset and IGeneralAsset //
    /////////////////////////////////////////////////////////////////////////

    /** returns the asset name */
    virtual string getName() const = 0;

    /** returns the asset name without the effect of any current treatment 
        ie the name that would have been obtained if currency treatment had 
        been CCY_TREATMENT_VANILLA. Default implementation returns getName() */
    virtual string getTrueName() const;

    //// Returns the ccy treatment for the asset. Default returns 
    //// CCY_TREATMENT_NONE
    virtual string getCcyTreatment() const;

    /** Provides default internal name for different ccy treatments */
    static string getSyntheticName(string baseName,
                                   string ccyTreatment,
                                   string ccyCode);

    /** Calculates the expected spot price of the asset at the given date.
        Do not use this repeatedly to calculate values over a set of
        dates (poor performance) - instead use other fwdValue method */
    virtual double fwdValue(const DateTime& date) const = 0;

    /** Quicker version of above when you have multiple dates */
    virtual void fwdValue(const DateTimeArray& dateList,
                          CDoubleArray&        result) const;

    /** Fancy sounding class that can be used to alter the way the
        forward value is calculated (see the relevant fwdValue
        method). Typically, the values will only have meaning to
        certain types of assets. Currently limited to ignoring dividends */
    class MARKET_DLL FwdValueAlgorithm{
    public:
        /** Constructor where you can specify whether dividends (inside an
            equity) are ignored. Note that typically it is sufficient to
            ignore all dividends rather than over an interval since it is
            only the return over the interval that you want */
        FwdValueAlgorithm(bool  ignoreDivs);
        
        /** whether dividends should be ignored */
        bool ignoreDivs() const;
    private:
        bool divsIgnored;
    };

    /** Calculates the expected spot prices of the asset at the given dates
        respecting any 'algorithmic' choices set in the FwdValueAlgorithm */
    virtual void fwdValue(const DateTimeArray&     dateList,
                          const FwdValueAlgorithm& algo,
                          CDoubleArray&            result) const = 0;

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const = 0;

    /** Returns an LN processed vol - which combines the vol market
        data with the instrument data in the volRequest */
    virtual CVolProcessedBS* getProcessedVol(
        const CVolRequestLN* volRequest) const;

    /** If supplied asset wrapper is using the market data cache, then
        retrieves, and makes if necessary, an asset composed of the
        underlying asset together with the requested currency treatment */
    static void getAssetMarketData(const IModel*            model, 
                                   const MarketData*        market,
                                   const string&            ccyTreatment,
                                   const string&            payOutYCName,
                                   MarketWrapper<CAsset>&   asset);

    /** Same as above but takes a yield curve wrapper for convenience */
    static void getAssetMarketData(const IModel*            model, 
                                   const MarketData*        market,
                                   const string&            ccyTreatment,
                                   const YieldCurveWrapper& payOutYC,
                                   MarketWrapper<CAsset>&   asset);

    /** Essentially does assetWrapper.getData(model, market) but switches
        the current 'domestic' yield curve in the model/market data fetcher
        before asking the newly retrieved asset to get its market data.
        This should be used by Assets which contain another asset but
        switch its currency eg ProtAsset, StruckAsset etc */
    static void getAssetInNewCurrency(const IModel*            model, 
                                      const MarketData*        market,
                                      MarketWrapper<CAsset>&   assetWrapper);

    // A simple implementation for single assets
    virtual IMultiFactors* asMultiFactors() const;

    /** record forwards at maturity*/
    virtual void recordFwdAtMat(OutputRequest*  request,
                                Results*        results,
                                const DateTime& maturityDate) const;

     /** Returns the spot value to use when populating a sample during a
        theta shift */
    virtual double getThetaSpotOnDate(const Theta *shift,
                                      const DateTime &date) const;

    /** Returns the sensitive strikes for this asset using the data in
        the vol request */
    virtual void getSensitiveStrikes(
        const CVolRequest*               volRequest,
        OutputNameConstSP                outputName,
        const SensitiveStrikeDescriptor& sensStrikeDesc,
        DoubleArraySP                    sensitiveStrikes) const;

    /** Calculates the expected spot price of the asset at the given date if
        the spot price had the given value spot on spotDate */
    virtual double fwdFwd(const DateTime& spotDate,
                          double          spot, 
                          const DateTime& fwdDate) const;

    /** Array version of fwdFwd */
    virtual void fwdFwd(const DateTime&      spotDate,
                        double               spot, 
                        const DateTimeArray& fwdDates,
                        DoubleArray&         results) const;

    /** Calculates the expected number of shares at given fwdDate,
        based on $1 notional. By default, related convexity adjustment
        is TRUE */
    virtual double expNumberShares( const DateTime& valueDate, 
                                    const DateTime& fwdDate,
                                    const bool& convAdjust) const;

    /** Can this asset physically settle? */
    virtual bool canPhysicallySettle() const;

    /** Retrieve, and make if necessary, an asset composed of the underlying
        asset together with the requested currency treatment. Resulting
        asset IS populated with market data */
    static CAsset* makeAssetUsingCache(const IModel*     model, 
                                       const MarketData* market,
                                       const string&     undAssetName,
                                       const string&     ccyTreatment,
                                       const string&     payOutYCName);

protected:
    CAsset(const CClassConstSP& clazz);
    static void acceptNameCollector(const CAsset*       asset,
                                    AssetNameCollector* collector);

private:
    static void load(CClassSP& clazz);
    CAsset(const CAsset& rhs);
    CAsset& operator=(const CAsset& rhs);
};

typedef CAsset Asset;
// smart pointer support
typedef smartPtr<Asset> AssetSP;
typedef smartConstPtr<Asset> AssetConstSP;
#ifndef QLIB_ASSET_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<Asset>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<Asset>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<Asset>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<Asset>);
#endif
typedef AssetSP CAssetSP;
typedef AssetConstSP CAssetConstSP;

// support for arrays of assets (note array of smart pointers)
typedef array<CAssetSP, CAsset> CAssetArray;
#ifndef QLIB_ASSET_CPP
EXTERN_TEMPLATE(class MARKET_DLL array<CAssetSP _COMMA_ CAsset>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL array<CAssetSP _COMMA_ CAsset>);
#endif
// support for smart pointer to array
typedef smartConstPtr<CAssetArray> CAssetArrayConstSP;
typedef smartPtr<CAssetArray> CAssetArraySP;
#ifndef QLIB_ASSET_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CAssetArray>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CAssetArray>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CAssetArray>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CAssetArray>);
#endif

// support for wrapper class
typedef MarketWrapper<CAsset> CAssetWrapper;
#ifndef QLIB_ASSET_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CAsset>);
EXTERN_TEMPLATE(IObjectSP MARKET_DLL FieldGetInLine<CAssetWrapper>(CAssetWrapper* t));
EXTERN_TEMPLATE(void MARKET_DLL FieldSetInLine<CAssetWrapper>(CAssetWrapper* t,
                                                   IObjectSP o));
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CAsset>);
INSTANTIATE_TEMPLATE(IObjectSP MARKET_DLL FieldGetInLine<CAssetWrapper>(CAssetWrapper* t));
INSTANTIATE_TEMPLATE(void MARKET_DLL FieldSetInLine<CAssetWrapper>(CAssetWrapper* t,
                                                        IObjectSP o));
#endif

/** specialisations of arrayObjectCast */
template <> class MARKET_DLL arrayObjectCast<CAssetWrapper>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const CAssetWrapper& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(CAssetWrapper& value);

    /** Turns the IObjectSP into a DateTime */
    static CAssetWrapper fromIObject(IObjectSP& value);
};
// arrays of wrappers (note array of structures)
typedef array<CAssetWrapper, CAssetWrapper> CAssetWrapperArray;
#ifndef QLIB_ASSET_CPP
EXTERN_TEMPLATE(class MARKET_DLL array<CAssetWrapper _COMMA_ CAssetWrapper>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL array<CAssetWrapper _COMMA_ CAssetWrapper>);
#endif

typedef smartPtr<CAssetWrapperArray> CAssetWrapperArraySP;
#ifndef QLIB_ASSET_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CAssetWrapperArray>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CAssetWrapperArray>);
#endif
DRLIB_END_NAMESPACE
#endif
