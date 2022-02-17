//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DividendCollector.hpp
//
//   Description : Collects dividends from assets
//
//   Author      : Mark A Robson
//
//   Date        : 25 Oct 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_DIVIDENDCOLLECTOR_HPP
#include "edginc/Collector.hpp"
#include "edginc/DividendList.hpp"
#include "edginc/Asset.hpp"
#include "edginc/OutputRequestUtil.hpp"

#define EDR_DIVIDENDCOLLECTOR_HPP
DRLIB_BEGIN_NAMESPACE
/** Collects dividends from objects, typically assets */
class MARKET_DLL DividendCollector: public CObject,
                         public virtual ICollector {
public:
    enum TDividendConversion {
        NONE = 0,
        DOLLAR_TO_YIELD,
        YIELD_TO_DOLLAR
    };
    
    static CClassConstSP const TYPE;

    /** Creates a dividend collector which collects divs between
        (startDate, endDate] from the supplied asset. Any yield divs
        are converted into dollar divs whilst any continuous divs will
        cause a failure. Struck dividends <= valDate are converted
        using fx spot. Yield dividends (in the interval) < valDate
        will cause a failure. Use getDividends() to retrieve the dividends
        or hasYieldDividends() to see if there were any yield
        dividends. Handles all types of assets */
    static DividendListSP divsBetweenDates(const Asset*    asset,
                                           const DateTime& valDate,
                                           const DateTime& startDate,
                                           const DateTime& endDate,
                                           const TDividendConversion conversion = YIELD_TO_DOLLAR);

    /** Returns the 'underlying' divs between (startDate, endDate] from
        the supplied asset in the currency of the 'underlying'. Any yield
        divs are left as is, whilst any continuous divs will cause a
        failure. Handles limited types of assets. Only provided for backwards
        compatibility - please use divsBetweenDates instead */
    static DividendListSP undDivsBetweenDates(const Asset*    asset,
                                              const DateTime& startDate,
                                              const DateTime& endDate,
                                              const TDividendConversion conversion = YIELD_TO_DOLLAR);

    /** Creates a dividend collector which will collect divs whose
        ex-div date is between (startDate, endDate] from asset
        supplied in the call to go(). Any yield divs are converted
        into dollar divs whilst any continuous divs will cause a
        failure. The optional payDateAdj can be used to modify
        dividend pay dates. Struck dividends with pay date <= valDate
        are converted using fx spot. Yield dividends (in the interval)
        < valDate will cause a failure. Use getDividends() to retrieve
        the dividends or hasYieldDividends() to see if there were any
        yield dividends. Handles all types of assets */
    DividendCollector(DividendList::IDivAdjuster*        divAdj, // optional
                      const DateTime&                    valDate,
                      const DateTime&                    startDate,
                      const DateTime&                    endDate,
                      OutputRequestUtil::KnownCashFlows* knownCFs = 0,
                      const TDividendConversion conversion = YIELD_TO_DOLLAR);

    /** Same as constructor above but additionally calls go(topAsset) */
    DividendCollector(const Asset*                    topAsset,
                      DividendList::IDivAdjuster*     divAdj, // optional
                      const DateTime&                 valDate,
                      const DateTime&                 startDate,
                      const DateTime&                 endDate,
                      const TDividendConversion conversion = YIELD_TO_DOLLAR);

    /** Collects dividends (using information specified in constructor) from
        the supplied asset) */
    void go(const Asset*    asset);

    /** Collects fx/ccy (using information specified in constructor) from
        the supplied asset) but use supplied divs*/
    void go(const Asset*    asset, DividendListSP divs);

    /** Returns the dividends that were found in the asset */
    DividendListSP getDividends() const;

    /** Indicates whether any yield dividends were found in the interval
        specified in the constructor */
    bool hasYieldDividends() const;

    /** Indicates whether any dollar dividends were found in the interval
        specified in the constructor */
    bool hasDollarDividends() const;

    /** Add the divs from the supplied dividends  between the start/end dates 
        to the collector. If fx asset is non null, it is used to scale
        dollar divs (by the forward rate of the asset) - used for struck
        assets */
    void addDivs(const DividendListConstSP& divsToAdd,
                 const FXAsset*             fxAsset); // optional

    /** Add the divs from the supplied asset between the start/end dates 
        to the collector. The dividends are scaled by the forward rate of 
        the [fx] asset (if fxAsset not null) */
    void addDivs(const Asset*               undAsset, // not struck
                 const FXAsset*             fxAsset);


    /* takes the divs from the asset and add them to the visitor,
       scale by weight first. Also turns yield divs into dollar ones
       and validates against continuous divs */
    void processComponentAsset(const Asset*   asset,
                               double         weight);

    /** returns the ISO code for the currency of the 'current' dividend.
        This method can be used in the adjust() method to look at the ISO
        code for each dividend in turn */
    const string& divCcyISOCode() const;

private:
    const DividendCollector* getCurrentCollector() const;

    DividendCollector(const DateTime& valDate,
                      const DateTime& startDate,
                      const DateTime& endDate,
                      const TDividendConversion conversion = YIELD_TO_DOLLAR);
    static void load(CClassSP& clazz);
    /// fields ///
    int                 conversion; // $unregistered

    DateTime            valDate;   // no scaling of divs in the past $unregistered
    DateTime            start;     /* get divs due to be paid > startDate, $unregistered */
    DateTime            end;       /* and <= endDate. $unregistered */
    double              weight;    /* scale divs added by this $unregistered */
    mutable DividendListSP divList;   /* holds collected divs $unregistered */
    const Asset*        asset;     /* presence denotes need to convert  // $unregistered
                                      yield divs to dollar divs -a reference */
    bool                hasYieldDivs;  // true => yield divs found $unregistered
    bool                hasDollarDivs; // true => dollar divs found $unregistered
    bool                undDivsMode; // don't scale struck $unregistered
    DividendList::IDivAdjuster* divAdj; /* adjust divs $unregistered */
    DividendCollector*  currentCol; /* either null or the one handling the  // $unregistered
                                       current div (used for call backs) */
    string              ccyISOCode; /* ISO for dividends $unregistered */
    /* access via getCurrentCollector. Optional - collects known cash flows */
    OutputRequestUtil::KnownCashFlows*  knownCFs;  // $unregistered
    bool                isOverrideDivs; // don't pull divs from underlying just the ccy/fx details $unregistered
    DividendListSP       overrideDivs;    // override divs $unregistered
};

DRLIB_END_NAMESPACE
#endif
