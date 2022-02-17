//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Equity.hpp
//
//   Description : Equity - forward price and spot
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------

#ifndef EQUITY_HPP
#define EQUITY_HPP

#include <string>
#include "edginc/Settlement.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/MuParallel.hpp"
#include "edginc/MuSpecial.hpp"
#include "edginc/MuPointwise.hpp"
#include "edginc/DividendList.hpp"
#include "edginc/BorrowCurve.hpp"
#include "edginc/Spot.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/SpotLevel.hpp"
#include "edginc/SpotShift.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/WrapperNameCollector.hpp"
#include "edginc/Theta.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/BorrowParallelShift.hpp"
#include "edginc/BorrowLevel.hpp"
#include "edginc/CreditCurve.hpp"
#include "edginc/DeltaSurface.hpp"
#include "edginc/CtsDivOverride.hpp"
#include "edginc/DeltaDDE.hpp"
#include "edginc/AssetHistoryContainer.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

class EquityBase;
class CriticalDateCollector;
/** Represents the simplest form of a stock. No volatility data nor any 
    support for currency treatment */
class MARKET_DLL Equity: public CObject, 
                         virtual public INamedObject,
                         virtual public IRestorableWithRespectTo<Spot>,
                         virtual public DeltaDDE::RestorableShift,
                         virtual public MuParallel::IShift,
                         virtual public MuSpecial::IShift,
                         virtual public MuPointwise::IShift,
                         virtual public SpotLevel::Shift,
                         virtual public SpotShift::Shift,
                         virtual public RhoBorrowParallel::RestorableShift,
                         virtual public RhoBorrowPointwise::IRestorableShift,
                         virtual public BorrowParallelShift::Shift,
                         virtual public BorrowLevel::Shift,
                         virtual public Theta::IShift,
                         virtual public DeltaSurface::RestorableShift,
                         virtual public CtsDivOverride::IShift {
public:
    static CClassConstSP const TYPE;
    friend class EquityHelper;

    enum {
        INCLUDE_ALL_DIVS = 0,
        IGNORE_ALL_DIVS,
        IGNORE_DOLLAR_DIVS
    } TDivTreat; // $unregistered

    /** Pull out the yield curve and value date from market data 
    and asset history */
    void getMarket(const IModel* model, const MarketData* market);

    /** Constructor */
    Equity(const string&       name, 
           const DateTime&     valueDate,
           const DateTime&     stockDate,
           double              spot, 
           const Settlement*   settle, 
           const YieldCurve*   yc,
           const DividendList* divList,
           const BorrowCurve*  borrowCurve);

    /** Pop2Object validation called from constructor and after deserialisation */
    virtual void validatePop2Object();

    /** Returns the stock name */
    string getName() const;

    /** Returns the spot price */
    virtual double spot() const;

    /** Returns fair value of stock price */
    double fairValue() const;

    /** Returns valueDate */
    DateTime getValueDate() const; 
    /** Returns the name (not the ISO code) of the yield curve used to
        grow the stock */
    string getYCName() const;

    /** Returns the ISO code (not the name) of the yield curve used to
        grow the stock */
    string getYCIsoCode() const;

    YieldCurveWrapper  getYC() const;

    /** returns smart pointer to market holiday object */
    HolidayConstSP getMarketHolidays() const;

    // the IMarketObservable interface for retrieving a single sample
    virtual double pastValue(const DateTime&             sampleDate,
                             const ObservationType*      obsType,
                             const ObservationSource*    source,
                             const FixingType*           fixType,
                             const IObservationOverride* overrides,
                             const SamplingConvention*   sampleRule) const;

    // IMarketObservable - retrieve a single observation date
    // Returns false if obs is to be omitted
    virtual bool observationDate(const DateTime&           sampleDate,
                                 const ObservationSource*  source,
                                 const SamplingConvention* sampleRule,
                                 DateTime*                 obsDate) const;

    // the IMarketObservable interface for retrieving past samples events
    virtual double addPastSampleEvent(const DateTime&             sampleDate,
                                    const ObservationType*      obsType,
                                    const ObservationSource*    source,
                                    const FixingType*           fixType,
                                    const IObservationOverride* overrides,
                                    const SamplingConvention*   sampleRule,
                                    PastSamplesCollector*        collector) const;

    // the IMarketObservable interface for 
    // is the given date a holiday for the relevant source
    virtual bool isHoliday(const DateTime& sampleDate,
                           const ObservationSource*   source) const;

    /** Calculated the forward price on the given data. Another one is
        required which takes an array of dates - to do */
    virtual double fwdValue(const DateTime& date) const;

    /** Calculates an array of forward prices which are assumed to be in 
        ascending order */
    virtual void fwdValue(const DateTimeArray& dateList,
                          CDoubleArray&        result) const;

    /** Calculates an array of forward prices on the supplied dates
        (which are assumed to be in ascending order) ignoring all dividends */
    void fwdValue(const DateTimeArray&             dateList,
                  const CAsset::FwdValueAlgorithm& algo,
                  CDoubleArray&                    result) const;

    /** Calculates an array of forward prices on the supplied dates
        (which are assumed to be in ascending order) as if the stock had
        the given value on the supplied date */
    void fwdFwd(const DateTime&      valueDate,
                const DateTime&      stockDate,
                double               stockPrice,
                const DateTimeArray& dateList,
                CDoubleArray&        result) const;

    /** Calculates the stock floor at given dates assuming no dollar dividends 
        are paid passed endDate. In the absence of borrows and dividend yields 
        (continous and discrete) the stock floor is the PV of future dividends
        till endDate. */
    void calcStockFloor(const DateTimeArray& dateList,
                        const DateTime&      endDate,
                        CDoubleArray&        result) const;

    /** Calculates spot price of the pseudo asset Z 
        that is associated with a dollar dividend paying equity S. See below for detail. */
    double getPseudoSpot(const DateTime& endDate) const;

    /** Calculates an array of forward prices which are assumed to be in 
        ascending order. The forward prices are those of the pseudo asset Z 
        that is associated with a dollar dividend paying equity S. See below for detail. */
    void calcPseudoFwdValue(const DateTimeArray& dateList,
                            const DateTime&      endDate,
                            CDoubleArray&        result) const;

    /** Given a vol base and vol request, returns a processed vol where the processed vol
        is that of the pseudo asset Z that is associated with a dollar dividend paying equity S. 
        The pseudo asset has the property that Z = S - L where L is the floor of the stock 
        price; it grows like S, except for the dollar dividend drops. 
        NB Only simple equities and started options are supported at present. */
    CVolProcessedBS* getPseudoProcessedVol(const CVolBase*        vol,
                                           const CVolRequestLN*   volRequest,
                                           const EquityBase*      asset,
                                           const DateTime&        endDate,
                                           const DateTimeArray&   divCritDates,
                                           bool                   isCall,
                                           int                    noExerciseWindow) const;

    /** Returns the settlement date for the given trade date */
    DateTime settles(const DateTime &tradeDate) const;

    /** Returns name identifying this object for Delta */
    virtual string sensName(const Spot*) const;
    /** Shifts the object using given shift (see Delta::Shift)*/
    virtual TweakOutcome sensShift(const PropertyTweak<Spot>& shift);
    /** Restores the object to its original form */
    void sensRestore(const PropertyTweak<Spot>& shift);

    /** Returns name identifying this object for DeltaDDE */
    virtual string sensName(DeltaDDE* shift) const;
    /** Shifts the object using given shift (see Delta::Shift)*/
    virtual bool sensShift(DeltaDDE* shift);
    /** Restores the object to its original form */
    virtual void sensRestore(DeltaDDE* shift);

    /** Returns name identifying this object for MU_PARALLEL */
    virtual string sensName(MuParallel* shift) const;
    /** Shifts the object using given shift. This is a wrapper for the
        DividendList MU_PARALLEL shift method */
    virtual bool sensShift(MuParallel* shift);

    /** Returns name identifying this object for MU_S */
    virtual string sensName(MuSpecial* shift) const;
    /** Shifts the object using given shift. This is a wrapper for the
        DividendList MU_S shift method */
    virtual bool sensShift(MuSpecial* shift);

    /** Returns name identifying this object for MU_POINTWISE */
    virtual string sensName(MuPointwise* shift) const;
    /** Shifts the object using given shift. This is a wrapper for the
        DividendList MU_POINTWISE shift method */
    virtual bool sensShift(MuPointwise* shift);
  
    /** Implements SpotLevel scenario */
    /** Returns name identifying this object for SpotLevel */
    virtual string sensName(SpotLevel* shift) const;
    /** Shifts the object using given shift (see SpotLevel::Shift)*/
    virtual bool sensShift(SpotLevel* shift);

    /** Returns name identifying this object for Delta */
    virtual string sensName(SpotShift* shift) const;
    /** Shifts the object using given shift (see Delta::Shift)*/
    virtual bool sensShift(SpotShift* shift);

    /** Returns the name of the borrow curve - used to determine
        whether to tweak the object */
    virtual string sensName(RhoBorrowParallel* shift)const;

    /** Shifts the object using given shift */    
    virtual bool sensShift(RhoBorrowParallel* shift);

    /** Restores the object to its original form */
    virtual void sensRestore(RhoBorrowParallel* shift);

    /** Returns the name of the borrow curve - used to determine whether 
        to tweak the object */
    virtual string sensName(RhoBorrowPointwise* shift)const;

    /** Return the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this  borrow curve */
    virtual ExpiryArrayConstSP sensExpiries(RhoBorrowPointwise* shift)const;
    
    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
    virtual bool sensShift(RhoBorrowPointwise* shift);

    /** Restores the object to its original form */
    virtual void sensRestore(RhoBorrowPointwise* shift);  

 
    /** Shifts the object using given shift (see Theta::Shift)*/
    virtual bool sensShift(Theta* shift);

    /** Returns the name of the borrow curve - used to determine
        whether to tweak the object */
    virtual string sensName(BorrowParallelShift* shift)const;

    /** Shifts the object using given shift */    
    virtual bool sensShift(BorrowParallelShift* shift);

    /** Returns the name of the borrow curve - used to determine
        whether to tweak the object */
    virtual string sensName(BorrowLevel* shift)const;

    /** Shifts the object using given shift */    
    virtual bool sensShift(BorrowLevel* shift);

    /** Returns name identifying this object for DeltaSurface */
    virtual string sensName(DeltaSurface* shift) const;
    /** Shifts the object using given shift (see DeltaSurface::Shift)*/
    virtual bool sensShift(DeltaSurface* shift);
    /** Restores the object to its original form */
    virtual void sensRestore(DeltaSurface* shift);

    // The pieces for CtsDivOverride
    virtual string sensName(CtsDivOverride* shift) const;
    virtual bool sensShift(CtsDivOverride* shift);

    /** Returns the dividends in the underlying equity */
    DividendListConstSP  getDivList() const;

    bool getWillZeroNegFwd() const;
    void setWillZeroNegFwd(bool flag);

    /** Given an equity, return an equity where the dividends have been converted according the 
        following rule. 
        If not a dollar div, leave as is.
        If dollar div and ex date <= 3Y, leave as is
        If dollar div and 3Y < ex date < 5Y, convert dollar div into 1 dollar div component
        and 1 div yield component where the dollar div component is equal to 
        proportion * initial dollar div amount
        and proportion varies linearly from 1 to 0 between 3Y and 5Y
        If dollar div and ex date >= 5Y, convert dollar div into yield (no dollar div component).
        NB The conversion into yield is done in such way that fwds at ex-div dates are unchanged. */
    static Equity* convertDollarDivs(const Equity& equity);

    // calculates the implied dividend yield between two dates (ignoring borrow costs)
    double dividendYield(const DateTime& startDate, const DateTime& endDate);

    const DateTime& getDivTransPeriodStartDate() const;   
    static DateTime calcDivTransPeriodStartDate(const DateTime& valueDate);

    const DateTime& getDivTransPeriodEndDate() const;   
    static DateTime calcDivTransPeriodEndDate(const DateTime& valueDate);

    void makeRiskyEquity(ICreditCurveSP creditSpreads,
                         const  DateTime *maturityDate=0);

    /** Replaces the div list with a new one */
    void setDivList(const DividendList* newDivList);

    /** get bc to allow theta shift for fwd caching in FastQuote */
    BorrowCurveSP getBorrow() const {return borrowCurve;};

    const Settlement& getSettlement() {return *settlement;}

protected:
    Equity();
    Equity(const Equity &rhs);
    Equity& operator=(const Equity& rhs);

    // allow child class
    Equity(const CClassConstSP& clazz) : CObject(clazz) {};

    static void acceptValueDateCollector(const Equity* equity, 
                                         CValueDateCollector* collector);

    static void acceptCriticalDateCollector(const Equity* equity, 
                                            CriticalDateCollector* collector);

    static void acceptWrapperNameCollector(const Equity* asset, 
                                           WrapperNameCollector* collector);

    void fwdFwd(const DateTime&      valueDate,
                const DateTime&      stockDate,
                double               stockPrice,
                int                  divTreat,
                const DateTimeArray& dateList,
                CDoubleArray&        result) const;

    void rollOverYieldDivs(const DateTime& rollDate, const Theta* shift);

    void validateValueDate();

    /** implies spot price at value date given the quoted price at the
        quote date */
    double implySpot(const DateTime& quoteDate, 
                     const double&   quotePrice ) const;      

    int pseudoShortVars(const DateTimeArray&     shortBenchDates,
                        const DateTime&          valueDate,
                        double                   spotPrice,
                        double                   strike,
                        const DateTime&          endDate,
                        const CVolProcessedBSSP& volProcessedBS,
                        CDoubleArray&            pseudoVars) const;

    /* Needed for optimization */
    struct MARKET_DLL FwdDiffFunc{
        const Equity& equity;
        Dividend&     dividend;
        const         DateTime& exDivDate;
        double        tgtExDivFwd;
        double        initAmount;

        FwdDiffFunc(const Equity& equity,
                    Dividend&     dividend,
                    double        tgtExDivFwd);

        double operator()(double amount) const;
    };
    // fields
    string             name;
    DateTime           valueDate;
    DateTime           stockDate;
    double             stockPrice;
    SettlementSP       settlement;
    YieldCurveWrapper  yc;
    DividendListSP     divList;
    BorrowCurveSP      borrowCurve;
    AssetHistoryContainerSP history;

    static const string divTransPeriodStart;
    DateTime            divTransPeriodStartDate;
    static const string divTransPeriodEnd;
    DateTime            divTransPeriodEndDate;

    bool               hasBorrow;  // indicates non-zero borrowCurve

    // willZeroNegFwd flag: if true, then don't throw
    // a negative forward, return 0 instead.  If false,
    // usual behavior: throw if forward computed is negative.
    bool willZeroNegFwd;
};

typedef smartConstPtr<Equity> EquityConstSP;
typedef smartPtr<Equity> EquitySP;
#ifndef QLIB_EQUITY_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<Equity>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<Equity>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<Equity>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<Equity>);
#endif

DRLIB_END_NAMESPACE

#endif
