//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PseoSimpleEquity.hpp
//
//   Description : pseudo equity - see below.
//
//   Date        : 30 Nov 2001
//
//
//----------------------------------------------------------------------------

#ifndef PSEUDO_SIMPLE_EQUITY_HPP
#define PSEUDO_SIMPLE_EQUITY_HPP

#include "edginc/SimpleEquity.hpp"
#include "edginc/EquityCache.hpp"

DRLIB_BEGIN_NAMESPACE
/** Implementation of the "pseudo" asset Z that is associated with a 
    dollar dividend paying SimpleEquity S. 
    The pseudo asset Z grows like S except that (i) it doesn't drop by the 
    dollar dividend amounts and (ii) it starts at a different level, today's 
    pseudo spot price, which is equal to S - L. L is the stock floor (of which 
    more below). 
    With any pseudo asset is associated an horizon time, or end date, beyond 
    which it is assumed that no dollar dividends are paid. In the absence of 
    borrows and dividend yields (continous and discrete) the stock floor L is 
    the PV of all future dollar dividends until the end date. At the end date, 
    we always have Z = S (i.e., L = 0 at end date). This means that if S pays 
    no dollar dividend beyond the end date then Z and S are undistinguishable 
    passed that date. 
    NB Most of the methods declared in this class are just wrappers around
    methods implemented in Equity. It may make more sense to move the 
    implementation of those methods here. (This needs to be investigated. The 
    currency protection case will probably tell us more about what to do about
    this issue.) */
class MARKET_DLL PseudoSimpleEquity: public CAsset {
public:
    static CClassConstSP const TYPE;
    friend class PseudoSimpleEquityHelper;

    /** returns the spot price */
    virtual double getSpot() const;

    /** Returns fair value of stock price */
    virtual double fairValue() const;

    /** returns the asset name */
    virtual string getName() const;

    /** returns the asset name */
    virtual string getTrueName() const;

    /** Returns the name (not the ISO code) of the yield curve used to
        grow the stock */
    string getYCName() const;

   /** Calculate the settlement date associated with a given trade date */
    virtual DateTime settleDate(const DateTime& tradeDate) const;

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

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    /** Given a vol base and a vol request, returns a processed vol where the 
        processed vol is that of the pseudo asset Z that is associated with a 
        dollar dividend paying equity S. The pseudo asset has the property that 
        Z = S - L where L is the floor of the stock price; it grows like S, 
        except for the dollar dividend drops. 
        NB Only simple equities and started options are supported at present. */
    virtual CVolProcessed * getProcessedVol(
        const CVolRequest* volRequest) const;
    
    /** Calculates the expected spot price of the asset at the given date */
    virtual double fwdValue(const DateTime& date) const;

    /** Calculates the expected spot price of the asset at each of the
        given dates */
    /** Calculates an array of forward prices which are assumed to be in 
        ascending order. The forward prices are those of the pseudo asset Z 
        that is associated with a dollar dividend paying equity S. See below for detail. */
    virtual void fwdValue(const DateTimeArray& dateList,
                          CDoubleArray&        result) const;

    /** Calculates the expected spot prices of the asset at the given dates
        respecting any 'algorithmic' choices set in the FwdValueAlgorithm */
    virtual void fwdValue(const DateTimeArray&     dateList,
                          const FwdValueAlgorithm& algo,
                          CDoubleArray&            result) const;

    /** Calculates the expected spot price of the asset at the given date if
        the spot price had the given value spot on spotDate */
    virtual double fwdFwd(const DateTime& spotDate,
                          double          spot, 
                          const DateTime& fwdDate) const;

    /** Calculates the stock floor at given dates assuming no dollar dividends 
        are paid passed horizonDate. In the absence of borrows and dividend yields 
        (continous and discrete) the stock floor is the PV of future dividends
        till horizonDate. */
    void calcStockFloor(const DateTimeArray& dateList,
                        CDoubleArray&        result) const;

    virtual void getSensitiveStrikes(const CVolRequest* volRequest,
                                     OutputNameConstSP outputName,
                                     const SensitiveStrikeDescriptor& sensStrikeDesc,
                                     DoubleArraySP sensitiveStrikes) const;

    /** The redirects everything to the SimpleEquity contained within */
    virtual bool accept(ICollector* collector) const;

    /** Create a PseudoSimpleEquity. asset must be a SimpleEquity. */
    static PseudoSimpleEquity* create(const CAsset*        asset,
                                      const DateTimeArray& divCritDates,
                                      bool                 isCall,
                                      int                  noExerciseWindow);

    /** return a pdf calculator */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;

private:
    PseudoSimpleEquity();
    PseudoSimpleEquity(const PseudoSimpleEquity& rhs);
    PseudoSimpleEquity& operator=(const PseudoSimpleEquity& rhs);

    PseudoSimpleEquity(const EquityBase*  simpleEquity,
                       const DateTimeArray& divCritDates,
                       bool                 isCall,
                       int                  noExerciseWindow);

    // fields
    DateTime        horizonDate;
    DateTimeArraySP divCritDates;
    EquityBaseConstSP parent;   // ref to the equity from which 
                                  // this PseudoSimpleEquity is "derived"
    bool            isCall;
    int             noExerciseWindow;
};

typedef smartPtr<PseudoSimpleEquity> PseudoSimpleEquitySP;
typedef smartConstPtr<PseudoSimpleEquity> PseudoSimpleEquityConstSP;

// support for wrapper class
typedef MarketWrapper<PseudoSimpleEquity> PseudoSimpleEquityWrapper;

DRLIB_END_NAMESPACE
#endif
