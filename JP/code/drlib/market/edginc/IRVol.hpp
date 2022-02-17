//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRVol.hpp
//
//   Description : Interest rate vol
//
//   Author      : Andrew J Swain
//
//   Date        : 7 November 2001
//
//
//----------------------------------------------------------------------------

#ifndef _IRVOL_HPP
#define _IRVOL_HPP

#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/VolatilityBS.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/IRCalib.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/Theta.hpp"
#include "edginc/IRVegaParallel.hpp"
#include "edginc/IRVegaPointwise.hpp"
#include "edginc/DayCountConvention.hpp"

DRLIB_BEGIN_NAMESPACE
class CriticalDateCollector;

class MARKET_DLL IRVol: public IRVolCommon,
             public virtual IVolatilityBS,
             public IRVegaParallel::RestorableShift,
             public IRVegaPointwise::IShift,
             public Theta::IShift {
public:
	friend class IrConverter;
	
    /** Enum VolType */
    enum VolType
    {
        BASE_VOL,   // caplet vols (tenors <= 12M)
        SWAP_VOL    // swaption vol surface (tenors >= 1Y)
    };

    /** Enum CalibType */
    enum CalibType
    {
        CAP,        // caplet
        CMS,        // constant maturity swap
        FIX,        // fixed maturity swap
        FMD         // fixed swap maturity date
    };

    static CClassConstSP const TYPE;

    /** Constructor - in general though object is created via data dict.
        Note that the metric must already be populated with its market
        data (ie its holidays) */
    IRVol(const string&           volName,
          const TimeMetric*       timeMetric,
          const DoubleMatrix&     vol,
          const ExpiryArray&      expiryArray,
          const ExpiryArray&      tenorArray,
          const DateTime&         baseDateTime,
          const HolidayWrapper    holidayWrapper,
          const BadDayConvention& badDayConv);
    
    /** Returns name of vol */
    virtual string getName() const;

    /** Combines market and instrument data together to give a
        Processed Vol */
    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const MarketObject*  yc) const;

    /** overrides default */
    virtual void validatePop2Object();

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift);

    /** Get methods. */
    virtual DateTime getBaseDate() const { return baseDate; }
    virtual const ExpiryArray&  getExpiries() const { return (*expiry); };
    virtual const ExpiryArray&  getTenors() const { return (*tenor); }
    virtual const DoubleMatrix& getVolMatrix() const { return matrix; }

    virtual ~IRVol();

    virtual string sensName(IRVegaParallel* shift) const;
    virtual bool sensShift(IRVegaParallel* shift);
    virtual void sensRestore(IRVegaParallel* shift);

    virtual string sensName(IRVegaPointwise* shift) const;
    virtual bool sensShift(IRVegaPointwise* shift);


protected:
    static void acceptValueDateCollector(const IRVol* vol,
                                         CValueDateCollector* collector);

    friend class IRVolHelper;
    friend class IRInterpVol;
    IRVol();
    IRVol(const IRVol& rhs);
    IRVol& operator=(const IRVol& rhs);

    void buildCache(bool removeHistoricExpiries);
    DateTime expiryToDate(const Expiry* expiryBM) const;

    DateTime calcSwapStart(int expiryIdx) const;
    DateTime calcSwapMaturity(const DateTime& swapStart,
                              int             tenorIdx) const;
    DateTime calcSwapMaturity(int  expiryIdx,
                              int  tenorIdx) const;

    IRGridPointAbsSP createIRGridPointAbs(int  expiryIdx,
                                          int  tenorIdx) const;

    /************* exported fields *************/
    string                name;    // name of the vol
    TimeMetricSP          metric;
    ExpiryArraySP         expiry;  // option maturity
    int                   spotOffset; // swap starts at option mat + spotOffset
    ExpiryArraySP         tenor;   // underlying maturity
    DoubleMatrix          matrix;       
    IRCalibWrapper        params;
    DateTime              baseDate;
    HolidayWrapper        hols;
    BadDayConventionSP    bdc;
    DayCountConventionSP  swapDayCount;
    MaturityPeriodSP      swapFrequency;

    /************* transient fields *************/
    DoubleArray   tradYears;// trading time, in years, between tenors
    DateTimeArray expiryDates; // cached - internally derived from expiry
    DateTimeArray tenorDates; // cached - internally derived from tenor

    /************* end of fields *************/
};

typedef smartConstPtr<IRVol> IRVolConstSP;
typedef smartPtr<IRVol> IRVolSP;
typedef MarketWrapper<IRVol> IRVolWrapper;

class MARKET_DLL IRVolSpot : public IRVolCommon {
public:
    static CClassConstSP const TYPE;

    /** Returns name of vol */
    virtual string getName() const;

    /** Combines market and instrument data together to give a
    Processed Vol */
    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
        const MarketObject*  yc) const;

    /** overrides default */
    virtual void validatePop2Object();

    virtual void getMarket(const IModel* model, const MarketData* market);

    virtual ~IRVolSpot();

    class MARKET_DLL IRVolSpotProcessed : public virtual CVolProcessed,
        public CObject
    {

    private:
        const IRVolSpot* vol;
    public:
        static CClassConstSP const TYPE;
        const DoubleArray& getSpotVols() const { return vol->spotVols; }
        const DateTimeArray& getSpotExpiryDates() const { return vol->spotExpiryDates; }
        IRVolSpotProcessed(const IRVolSpot* _vol) : CObject(TYPE), vol(_vol) {}
        IRVolSpotProcessed() : CObject(TYPE), vol(0) {} // for vol request

        virtual string getName() const { return vol->name; }

        virtual double calcTradingTime(const DateTime &date1, 
            const DateTime &date2) const;
        virtual TimeMetricConstSP GetTimeMetric()const;
        static void load(CClassSP& clazz);

    };
    friend class IRVolSpotProcessed;

protected:
    static void acceptValueDateCollector(const IRVolSpot* vol,
        CValueDateCollector* collector);
    static void acceptCriticalDateCollector(
        const IRVolSpot*           vol,
        CriticalDateCollector* collector) ;

    IRVolSpot();
    IRVolSpot(const IRVolSpot& rhs);
    IRVolSpot& operator=(const IRVolSpot& rhs);

    static IObject* defaultIRVolSpot();
    static void load(CClassSP& clazz);


    /************* exported fields *************/
    string                name;    // name of the vol
    DateTime              baseDate;

    DoubleArray           spotVols;
    DateTimeArray         spotExpiryDates; // cached - internally derived from expiry

    IRCalibWrapper        params;

    /************* end of fields *************/

};

DRLIB_END_NAMESPACE

#endif
