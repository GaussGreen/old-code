//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : MQQuasiIRVol.hpp
//
//   Description : Base Class for IR Vol for MultiQ quasi pricing
//
//   Author      : Keith Jia
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------

#ifndef _MQQUASIIRVOL_HPP
#define _MQQUASIIRVOL_HPP

#include "edginc/IRVolBase.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IDistribution1D);
FORWARD_DECLARE(MultiQDistribution);


class MARKET_DLL MQQuasiIRVol: public IRVolCommon
{
public:
    static CClassConstSP const TYPE;
    virtual ~MQQuasiIRVol();

    virtual void validatePop2Object();

    //-----------------------------------------------------
    // MarketObject methods
    //-----------------------------------------------------
    virtual string getName() const;

    /** set valueDate from market */
    virtual void getMarket(const IModel *model, 
                           const MarketData *market);

public:
    virtual TimeMetricConstSP GetTimeMetric() const;
    virtual HolidayConstSP getHoliday() const;

protected:
    /** create MultiQDistribution mq */
    virtual MultiQDistribution* CreateMQ(ExpiryConstSP tenor,
                                         double parFwdSprd, 
                                         const DateTime& resetDate) const;

protected:
    MQQuasiIRVol(CClassConstSP clazz = TYPE);

    string           name;
    TimeMetricSP     metric;

    //market data inputs
    DateTime            valueDate;
    int                 spotOffset;
    ExpiryArraySP       atmTenors;
    DateTimeArraySP     atmMatDates;
    CDoubleMatrixSP     atmVols;           //at the money vols, at the above maturiries and tenors.

    ExpiryArraySP       otmTenors;         //otm Tenors (z-axis)
    DateTimeArraySP     otmMatDates;       //maturities (y-axis).
    CDoubleArraySP      strikeBoundaries;  //strike ratios or deltas (x-axis)
    DoubleMatrixArraySP impliedVols;       //vols or Q's cube 
    bool                isInputQs;         //true if deltas and Qs are inputs.

    HolidayWrapper  hols;
    string          swapDCC;
    string          swapBDC;

    //transient fields

    /** maturities in terms of DateTime for interpolation, still market data.
     ** reset date comes from instrument, and tenorDate is determined by reset date.
     */
    mutable ExpiryConstSP        tenor;             //underlying tenor.
    mutable DateTime             spotDate;
    DayCountConventionSP dcc;
    BadDayConventionSP   bdc;

    /** info from instrument for specific mq creation */
    mutable double atmBSVol;
    mutable DateTime resetDate;

    /** MQ object */
    //mutable IDistribution1DSP mq;    


private:
    static void load(CClassSP& clazz);
};


DECLARE(MQQuasiIRVol);
typedef MarketWrapper<MQQuasiIRVol> MQQuasiIRVolWrapper;


DRLIB_END_NAMESPACE

#endif
