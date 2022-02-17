#ifndef EDR_SRMFXVOLBASE_HPP
#define EDR_SRMFXVOLBASE_HPP

#include "edginc/FXVolBase.hpp"
//#include "edginc/VolProcessed.hpp"
//#include "edginc/Calibrator.hpp"
//#include "edginc/IDynamicsParameter.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/** Not clear whether this vol will be reused by the equity components. So
    for now just put here directly. It is market data so you could argue it
    ought to be in the market directory ... */
class MARKET_DLL SRMFXVolBase : public FXVolBase
{
public:
    /** The type of original FX vol - in particular models (via the market
        data fetcher) will need to select a particular SRM type of fx vol */
    static CClassConstSP const VOL_TYPE;

    //static CClassConstSP const TYPE; // inside SRMFX class for ease
    virtual ~SRMFXVolBase(){}

    /** force all times of day to be the same as today - unfortunate hack */
    virtual void forceTimeOfDay(void) = 0;

    /** returns merged benchmark dates (spot vol and comp vol) */
    virtual DateTimeArray getVolBmDates() const =0;

    // Accessors

    virtual const DateTimeArray& getSpotVolDate() const =0;
    virtual const DoubleArray& getSpotVol() const =0;
    virtual const DoubleArray& getSmileA1() const =0;
    virtual const DoubleArray& getSmileA2() const =0;
    virtual const DoubleArray& getSmileA3() const =0;
    virtual const DateTimeArray& getSmileDate() const =0;

    virtual void setToday( const DateTime& date ) =0;

protected:
    SRMFXVolBase(CClassConstSP type): FXVolBase(type) {}

private:
    static void load(CClassSP& clazz);

};

DECLARE(SRMFXVolBase);

DRLIB_END_NAMESPACE
#endif


