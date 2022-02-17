//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergySwapScheduleBase.hpp
//
//   Description : Energy Swap Schedule
//
//   Author      : Sean Chen
//
//   Date        : Aug. 29, 2005
//
//----------------------------------------------------------------------------

#ifndef _EnergySwapScheduleBase_H_
#define _EnergySwapScheduleBase_H_

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/Object.hpp"
#include "edginc/EnergyUnderlyer.hpp"
#include "edginc/EnergyDWMYPeriod.hpp"
#include "edginc/EnergyTermPeriod.hpp"
#include "edginc/EnergyStreamSchedule.hpp"


DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL EnergySwapScheduleBase : public CObject
{

public:

    static CClassConstSP const TYPE;

    friend class EnergySwapScheduleBaseHelper;
    friend class GetEnergySwapScheduleDetailsAddin;

    virtual ~EnergySwapScheduleBase();

    virtual void validatePop2Object();

  //  virtual void getMarket(const IModel* model, const CMarketDataSP  market);

    DateTime getValueDate() const;
    DateTime getStartDate() const;
    DateTime getEndDate() const;

    double calculatePV(const EnergyFuturesCurveSP& futuresCurve,
                       const YieldCurveSP& yieldCurve,
                       double rate,
                       double notionalAmount,
                       const string& notionalType,
                       const string& dealDerection,
                       double pastAverage,
                       bool pastAverageInclToday);  

    ObjectArraySP getFixedLegDetails() const;
    ObjectArraySP getFloatingLegDetails() const;
    EnergyUnderlyerConstSP getEnergyUnderlyer() const;

	EnergyStreamScheduleFloatingConstSP getFloatingEnergyStreamSchedule() const;
	EnergyStreamScheduleFixedConstSP getFixedEnergyStreamSchedule() const;

protected:

    EnergySwapScheduleBase(CClassConstSP clazz);
    EnergySwapScheduleBase();

    DateTime valueDate;                         // Value Date
    EnergyUnderlyerWrapper energyUnderlyer;     // Underlyer Wrapper
    EnergyTermPeriodSP dealPeriod;                // Handle to EnergyDealPeriod

    EnergyDWMYPeriod payFrequencyFixed;       //Fixed Payment Frequency, e.g. nD,M,W,Y $unregistered
    EnergyDWMYPeriod payFrequencyFloating;    //Floating Payment Frequency $unregistered


    // Averaging stuff...
    EnergyDWMYPeriod avgPeriod;     // Averaging Period nD,M,W,Y. Default, whole coupon period
    EnergyDWMYPeriod avgFrequency;  // Averaging Frequency
    int avgDaysB4End;   // Average Period Starts n days before END. 0(default), 1,2,...,

    int settleDays;               // Days to Settlement
    int nearbyRel;           // Relative Nearby Contract
    string nearbyLabel;
    EnergyContractLabel nearbyAbsLabel;           // Absolute Nearby Contract, eg, Sep08.

    DateTime startDate;           // Deal Start Date
    DateTime endDate;             // Deal End Date

    EnergyStreamScheduleFixedSP fixedEnergyStreamScheduleSP;
    EnergyStreamScheduleFloatingSP floatingEnergyStreamScheduleSP;


};

typedef MarketWrapper<EnergySwapScheduleBase> EnergySwapScheduleBaseWrapper;
typedef smartPtr<EnergySwapScheduleBase> EnergySwapScheduleBaseSP;
typedef smartConstPtr<EnergySwapScheduleBase> EnergySwapScheduleBaseConstSP;


DRLIB_END_NAMESPACE

#endif
