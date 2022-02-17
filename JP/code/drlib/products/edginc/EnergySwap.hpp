//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergySwap.hpp
//
//   Description : Energy Swap Schedule
//
//   Author      : Sean Chen
//
//   Date        : Aug. 29, 2005
//
//----------------------------------------------------------------------------

#ifndef _EnergySwap_H_
#define _EnergySwap_H_

#include "edginc/config.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/FirstSensDate.hpp"

#include "edginc/ClosedFormEnergy.hpp"
#include "edginc/EnergyFuturesCurve.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/EnergyStreamSchedule.hpp"
#include "edginc/EnergySwapScheduleBase.hpp"


DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL EnergySwap : public CInstrument, 
                   public virtual ClosedFormEnergy::IIntoProduct,
                   public LastSensDate,
                   public FirstSensDate

{

public:

    static CClassConstSP const TYPE;

    friend class EnergySwapHelper;

    virtual ~EnergySwap();

    virtual void validatePop2Object();

    virtual void GetMarket(const IModel* model, const CMarketDataSP  market);

    EnergySwapScheduleBaseSP getSwapSchedule() const;
    string getNotionalType() const;
    EnergyFuturesCurveConstSP getFuturesCurve() const;
    YieldCurveConstSP getYieldCurve() const;

    DateTime getValueDate() const;

    double price();

	double getStrike() const;
	double getNotional() const;

    ClosedFormEnergy::IProduct* createProduct(ClosedFormEnergy* model) const;

    DateTime endDate(const Sensitivity* sensControl) const;
    DateTime beginDate(const SensControl* sensControl) const; 
    string getCcy() const;
    void Validate();
    string discountYieldCurveName() const;

    void addMoreToResults(Control* control, Results* results);

private:

    EnergySwap();

    EnergySwapScheduleBaseSP swapScheduleSP;
    EnergyFuturesCurveWrapper futuresCurve;
    YieldCurveWrapper yieldCurve;
    double rate;
    double notionalAmount;
    string notionalType;
    string dealDirection;
    double pastAverage;           // on coupons > today
    bool pastAverageInclToday;    // Includes Today's fixing?

    ObjectArraySP outputFixed; // $unregistered
    ObjectArraySP outputFloating; // $unregistered
    
};

typedef MarketWrapper<EnergySwap> EnergySwapWrapper;
typedef smartPtr<EnergySwap> EnergySwapSP;
typedef smartConstPtr<EnergySwap> EnergySwapConstSP;


DRLIB_END_NAMESPACE

#endif
