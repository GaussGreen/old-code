//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyFutureSwapSchedule.hpp
//
//   Description : Energy Swap Schedule
//
//   Author      : Sean Chen
//
//   Date        : Aug. 29, 2005
//
//----------------------------------------------------------------------------

#ifndef _EnergyFutureSwapSchedule_H_
#define _EnergyFutureSwapSchedule_H_

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/Object.hpp"
#include "edginc/EnergyUnderlyer.hpp"
#include "edginc/EnergyDWMYPeriod.hpp"
#include "edginc/EnergyTermPeriod.hpp"
#include "edginc/EnergySwapScheduleBase.hpp"


DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL EnergyFutureSwapSchedule : public EnergySwapScheduleBase
{

public:

    static CClassConstSP const TYPE;

    friend class EnergyFutureSwapScheduleHelper;

    virtual ~EnergyFutureSwapSchedule();

    virtual void validatePop2Object();

    //virtual void getMarket(const IModel* model, const CMarketDataSP  market);


private:

    EnergyFutureSwapSchedule();

    string scheduleType;       //Payment schedule following INDEX or OPTION rule

};

typedef MarketWrapper<EnergyFutureSwapSchedule> EnergyFutureSwapScheduleWrapper;
typedef smartPtr<EnergyFutureSwapSchedule> EnergyFutureSwapScheduleSP;
typedef smartConstPtr<EnergyFutureSwapSchedule> EnergyFutureSwapScheduleConstSP;


DRLIB_END_NAMESPACE

#endif
