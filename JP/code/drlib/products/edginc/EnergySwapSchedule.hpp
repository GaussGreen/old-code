//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergySwapSchedule.hpp
//
//   Description : Energy Swap Schedule
//
//   Author      : Sean Chen
//
//   Date        : Aug. 29, 2005
//
//----------------------------------------------------------------------------

#ifndef _ENERGYSWAPSCHEDULE_H_
#define _ENERGYSWAPSCHEDULE_H_

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/Object.hpp"
#include "edginc/EnergyUnderlyer.hpp"
#include "edginc/EnergyDWMYPeriod.hpp"
#include "edginc/EnergyTermPeriod.hpp"
#include "edginc/EnergySwapScheduleBase.hpp"
#include "edginc/EnergyStreamSchedule.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL EnergySwapSchedule : public EnergySwapScheduleBase
{

public:

    static CClassConstSP const TYPE;

    friend class EnergySwapScheduleHelper;

    virtual ~EnergySwapSchedule();

    virtual void validatePop2Object();


private:

    EnergySwapSchedule();
 
    int rollDay;                  // Roll Day Overriding


};

typedef MarketWrapper<EnergySwapSchedule> EnergySwapScheduleWrapper;
typedef smartPtr<EnergySwapSchedule> EnergySwapScheduleSP;
typedef smartConstPtr<EnergySwapSchedule> EnergySwapScheduleConstSP;


DRLIB_END_NAMESPACE

#endif
