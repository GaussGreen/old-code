//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyTermPeriodDates.hpp
//
//   Description : Defines a deal period with start/end dates for energy. 
//
//   Author      : Sean Chen
//
//   Date        : 26 Sept. 2005
//
//----------------------------------------------------------------------------

#ifndef _EnergyTermPeriodDates_H_
#define _EnergyTermPeriodDates_H_

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/EnergyTermPeriod.hpp"

#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE


class TOOLKIT_DLL EnergyTermPeriodDates : public EnergyTermPeriod
{
public:

	static CClassConstSP const TYPE;

    friend class EnergyTermPeriodDatesHelper;

    //EnergyTermPeriodDates(const EnergyTermPeriodDates& thePeriod);

    void validatePop2Object();

    virtual ~EnergyTermPeriodDates();

private:

    EnergyTermPeriodDates();

    DateTime dateStart;
    DateTime dateEnd;

};

typedef smartPtr<EnergyTermPeriodDates> EnergyTermPeriodDatesSP;
typedef smartConstPtr<EnergyTermPeriodDates> EnergyTermPeriodDatesConstSP;

DRLIB_END_NAMESPACE

#endif
